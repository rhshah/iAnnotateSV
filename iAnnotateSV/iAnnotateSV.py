"""
Created on 25/11/2014.

@author: Ronak H Shah

"""

import polars as pl
import typer
import time
import helper as hp
import AnnotateEachBreakpoint as aeb
import PredictFunction as pf
import FindCanonicalTranscript as fct
import AddExternalAnnotations as aea
import AnnotationForKinaseDomain as kda
import VisualizeSV as vsv
from models import *
import os
import sys
import logging
from rich import print
from rich.logging import RichHandler

app = typer.Typer()

FORMAT = "%(message)s"
logging.basicConfig(
    level="INFO", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)

log = logging.getLogger("rich")


@app.command()
def main(
    ref_file_version: str = typer.Option(
        ...,
        "--refFileVersion",
        "-r",
        help="Which human reference file to be used, hg18, hg19 or hg38",
    ),
    ref_file: str = typer.Option(
        None,
        "--refFile",
        "-rf",
        help="Human reference file location to be used",
    ),
    output_file_prefix: str = typer.Option(
        ...,
        "--outputFilePrefix",
        "-ofp",
        help="Prefix for the output file",
    ),
    output_dir: str = typer.Option(
        ...,
        "--outputDir",
        "-o",
        help="Full Path to the output dir",
    ),
    sv_file: str = typer.Option(
        ...,
        "--svFile",
        "-i",
        help="Location of the structural variants file to annotate",
    ),
    distance: int = typer.Option(
        3000,
        "--distance",
        "-d",
        help="Distance used to extend the promoter region",
    ),
    auto_select: bool = typer.Option(
        True,
        "--autoSelect",
        "-a",
        help="Auto Select which transcript to be used[default]",
    ),
    canonical_transcripts: str = typer.Option(
        None,
        "--canonicalTranscripts",
        "-c",
        help="Location of canonical transcript list for each gene. Use only if you want the output for specific transcripts for each gene.",
    ),
    plot_sv: bool = typer.Option(
        False,
        "--plotSV",
        "-p",
        help="Plot the structural variant in question",
    ),
    uniprot_file: str = typer.Option(
        None,
        "--uniprotFile",
        "-u",
        help="Location of UniProt list contain information for protein domains. Use only if you want to plot the structural variant",
    ),
    repeat_file: str = typer.Option(
        None,
        "--repeatFile",
        "-rr",
        help="Location of the Repeat Region Bed File",
    ),
    dgv_file: str = typer.Option(
        None,
        "--dgvFile",
        "-dgv",
        help="Location of the Database of Genomic Variants Bed File",
    ),
    cosmic_consensus_file: str = typer.Option(
        None,
        "--cosmicConsensusFile",
        "-cc",
        help="Location of the Cosmic Consensus TSV file",
    ),
    cosmic_counts_file: str = typer.Option(
        None,
        "--cosmicCountsFile",
        "-cct",
        help="Location of the Cosmic Counts TSV file",
    ),
    verbose: bool = typer.Option(
        True,
        "--verbose",
        "-v",
        help="make lots of noise [default]",
    ),
):
    """
    Annotate SV based on a specific human reference
    """

    start_time = time.time()

    # Get current location
    this_dir, this_filename = os.path.split(__file__)

    # Check if file for canonical transcript is given or not
    if canonical_transcripts:
        auto_select = False

    if ref_file_version in {'hg18', 'hg19', 'hg38'}:
        if not ref_file:
            ref_file = f"{ref_file_version}.sv.table.txt"
            ref_file = os.path.join(this_dir, "data/references", ref_file)

        if repeat_file:
            rr_path = repeat_file
        else:
            rr_filename = f"{ref_file_version}_repeatRegion.tsv"
            rr_path = os.path.join(this_dir, "data/repeat_region", rr_filename)
            repeat_file = rr_path

        if dgv_file:
            dgv_path = dgv_file
        else:
            dgv_filename = f"{ref_file_version}_DGv_Annotation.tsv"
            dgv_path = os.path.join(
                this_dir, "data/database_of_genomic_variants", dgv_filename)
            dgv_file = dgv_path

        if cosmic_consensus_file:
            cc_path = cosmic_consensus_file
        else:
            cc_filename = "cancer_gene_census.tsv"
            cc_path = os.path.join(this_dir, "data/cosmic", cc_filename)
            cosmic_consensus_file = cc_path

        if cosmic_counts_file:
            cct_path = cosmic_counts_file
        else:
            cct_filename = "cosmic_fusion_counts.tsv"
            cct_path = os.path.join(this_dir, "data/cosmic", cct_filename)
            cosmic_counts_file = cct_path

        if not uniprot_file:
            up_filename = f"{ref_file_version}.uniprot.spAnnot.table.txt"
            uniprot_file = str(os.path.join(
                this_dir, "data/UcscUniprotdomainInfo", up_filename))
        uniprot_path = uniprot_file
        all_canonical_transcripts_path = str(os.path.join(
            this_dir, "data/canonicalInfo/canonical_transcripts.txt"))
    else:
        log.fatal(
            "iAnnotateSV: Please enter correct reference file version. Values can be: hg18 or hg19 or hg38")
        sys.exit(1)

    log.info(f"Reading reference file: {ref_file}")
    refDF = hp.ReadFile(ref_file)
    NewRefDF = hp.ExtendPromoterRegion(refDF, distance)

    log.info(f"Reading SV file: {sv_file}")
    svDF = hp.ReadFile(sv_file)

    log.info("Processing SVs...")
    annDF = processSV(svDF, NewRefDF, auto_select, canonical_transcripts, all_canonical_transcripts_path, uniprot_path, ref_file, verbose)

    plotDF = annDF.clone()

    # Print to TSV file
    out_file_prefix_path = os.path.join(output_dir, f"{output_file_prefix}_functional.txt")
    log.info(f"Writing functional annotations to: {out_file_prefix_path}")
    annDF.write_csv(out_file_prefix_path, separator='\t')

    # Add External Annotations
    log.info("Adding External Annotations...")
    make_command_line_for_aea = (
        f"-r {rr_path} -d {dgv_path} -c {cosmic_consensus_file} -cct {cosmic_counts_file} "
        f"-s {out_file_prefix_path} -ofp {output_file_prefix}_Annotated -o {output_dir}"
    )
    try:
        aea.main(make_command_line_for_aea)
    except Exception as e:
        log.error(f"Error in AddExternalAnnotations: {e}")

    # Plot if required
    if plot_sv:
        log.info("Plotting Each Structural Variants")
        plotSV(plotDF, NewRefDF, uniprot_path, output_dir, output_file_prefix, verbose)

    end_time = time.time()
    log.info(f"Finished Running the Annotation Process!!! Elapsed time: {end_time - start_time:.2f} seconds")


def processSV(svDF, refDF, auto_select, canonical_transcripts, all_canonical_transcripts_path, uniprot_path, ref_file, verbose):
    log.info("Processing Each Structural Variants...")

    # Read Canonical Transcript if the file is given in the cmdline
    if canonical_transcripts:
        ctDict = hp.ReadTranscriptFile(canonical_transcripts)
        log.info(f"Using canonical transcripts from: {canonical_transcripts}")
    else:
        ctDict = None
        log.info("Not using canonical transcripts.")

    def annotate_row(row):
        chr1, chr2, pos1, pos2, str1, str2 = row['chr1'], row['chr2'], row['pos1'], row['pos2'], row['str1'], row['str2']
        log.debug(f"Annotating SV: {chr1}:{pos1}:{str1} - {chr2}:{pos2}:{str2}")

        b1, b2 = (None,) * 2
        gene1, transcript1, site1, gene2, transcript2, site2, fusionFunction = ("-",) * 7  # Default values

        try:
            if auto_select:
                log.debug(f"Auto-selecting transcripts.")
                (gene1, transcript1, site1, zone1, strand1, intronnum1,
                 intronframe1) = aeb.AnnotateEachBreakpoint(chr1, pos1, str1, refDF, auto_select)
                (gene2, transcript2, site2, zone2, strand2, intronnum2,
                 intronframe2) = aeb.AnnotateEachBreakpoint(chr2, pos2, str2, refDF, auto_select)
                ann1S = pl.Series([gene1, transcript1, site1, zone1, strand1, str1, intronnum1, intronframe1],
                                  name='ann1')
                ann2S = pl.Series([gene2, transcript2, site2, zone2, strand2, str2, intronnum2, intronframe2],
                                  name='ann2')
                fusionFunction = pf.PredictFunctionForSV(ann1S, ann2S)
                log.debug(f"Fusion function: {fusionFunction}")
            else:
                log.debug(f"Using canonical transcripts.")
                try:
                    (gene1List, transcript1List, site1List, zone1List, strand1List, intronnum1List,
                     intronframe1List) = aeb.AnnotateEachBreakpoint(chr1, pos1, str1, refDF, auto_select)
                    (gene1, transcript1, site1, zone1, strand1, intronnum1, intronframe1) = fct.FindCT(
                        gene1List, transcript1List, site1List, zone1List, strand1List, intronnum1List, intronframe1List, ctDict)
                except (IntergenicError, ChrError) as b1:
                    log.info(f"iAnnotateSV: {str(b1)}")
                    (gene1, transcript1, site1, zone1, strand1, intronnum1, intronframe1) = ("-",) * 7
                try:
                    (gene2List, transcript2List, site2List, zone2List, strand2List, intronnum2List,
                     intronframe2List) = aeb.AnnotateEachBreakpoint(chr2, pos2, str2, refDF, auto_select)
                    (gene2, transcript2, site2, zone2, strand2, intronnum2, intronframe2) = fct.FindCT(
                        gene2List, transcript2List, site2List, zone2List, strand2List, intronnum2List, intronframe2List, ctDict)
                except (IntergenicError, ChrError) as b2:
                    log.info(f"iAnnotateSV: {str(b2)}")
                    (gene2, transcript2, site2, zone2, strand2, intronnum2, intronframe2) = ("-",) * 7
                ann1S = pl.Series([gene1, transcript1, site1, zone1, strand1, str1, intronnum1, intronframe1],
                                  name='ann1')
                ann2S = pl.Series([gene2, transcript2, site2, zone2, strand2, str2, intronnum2, intronframe2],
                                  name='ann2')
                if not any([b1, b2]):
                    fusionFunction = pf.PredictFunctionForSV(ann1S, ann2S)
                else:
                    fusionFunction = "-"
                log.debug(f"Fusion function: {fusionFunction}")
        except Exception as e:
            log.error(f"Error processing row: {e}")
            # Return default values in case of an error
            return [chr1, pos1, str1, chr2, pos2, str2, gene1, transcript1, site1, gene2, transcript2, site2, fusionFunction]

        return [chr1, pos1, str1, chr2, pos2, str2, gene1, transcript1, site1, gene2, transcript2, site2, fusionFunction]

    # Annotate each row
    results = []
    for row in svDF.iter_rows(named=True):
        results.append(annotate_row(row))

    # Convert results back to a Polars DataFrame
    annDF = pl.DataFrame(results,
                         schema=[
                             ("chr1", pl.Utf8), ("pos1", pl.Int64), ("str1", pl.Utf8),
                             ("chr2", pl.Utf8), ("pos2", pl.Int64), ("str2", pl.Utf8),
                             ("gene1", pl.Utf8), ("transcript1", pl.Utf8), ("site1", pl.Utf8),
                             ("gene2", pl.Utf8), ("transcript2", pl.Utf8), ("site2", pl.Utf8),
                             ("fusion", pl.Utf8)
                         ],
                         )
    log.info("Finished processing each structural variants")

    if canonical_transcripts:
        annDF = kda.run(annDF, ref_file, canonical_transcripts,
                        all_canonical_transcripts_path, uniprot_path, verbose)
        return annDF
    else:
        return annDF


def plotSV(svDF, refDF, uniprot_path, output_dir, output_file_prefix, verbose):
    log.info("Will now try to plot Each Structural Variants")
    upDF = None
    if(os.path.isfile(uniprot_path)):
        upDF = hp.ReadFile(uniprot_path)
    else:
        log.fatal(
            "iAnnotateSV: %s file does not exist!!, Please use it to plot structural variants",
            uniprot_path)
        sys.exit(1)

    vsv.VisualizeSV(svDF, refDF, upDF, output_dir, output_file_prefix, verbose)


if __name__ == "__main__":
    app()