"""
Created on 12/23/2015.

@Ronak Shah

"""

import sys
import time
import logging
import AnnotateForRepeatRegion as afr
import AnnotateForCosmic as afc
import AnnotateForDGv as afd
import polars as pl
import typer
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
    repeat_file: str = typer.Option(
        ...,
        "--repeatFile",
        "-r",
        help="Location of the Repeat Region Bed File",
    ),
    dgv_file: str = typer.Option(
        ...,
        "--dgvFile",
        "-d",
        help="Location of the Database of Genomic Variants Bed File",
    ),
    cosmic_consensus_file: str = typer.Option(
        ...,
        "--cosmicConsensusFile",
        "-c",
        help="Location of the Cosmic Consensus TSV file",
    ),
    cosmic_counts_file: str = typer.Option(
        ...,
        "--cosmicCountsFile",
        "-cct",
        help="Location of the Cosmic Counts TSV file",
    ),
    sv_file: str = typer.Option(
        ...,
        "--svFile",
        "-s",
        help="Location of the structural variant file to be annotated",
    ),
    output_file_prefix: str = typer.Option(
        ...,
        "--outputFilePrefix",
        "-ofp",
        help="Full path with prefix name for the output file",
    ),
    output_dir: str = typer.Option(
        ...,
        "--outputDir",
        "-o",
        help="Full Path to the output dir",
    ),
    verbose: bool = typer.Option(
        True,
        "--verbose",
        "-v",
        help="make lots of noise [default]",
    ),
):
    """
    Add External Annotations to the Structural Variants.
    """
    start_time = time.time()

    out_file_txt = f"{output_dir}/{output_file_prefix}.txt"
    out_file_exl = f"{output_dir}/{output_file_prefix}.xlsx"
    out_file_json = f"{output_dir}/{output_file_prefix}.json"

    if verbose:
        log.info(f"Reading {sv_file}...")
    try:
        data = ReadSVFile(sv_file, verbose)
        if verbose:
            log.info(f"Finished Reading {sv_file}")
    except Exception as e:
        log.error(f"Error reading SV file: {e}")
        sys.exit(1)

    if verbose:
        log.info(f"Reading {repeat_file}...")
    try:
        repeat_region_dict = afr.ReadRepeatFile(repeat_file, verbose)
        if verbose:
            log.info(f"Finished Reading {repeat_file}")
    except Exception as e:
        log.error(f"Error reading repeat region file: {e}")
        sys.exit(1)

    if verbose:
        log.info(f"Reading {dgv_file}...")
    try:
        dgv_dict = afd.ReadDGvFile(dgv_file, verbose)
        if verbose:
            log.info(f"Finished Reading {dgv_file}")
    except Exception as e:
        log.error(f"Error reading DGV file: {e}")
        sys.exit(1)

    data = data.with_columns([
        pl.lit("-").alias("Cosmic_Fusion_Counts"),
        pl.lit("-").alias("repName-repClass-repFamily:-site1"),
        pl.lit("-").alias("repName-repClass-repFamily:-site2"),
        pl.lit("-").alias("CC_Chr_Band"),
        pl.lit("-").alias("CC_Tumour_Types(Somatic)"),
        pl.lit("-").alias("CC_Cancer_Syndrome"),
        pl.lit("-").alias("CC_Mutation_Type"),
        pl.lit("-").alias("CC_Translocation_Partner"),
        pl.lit("-").alias("DGv_Name-DGv_VarType-site1"),
        pl.lit("-").alias("DGv_Name-DGv_VarType-site2")
    ])

    for count, row in enumerate(data.iter_rows(named=True)):
        sv_chr1 = row['chr1']
        sv_pos1 = row['pos1']
        sv_chr2 = row['chr2']
        sv_pos2 = row['pos2']
        sv_gene1 = row['gene1']
        sv_gene2 = row['gene2']
        if verbose:
            log.info(
                "Processing Record: %s\t%s\t%s\t%s\t%s\t%s",
                sv_chr1,
                sv_pos1,
                sv_chr2,
                sv_pos2,
                sv_gene1,
                sv_gene2,
            )

        # Repeat Region Data
        try:
            rr_loc1, rr_loc2 = afr.AnnotateRepeatRegion(
                verbose, count, row, repeat_region_dict)
        except Exception as e:
            log.error(f"Error annotating repeat region: {e}")
            rr_loc1, rr_loc2 = [], []

        data = data.with_columns([
            pl.when(pl.lit(True))
            .then(pl.lit("<=>".join(rr_loc1)))
            .otherwise(pl.col("repName-repClass-repFamily:-site1"))
            .alias("repName-repClass-repFamily:-site1")
        ])
        data = data.with_columns([
            pl.when(pl.lit(True))
            .then(pl.lit("<=>".join(rr_loc2)))
            .otherwise(pl.col("repName-repClass-repFamily:-site2"))
            .alias("repName-repClass-repFamily:-site2")
        ])

        # Cosmic Consensus Data
        try:
            cc_SV = afc.AnnotateFromCosmicCensusFile(
                cosmic_consensus_file, verbose, count, row)
            cct_SV = afc.AnnotateFromCosmicFusionCountsFile(
                cosmic_counts_file, verbose, count, row)
        except Exception as e:
            log.error(f"Error annotating cosmic data: {e}")
            cc_SV, cct_SV = [], None

        ccA, ccB, ccC, ccD, ccE = ([] for _ in range(5))
        if cc_SV:
            for cc in cc_SV:
                ccData = cc.split('\t')
                ccA.append(ccData[0])
                ccB.append(ccData[1])
                ccC.append(ccData[2])
                ccD.append(ccData[3])
                ccE.append(ccData[4])

        data = data.with_columns([
            pl.lit(cct_SV if cct_SV else "-").alias("Cosmic_Fusion_Counts"),
            pl.lit("<=>".join(ccA)).alias("CC_Chr_Band"),
            pl.lit("<=>".join(ccB)).alias("CC_Tumour_Types(Somatic)"),
            pl.lit("<=>".join(ccC)).alias("CC_Cancer_Syndrome"),
            pl.lit("<=>".join(ccD)).alias("CC_Mutation_Type"),
            pl.lit("<=>".join(ccE)).alias("CC_Translocation_Partner")
        ])

        # DGvData
        try:
            dgv_loc1, dgv_loc2 = afd.AnnotateDGv(verbose, count, row, dgv_dict)
        except Exception as e:
            log.error(f"Error annotating DGV data: {e}")
            dgv_loc1, dgv_loc2 = [], []

        data = data.with_columns([
            pl.when(pl.lit(True))
            .then(pl.lit("<=>".join(dgv_loc1)))
            .otherwise(pl.col("DGv_Name-DGv_VarType-site1"))
            .alias("DGv_Name-DGv_VarType-site1")
        ])
        data = data.with_columns([
            pl.when(pl.lit(True))
            .then(pl.lit("<=>".join(dgv_loc2)))
            .otherwise(pl.col("DGv_Name-DGv_VarType-site2"))
            .alias("DGv_Name-DGv_VarType-site2")
        ])

    # Print to TSV file
    try:
        df = data.to_pandas()
        df.to_csv(out_file_txt, sep='\t', index=False)
        log.info(f"Wrote output to: {out_file_txt}")
    except Exception as e:
        log.error(f"Error writing TSV file: {e}")

    # Print to Json
    try:
        df.to_json(out_file_json)
        log.info(f"Wrote output to: {out_file_json}")
    except Exception as e:
        log.error(f"Error writing JSON file: {e}")

    # Print to Excel
    try:
        df.to_excel(out_file_exl, sheet_name='Annotated_SVs', index=False)
        log.info(f"Wrote output to: {out_file_exl}")
    except Exception as e:
        log.error(f"Error writing Excel file: {e}")

    end_time = time.time()
    log.info(f"Elapsed time was {end_time - start_time:.2f} seconds")


def ReadSVFile(filename, verbose):
    """Reads a structural variant file and returns a Polars DataFrame."""
    if verbose:
        log.info("Reading Structural Variant File")
    try:
        data = pl.read_csv(filename, sep='\t', infer_schema_length=1000)
        return data
    except Exception as e:
        log.error(f"Error reading SV file: {e}")
        raise


if __name__ == "__main__":
    app()