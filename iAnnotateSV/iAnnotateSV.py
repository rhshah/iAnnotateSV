import argparse
import time
import pandas as pd
import os
import sys
import logging
import coloredlogs

from . import helper as hp
from . import AnnotateEachBreakpoint as aeb
from . import PredictFunction as pf
from . import FindCanonicalTranscript as fct
from . import AddExternalAnnotations as aea
from . import AnnotationForKinaseDomain as kda
from . import VisualizeSV as vsv
from .models import *


def main(command=None):
    """
    The `main` function in this Python script annotates structural variants based on a specific human
    reference, with options for verbose output, selecting reference file version, output file settings,
    and additional annotations.
    
    :param command: The `main` function you provided is a Python script that defines a command-line
    interface using the `argparse` module. It sets up various command-line arguments for annotating
    structural variants based on a specific human reference
    """
    parser = argparse.ArgumentParser(
        prog="iAnnotateSV.py",
        description="Annotate SV based on a specific human reference",
        usage="%(prog)s [options]",
    )

    arguments = [
        (
            "-v",
            "--verbose",
            {"action": "store_true", "help": "make lots of noise [default]"},
        ),
        (
            "-r",
            "--refFileVersion",
            {
                "required": True,
                "metavar": "hg19",
                "help": "Which human reference file to be used, hg18,hg19 or hg38",
            },
        ),
        (
            "-rf",
            "--refFile",
            {
                "metavar": "hg19.sv.table.txt",
                "help": "Human reference file location to be used",
            },
        ),
        (
            "-ofp",
            "--outputFilePrefix",
            {"required": True, "metavar": "test", "help": "Prefix for the output file"},
        ),
        (
            "-o",
            "--outputDir",
            {
                "required": True,
                "metavar": "/somedir",
                "help": "Full Path to the output dir",
            },
        ),
        (
            "-i",
            "--svFile",
            {
                "required": True,
                "metavar": "svfile.txt",
                "help": "Location of the structural variants file to annotate",
            },
        ),
        (
            "-d",
            "--distance",
            {
                "default": 3000,
                "metavar": "3000",
                "help": "Distance used to extend the promoter region",
            },
        ),
        (
            "-a",
            "--autoSelect",
            {
                "action": "store_true",
                "default": True,
                "help": "Auto Select which transcript to be used[default]",
            },
        ),
        (
            "-c",
            "--canonicalTranscripts",
            {
                "metavar": "canonicalExons.txt",
                "help": "Location of canonical transcript list for each gene.",
            },
        ),
        (
            "-p",
            "--plotSV",
            {"action": "store_true", "help": "Plot the structural variant in question"},
        ),
        (
            "-u",
            "--uniprotFile",
            {
                "metavar": "uniprot.txt",
                "help": "Location of UniProt list contain information for protein domains.",
            },
        ),
        (
            "-rr",
            "--repeatFile",
            {
                "metavar": "RepeatRegionFile.tsv",
                "help": "Location of the Repeat Region Bed File",
            },
        ),
        (
            "-dgv",
            "--dgvFile",
            {
                "metavar": "DGvFile.tsv",
                "help": "Location of the Database of Genomic Variants Bed File",
            },
        ),
        (
            "-cc",
            "--cosmicConsensusFile",
            {
                "metavar": "CosmicConsensus.tsv",
                "help": "Location of the Cosmic Consensus TSV file",
            },
        ),
        (
            "-cct",
            "--cosmicCountsFile",
            {
                "metavar": "cosmic_fusion_counts.tsv",
                "help": "Location of the Cosmic Counts TSV file",
            },
        ),
    ]

    for arg in arguments:
        parser.add_argument(arg[0], arg[1], **arg[2])

    args = parser.parse_args(command.split()) if command else parser.parse_args()

    loggeroutput = os.path.join(
        args.outputDir, args.outputFilePrefix + "_iAnnotateSV.log"
    )
    logging.basicConfig(
        filename=loggeroutput,
        filemode="w",
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logging.DEBUG,
    )
    coloredlogs.install(level="DEBUG")

    this_dir = os.path.dirname(os.path.abspath(__file__))

    default_files = {
        "refFile": os.path.join(
            this_dir, "data/references", args.refFileVersion + ".sv.table.txt"
        ),
        "rrFilename": os.path.join(
            this_dir, "data/repeat_region", args.refFileVersion + "_repeatRegion.tsv"
        ),
        "dgvFilename": os.path.join(
            this_dir,
            "data/database_of_genomic_variants",
            args.refFileVersion + "_DGv_Annotation.tsv",
        ),
        "ccFilename": os.path.join(this_dir, "data/cosmic", "cancer_gene_census.tsv"),
        "cctFilename": os.path.join(
            this_dir, "data/cosmic", "cosmic_fusion_counts.tsv"
        ),
        "uniprot": os.path.join(
            this_dir,
            "data/UcscUniprotdomainInfo",
            args.refFileVersion + ".uniprot.spAnnot.table.txt",
        ),
    }

    for key, default_path in default_files.items():
        if not getattr(args, key):
            setattr(args, key, default_path)

    if args.canonicalTranscripts:
        args.autoSelect = False

    if args.refFileVersion not in ["hg18", "hg19", "hg38"]:
        logging.fatal(
            "iAnnotateSV: Please enter correct reference file version. Values can be: hg18 or hg19 or hg38"
        )
        sys.exit()

    refDF = hp.ReadFile(args.refFile)
    NewRefDF = hp.ExtendPromoterRegion(refDF, args.distance)
    svDF = hp.ReadFile(args.svFile)
    annDF = processSV(svDF, NewRefDF, args)
    plotDF = annDF.copy()

    outFilePrefixPath = os.path.join(
        args.outputDir, args.outputFilePrefix + "_functional.txt"
    )
    annDF.to_csv(outFilePrefixPath, sep="\t", index=False)

    if args.verbose:
        logging.info("iAnnotateSV: Adding External Annotations...")
    makeCommandLineForAEA = (
        f"-r {args.rrFilename} -d {args.dgvFilename} -c {args.ccFilename} -cct {args.cctFilename} "
        f"-s {outFilePrefixPath} -ofp {args.outputFilePrefix}_Annotated -o {args.outputDir}"
    )
    aea.main(makeCommandLineForAEA)

    if args.plotSV:
        if args.verbose:
            logging.info("iAnnotateSV: Plotting Each Structural Variants")
        plotSV(plotDF, NewRefDF, args.uniprot, args)

    if args.verbose:
        logging.info("iAnnotateSV: Finished Running the Annotation Process!!!")


def processSV(svDF, refDF, args):
    """
    The function `processSV` annotates structural variants in a DataFrame based on provided reference
    data and arguments, with an option to use canonical transcripts for annotation.
    
    :param svDF: `svDF` is a DataFrame containing structural variant data with columns such as "chr1",
    "pos1", "str1", "chr2", "pos2", "str2", etc. It seems to represent information about structural
    variants in a genomic context
    :param refDF: `refDF` seems to be a reference DataFrame containing information about genomic
    positions and transcripts. It is likely used for annotating structural variants in the `processSV`
    function. If you have any specific questions or need further assistance with this code snippet, feel
    free to ask!
    :param args: The `args` parameter in the `processSV` function seems to be a namespace object that
    contains various arguments or options passed to the function. These arguments could be used to
    control the behavior of the function based on the user's input
    :return: The function `processSV` returns either the annotated DataFrame `svDF` if
    `args.canonicalTranscripts` is not provided, or the result of running a function `kda.run` on the
    annotated DataFrame `annDF` if `args.canonicalTranscripts` is provided.
    """
    if args.verbose:
        logging.info("iAnnotateSV: Processing Each Structural Variants...")

    if args.canonicalTranscripts:
        ctDict = hp.ReadTranscriptFile(args.canonicalTranscripts)

    annDF = pd.DataFrame(
        columns=[
            "chr1",
            "pos1",
            "str1",
            "chr2",
            "pos2",
            "str2",
            "gene1",
            "transcript1",
            "site1",
            "gene2",
            "transcript2",
            "site2",
            "fusion",
        ]
    )

    def annotate_row(row):
        chr1, chr2 = str(row["chr1"]), str(row["chr2"])
        pos1, pos2 = int(row["pos1"]), int(row["pos2"])
        str1, str2 = int(row["str1"]), int(row["str2"])

        if args.autoSelect:
            gene1, transcript1, site1 = annotate_breakpoint(
                chr1, pos1, str1, refDF, args.autoSelect
            )
            gene2, transcript2, site2 = annotate_breakpoint(
                chr2, pos2, str2, refDF, args.autoSelect
            )
        else:
            gene1, transcript1, site1 = find_canonical_transcript(
                chr1, pos1, str1, refDF, ctDict
            )
            gene2, transcript2, site2 = find_canonical_transcript(
                chr2, pos2, str2, refDF, ctDict
            )

        fusionFunction = pf.PredictFunctionForSV(
            gene1, transcript1, site1, gene2, transcript2, site2
        )
        return pd.Series(
            [
                chr1,
                pos1,
                str1,
                chr2,
                pos2,
                str2,
                gene1,
                transcript1,
                site1,
                gene2,
                transcript2,
                site2,
                fusionFunction,
            ],
            index=annDF.columns,
        )

    annDF = svDF.apply(annotate_row, axis=1)

    if args.canonicalTranscripts:
        svDF = kda.run(
            annDF,
            args.refFile,
            args.canonicalTranscripts,
            args.allCanonicalTranscriptsPath,
            args.uniprot,
            args.verbose,
        )
        return svDF
    else:
        return annDF


def annotate_breakpoint(chr, pos, strand, refDF, autoSelect):
    """
    The function `annotate_breakpoint` takes chromosome, position, strand, reference dataframe, and
    autoSelect as input parameters, and returns gene, transcript, and site information after annotating
    the breakpoint.
    
    :param chr: The `chr` parameter likely refers to the chromosome on which the breakpoint is located.
    It could be a string representing the chromosome number or identifier
    :param pos: The `pos` parameter in the `annotate_breakpoint` function likely refers to the position
    of a breakpoint on a chromosome. This position is a numerical value that indicates the location of a
    genetic event, such as a DNA rearrangement or mutation, along the chromosome. The `pos` parameter is
    used
    :param strand: The `strand` parameter in the `annotate_breakpoint` function is used to specify the
    orientation of the DNA strand at the breakpoint location. It indicates whether the DNA strand is the
    forward strand (+) or the reverse strand (-) at that particular genomic position. This information
    is important for correctly annotating
    :param refDF: refDF is likely a reference DataFrame that contains information about genes,
    transcripts, sites, zones, intron numbers, and intron frames. This DataFrame is used to annotate
    breakpoints in the code snippet provided
    :param autoSelect: AutoSelect is a boolean parameter that determines whether to automatically select
    the best annotation for the breakpoint or not. If set to True, the function will automatically
    choose the best annotation based on certain criteria. If set to False, the user will have to
    manually select the annotation
    :return: The function `annotate_breakpoint` is returning the variables `gene`, `transcript`, and
    `site`.
    """
    (
        gene,
        transcript,
        site,
        zone,
        strand,
        intronnum,
        intronframe,
    ) = aeb.AnnotateEachBreakpoint(chr, pos, strand, refDF, autoSelect)
    return gene, transcript, site


def find_canonical_transcript(chr, pos, strand, refDF, ctDict):
    """
    The function `find_canonical_transcript` takes chromosome, position, strand, reference data frame,
    and a dictionary as input, annotates breakpoints, and returns gene, transcript, and site
    information.
    
    :param chr: The `chr` parameter in the `find_canonical_transcript` function is typically a string
    representing the chromosome on which the genetic variant is located. It is used to specify the
    chromosome for which the canonical transcript needs to be found
    :param pos: It seems like you were about to provide the explanation for the 'pos' parameter in the
    function 'find_canonical_transcript', but the explanation is missing. Could you please provide more
    details or complete the information about the 'pos' parameter so that I can assist you further?
    :param strand: The `strand` parameter in the `find_canonical_transcript` function is used to specify
    the orientation of the genomic region being analyzed. It indicates whether the region is located on
    the positive strand (+) or the negative strand (-) of the DNA molecule. This information is
    important for correctly identifying genes
    :param refDF: RefDF seems to be a reference DataFrame that is used in the function to annotate
    breakpoints. It is likely to contain information related to genomic positions, genes, transcripts,
    and other relevant data needed for the annotation process
    :param ctDict: The `ctDict` parameter in the `find_canonical_transcript` function is likely a
    dictionary that stores information about canonical transcripts. This dictionary is used in the
    function to find the canonical transcript based on the input chromosome, position, and strand
    :return: The function `find_canonical_transcript` returns the gene, transcript, and site as a tuple.
    """
    try:
        (
            geneList,
            transcriptList,
            siteList,
            zoneList,
            strandList,
            intronnumList,
            intronframeList,
        ) = aeb.AnnotateEachBreakpoint(chr, pos, strand, refDF, False)
        gene, transcript, site, zone, strand, intronnum, intronframe = fct.FindCT(
            geneList,
            transcriptList,
            siteList,
            zoneList,
            strandList,
            intronnumList,
            intronframeList,
            ctDict,
        )
    except (IntergenicError, ChrError) as e:
        logging.info("iAnnotateSV: " + str(e))
        gene, transcript, site, zone, strand, intronnum, intronframe = (
            "-",
            "-",
            "-",
            "-",
            "-",
            "-",
            "-",
        )
    return gene, transcript, site


def plotSV(svDF, refDF, uniprotPath, args):
    """
    The function `plotSV` plots each structural variant using input dataframes and a UniProt file.
    
    :param svDF: `svDF` is likely a DataFrame containing information about structural variants
    :param refDF: The `refDF` parameter in the `plotSV` function likely refers to a DataFrame containing
    reference data that is used for plotting structural variants. This data may include information such
    as genomic coordinates, gene annotations, or other relevant details needed for visualization and
    analysis of structural variants
    :param uniprotPath: The `uniprotPath` parameter in the `plotSV` function is a file path that should
    point to a file containing UniProt data. This file is used to plot structural variants in the
    function
    :param args: The `args` parameter in the `plotSV` function likely contains various arguments or
    options that control the behavior of the function. These arguments could include settings for
    verbosity (such as `verbose`), file paths, visualization options, or any other parameters that
    affect how the structural variants are plotted
    """
    if args.verbose:
        logging.info("iAnnotateSV: Will now try to plot Each Structural Variants")
    if not os.path.isfile(uniprotPath):
        logging.fatal(
            "iAnnotateSV: %s file does not exist!!, Please use it to plot structural variants",
            uniprotPath,
        )
        sys.exit()

    upDF = hp.ReadFile(uniprotPath)
    vsv.VisualizeSV(svDF, refDF, upDF, args)


if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    logging.info("iAnnotateSV: Elapsed time was %g seconds", (end_time - start_time))
