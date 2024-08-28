"""
Created on 12/23/2015.

@Ronak Shah

"""

import sys
import argparse
import pandas as pd
import time
import logging
import AnnotateForRepeatRegion as afr
import AnnotateForCosmic as afc
import AnnotateForDGv as afd
import coloredlogs

# Initialize logging
coloredlogs.install(level="DEBUG")


def main(command=None):
    """
    The `main` function in this Python script parses command line arguments, reads input files,
    annotates structural variants with external data, and outputs the annotated data to text, JSON, and
    Excel files.

    :param command: The `main` function you provided is a script that takes command-line arguments using
    the `argparse` module in Python. It defines several command-line arguments such as input file
    locations, output file paths, and verbosity level
    """
    parser = argparse.ArgumentParser(
        prog="AddExternalAnnotations.py",
        description="Add External Annotation to the Structural Variants",
        usage="%(prog)s [options]",
    )
    parser.add_argument(
        "-r",
        "--repeatFile",
        action="store",
        dest="rrFilename",
        required=True,
        metavar="RepeatRegionFile.tsv",
        help="Location of the Repeat Region Bed File",
    )
    parser.add_argument(
        "-d",
        "--dgvFile",
        action="store",
        dest="dgvFilename",
        required=True,
        metavar="DGvFile.tsv",
        help="Location of the Database of Genomic Variants Bed File",
    )
    parser.add_argument(
        "-c",
        "--cosmicConsensusFile",
        action="store",
        dest="ccFilename",
        required=True,
        metavar="CosmicConsensus.tsv",
        help="Location of the Cosmic Consensus TSV file",
    )
    parser.add_argument(
        "-cct",
        "--cosmicCountsFile",
        action="store",
        dest="cctFilename",
        required=True,
        metavar="cosmic_fusion_counts.tsv",
        help="Location of the Cosmic Counts TSV file",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        default=True,
        help="make lots of noise [default]",
    )
    parser.add_argument(
        "-s",
        "--svFile",
        action="store",
        dest="svFilename",
        required=True,
        metavar="SVfile.txt",
        help="Location of the structural variant file to be annotated",
    )
    parser.add_argument(
        "-ofp",
        "--outputFilePrefix",
        action="store",
        dest="outFilePrefix",
        required=True,
        metavar="AnnotatedSV",
        help="Full path with prefix name for the output file",
    )
    parser.add_argument(
        "-o",
        "--outputDir",
        action="store",
        dest="outDir",
        required=True,
        metavar="/somedir",
        help="Full Path to the output dir",
    )
    args = parser.parse_args(command.split()) if command else parser.parse_args()

    outFileTxt = f"{args.outDir}/{args.outFilePrefix}.txt"
    outFileExl = f"{args.outDir}/{args.outFilePrefix}.xlsx"
    outFileJson = f"{args.outDir}/{args.outFilePrefix}.json"

    data = ReadSVFile(args)
    repeatRegionDict = afr.ReadRepeatFile(args.rrFilename, args.verbose)
    dgvDict = afd.ReadDGvFile(args.dgvFilename, args.verbose)

    initialize_columns(data)

    for count, row in data.iterrows():
        log_info(
            args.verbose,
            f"Processing Record: {row['chr1']}\t{row['pos1']}\t{row['chr2']}\t{row['pos2']}\t{row['gene1']}\t{row['gene2']}",
        )

        # Repeat Region Data
        rr_loc1, rr_loc2 = afr.AnnotateRepeatRegion(
            args.verbose, count, row.copy(), repeatRegionDict
        )
        data.loc[count, "repName-repClass-repFamily:-site1"] = "<=>".join(rr_loc1)
        data.loc[count, "repName-repClass-repFamily:-site2"] = "<=>".join(rr_loc2)

        # Cosmic Consensus Data
        cc_SV = afc.AnnotateFromCosmicCensusFile(
            args.ccFilename, args.verbose, count, row.copy()
        )
        cct_SV = afc.AnnotateFromCosmicFusionCountsFile(
            args.cctFilename, args.verbose, count, row.copy()
        )
        ccA, ccB, ccC, ccD, ccE = ([] for _ in range(5))
        for cc in cc_SV:
            ccData = cc.split("\t")
            ccA.append(ccData[0])
            ccB.append(ccData[1])
            ccC.append(ccData[2])
            ccD.append(ccData[3])
            ccE.append(ccData[4])
        data.loc[count, "Cosmic_Fusion_Counts"] = cct_SV
        data.loc[count, "CC_Chr_Band"] = "<=>".join(ccA)
        data.loc[count, "CC_Tumour_Types(Somatic)"] = "<=>".join(ccB)
        data.loc[count, "CC_Cancer_Syndrome"] = "<=>".join(ccC)
        data.loc[count, "CC_Mutation_Type"] = "<=>".join(ccD)
        data.loc[count, "CC_Translocation_Partner"] = "<=>".join(ccE)

        # DGvData
        dgv_loc1, dgv_loc2 = afd.AnnotateDGv(args.verbose, count, row.copy(), dgvDict)
        data.loc[count, "DGv_Name-DGv_VarType-site1"] = "<=>".join(dgv_loc1)
        data.loc[count, "DGv_Name-DGv_VarType-site2"] = "<=>".join(dgv_loc2)

    # Print to TSV file
    data.to_csv(outFileTxt, sep="\t", index=False)
    # Print to Json
    data.to_json(outFileJson)
    # Print to Excel
    data.to_excel(outFileExl, sheet_name="Annotated_SVs", index=False)


def ReadSVFile(args):
    """
    The function `ReadSVFile` reads a structural variant file, checks if it is empty, and saves the data
    to different file formats if not empty.

    :param args: The `args` parameter in the `ReadSVFile` function seems to be an object or dictionary
    containing various attributes such as `svFilename`, `outDir`, and `outFilePrefix`. These attributes
    are likely used to specify file paths and names for input and output files in the function
    :return: The function `ReadSVFile` is returning the data read from the Structural Variant File in
    the form of a pandas DataFrame.
    """
    logging.info("Reading Structural Variant File")
    data = pd.read_csv(args.svFilename, sep="\t", header=0)
    if data.empty:
        logging.warning(
            f"File {args.svFilename} does not have any structural variants to annotate."
        )
        outFileTxt = f"{args.outDir}/{args.outFilePrefix}.txt"
        outFileExl = f"{args.outDir}/{args.outFilePrefix}.xlsx"
        outFileJson = f"{args.outDir}/{args.outFilePrefix}.json"
        data.to_csv(outFileTxt, sep="\t", index=False)
        data.to_excel(outFileExl, sheet_name="Annotated_SVs", index=False)
        data.to_json(outFileJson)
        sys.exit()
    return data


def initialize_columns(data):
    """
    The function `initialize_columns` initializes specific columns in a given dataset with a default
    value of "-".

    :param data: The `initialize_columns` function takes a `data` parameter, which is presumably a data
    structure like a dictionary or a DataFrame where you want to initialize specific columns with a
    default value "-"
    """
    columns = [
        "Cosmic_Fusion_Counts",
        "repName-repClass-repFamily:-site1",
        "repName-repClass-repFamily:-site2",
        "CC_Chr_Band",
        "CC_Tumour_Types(Somatic)",
        "CC_Cancer_Syndrome",
        "CC_Mutation_Type",
        "CC_Translocation_Partner",
        "DGv_Name-DGv_VarType-site1",
        "DGv_Name-DGv_VarType-site2",
    ]
    for col in columns:
        data[col] = "-"


def log_info(verbose, message):
    """
    The function `log_info` logs a message if the `verbose` flag is set to true.

    :param verbose: A boolean parameter that determines whether to display detailed information or not.
    If set to True, it will log the message; if set to False, it will not log the message
    :param message: The `message` parameter in the `log_info` function is a string that represents the
    information or message that you want to log. This message will be logged using the `logging.info`
    method if the `verbose` parameter is set to `True`
    """
    if verbose:
        logging.info(message)


if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    print(f"Elapsed time was {end_time - start_time} seconds")
