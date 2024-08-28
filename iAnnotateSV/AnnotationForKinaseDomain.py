"""
Created on 01/09/2018
@Ronak Shah

"""
import os
import sys
import pandas as pd
import logging
import coloredlogs
import re
import helper as hp

coloredlogs.install(level="DEBUG")


def read_file_or_exit(path, message):
    """
    The function reads a file at a specified path or exits the program with a message if the file does
    not exist.

    :param path: The `path` parameter in the `read_file_or_exit` function is a string that represents
    the file path of the file that you want to read
    :param message: The `message` parameter in the `read_file_or_exit` function is a string that
    contains the error message to be logged if the file does not exist at the specified `path`. This
    message will be critical level logged if the `verbose` flag is set to True. If the file does not
    :return: the result of calling `hp.ReadFile(path)` if the file at the specified `path` exists.
    """
    if os.path.isfile(path):
        return hp.ReadFile(path)
    if verbose:
        logging.critical(message)
    sys.exit(1)


def get_transcript(gene):
    """
    The function `get_transcript` retrieves the transcript information for a given gene from a
    DataFrame.

    :param gene: The `get_transcript` function you provided seems to be designed to retrieve the
    transcript information for a given gene from two different dataframes (`ctDF` and `allctDF`). If the
    gene is not found in the first dataframe, it attempts to find it in the second dataframe
    :return: The `get_transcript` function is designed to return the transcript associated with a given
    gene. It first tries to locate the gene in the `ctDF` DataFrame and retrieve the corresponding
    transcript. If the gene is not found in `ctDF`, it then tries to find the gene in the `allctDF`
    DataFrame and retrieve the transcript from there. If the gene is not found in either
    """
    try:
        return ctDF.loc[ctDF.Gene == gene, "Transcripts"].iat[0]
    except IndexError:
        try:
            return allctDF.loc[allctDF.Gene == gene, "Transcripts"].iat[0]
        except IndexError:
            return None


def run(svDFA, refPath, ctPath, allctPath, upPath, verbose):
    """
    The function `run` reads input files, annotates kinase domains for structural variants, and returns
    the annotated data frame.

    :param svDFA: It looks like the `run` function you provided is responsible for annotating kinase
    domains in a DataFrame containing structural variant data. The function reads several input files,
    processes the data, and annotates kinase domains based on the information provided
    :param refPath: The `refPath` parameter in the `run` function represents the location of the
    reference-based annotation file. This file is used to retrieve information needed for annotating
    kinase domains in the input data. The function reads this file and processes its contents to
    annotate kinase domains in the structural variant data
    :param ctPath: The `ctPath` parameter in the `run` function refers to the location of the
    assay-specific canonical transcript file. This file is used in the function to read data from and
    perform annotations related to kinase domains
    :param allctPath: The `allctPath` parameter in the `run` function is used to specify the location of
    the file containing all canonical transcript information. This file is read using the
    `read_file_or_exit` function to ensure that the path is correct. The data from this file is stored
    in the `all
    :param upPath: The `upPath` parameter is the location of the Uniprot Annotation file. This file
    contains information about protein sequences, functions, and domains. The function `run` reads this
    file and uses the information to annotate kinase domains in a given dataset of structural variations
    :param verbose: The `verbose` parameter in the `run` function is a boolean flag that controls
    whether additional information and logging messages should be displayed during the execution of the
    function. If `verbose` is set to `True`, the function will log messages to provide updates on the
    progress of the annotation process. If
    :return: The function `run` is returning a DataFrame `svDF` after performing certain operations on
    it. The DataFrame `svDF` is being modified within the function by inserting two new columns
    "kinase_domain1" and "kinase_domain2" and populating them with values based on the input data and
    calculations performed within the function. The function iterates over the rows of the input
    DataFrame `
    """

    upDF = read_file_or_exit(
        upPath,
        "iAnnotateSV::AnnotationForKinaseDomain: Location of Uniprot Annotation file is incorrect!!!",
    )
    ctDF = read_file_or_exit(
        ctPath,
        "iAnnotateSV::AnnotationForKinaseDomain: Location of assay specific canonical transcript file is incorrect!!!",
    )
    allctDF = read_file_or_exit(
        allctPath,
        "iAnnotateSV::AnnotationForKinaseDomain: Location of all canonical transcript file is incorrect!!!",
    )
    refDF = read_file_or_exit(
        refPath,
        "iAnnotateSV::AnnotationForKinaseDomain: Location of reference based annotation file is incorrect!!!",
    )
    refDF.columns = refDF.columns.str.replace("#", "")

    svDF = svDFA.copy()
    svDF.insert(loc=9, column="kinase_domain1", value=None)
    svDF.insert(loc=13, column="kinase_domain2", value=None)

    for count, row in svDFA.iterrows():
        if verbose:
            logging.info(
                "iAnnotateSV::AnnotateForKinaseDomain: Checking Entry %d in Uniprot data",
                count,
            )

        chr1 = (
            f"chr{row.loc['chr1']}"
            if not str(row.loc["chr1"]).startswith("chr")
            else str(row.loc["chr1"])
        )
        chr2 = (
            f"chr{row.loc['chr2']}"
            if not str(row.loc["chr2"]).startswith("chr")
            else str(row.loc["chr2"])
        )
        pos1, pos2 = int(row.loc["pos1"]), int(row.loc["pos2"])
        gene1, gene2 = str(row.loc["gene1"]), str(row.loc["gene2"])
        site1, site2 = str(row.loc["site1"]), str(row.loc["site2"])

        transcript1 = get_transcript(gene1)
        transcript2 = get_transcript(gene2)
        fusion = str(row.loc["fusion"])

        kanno1 = kanno2 = None

        if fusion != "-":
            fusionevent = re.search(r"\{(.*)\}", fusion)
            if fusionevent:
                eventType = fusionevent.group(1)
                if ":" in eventType:
                    egene1, egene2 = eventType.split(":")
                    if transcript1:
                        kanno1 = getKinaseInfo(
                            chr1, pos1, gene1, egene1, egene2, transcript1, refDF, upDF
                        )
                    if transcript2:
                        kanno2 = getKinaseInfo(
                            chr2, pos2, gene2, egene1, egene2, transcript2, refDF, upDF
                        )
            svDF.at[count, "kinase_domain1"] = kanno1
            svDF.at[count, "kinase_domain2"] = kanno2

    return svDF


def processData(chrom, transcript, refDF, upDF):
    """
    The function `processData` retrieves specific data based on chromosome and transcript information
    from reference and upstream dataframes.

    :param chrom: Chrom is a variable representing the chromosome number or identifier. It is used to
    specify the chromosome for which the data processing is being done in the given function
    `processData`
    :param transcript: It looks like the code snippet you provided is a function `processData` that
    takes in several parameters to process some data related to transcripts and reference dataframes.
    However, you didn't provide the value for the `transcript` parameter
    :param refDF: refDF is a DataFrame containing reference data with columns like "name", "chrom",
    "txStart", and "txEnd". It seems to store information about transcripts
    :param upDF: The `upDF` parameter in the `processData` function seems to be a DataFrame containing
    genomic data related to transcripts. The function processes this data based on the provided `chrom`
    and `transcript` parameters along with reference data stored in `refDF`
    :return: The function `processData` returns three values: `up_recordIndex`, `max_len`, and
    `min_len`.
    """
    transcripts = refDF[refDF["name"] == transcript]
    transcriptIdx = (
        getValueOrDefault(transcripts[transcripts["chrom"] == chrom].index, 0)
        if len(transcripts) > 1
        else getValueOrDefault(refDF[refDF["name"] == transcript].index, 0)
    )

    refTxSt = int(refDF.at[transcriptIdx, "txStart"])
    refTxEn = int(refDF.at[transcriptIdx, "txEnd"])

    up_idxList = upDF[upDF["#chrom"] == chrom].index
    up_recordIndex = [
        index
        for index in up_idxList
        if upDF.at[index, "chromStart"] >= refTxSt
        and upDF.at[index, "chromEnd"] <= refTxEn
        and upDF.at[index, "annotationType"] == "domain"
    ]

    allMaxVal = [max(refTxEn, upDF.at[val, "chromEnd"]) for val in up_recordIndex]
    allMinVal = [min(refTxSt, upDF.at[val, "chromStart"]) for val in up_recordIndex]

    max_len = max(allMaxVal) if allMaxVal else refTxEn
    min_len = max(allMinVal) if allMinVal else refTxSt

    return up_recordIndex, max_len, min_len


def getKinaseInfo(chrom, pos, gene, egene1, egene2, transcript, refDF, upDF):
    """
    The function `getKinaseInfo` retrieves information about kinase domains based on genomic coordinates
    and gene annotations.

    :param chrom: Chromosome number or identifier where the gene is located
    :param pos: The `pos` parameter in the `getKinaseInfo` function represents the position on the
    chromosome where you want to check for the presence of a kinase domain. This function takes various
    parameters related to genomic data and processes them to determine if the specified position falls
    within a protein kinase domain. If it
    :param gene: The `gene` parameter in the `getKinaseInfo` function represents the gene for which you
    want to retrieve kinase domain information. This gene will be compared with `egene1` and `egene2` to
    determine if it matches either of them in order to check for the presence of
    :param egene1: `egene1` is a variable representing a gene name
    :param egene2: egene2 is a variable representing the second gene in a pair of genes
    :param transcript: Transcript is a parameter that represents the transcript associated with the gene
    for which you are retrieving kinase information. It is used to identify the specific transcript in
    the reference data frame and upstream data frame for further processing in the `getKinaseInfo`
    function
    :param refDF: The function `getKinaseInfo` takes several parameters including `chrom`, `pos`,
    `gene`, `egene1`, `egene2`, `transcript`, `refDF`, and `upDF`. It processes the data using the
    `processData` function and then checks if a
    :param upDF: The `upDF` parameter in the `getKinaseInfo` function likely refers to a DataFrame
    containing information about protein domains, including the start and end positions of each domain
    on the chromosome. This DataFrame is used to determine if a given position (`pos`) falls within a
    protein kinase domain for a
    :return: The function `getKinaseInfo` returns the result of checking if the provided position falls
    within a kinase domain for a given gene. If the position is within the kinase domain, it returns
    whether it is fully included, partially included, or not included in the kinase domain. If the
    position is not within a kinase domain or if the necessary data is not available, it returns `None`.
    """
    domainIdx, maxLen, minLen = processData(chrom, transcript, refDF, upDF)
    if domainIdx is None:
        return None

    strand = refDF.at[refDF[refDF["name"] == transcript].index[0], "strand"]
    kanno = None

    def check_kinase_domain(pos, chromStart, chromEnd):
        if pos > chromEnd:
            return "Kinase Domain Included"
        if chromStart <= pos <= chromEnd:
            return "Partial Kinase Domain Included"
        return "Kinase Domain Not Included"

    if strand == "+":
        if egene1 == gene:
            for val in domainIdx:
                if "Protein kinase" in upDF.at[val, "name"]:
                    return check_kinase_domain(
                        pos, upDF.at[val, "chromStart"], upDF.at[val, "chromEnd"]
                    )
        if egene2 == gene:
            for val in domainIdx:
                if "Protein kinase" in upDF.at[val, "name"]:
                    return check_kinase_domain(
                        pos, upDF.at[val, "chromStart"], upDF.at[val, "chromEnd"]
                    )
    else:
        if egene1 == gene:
            for val in domainIdx:
                if "Protein kinase" in upDF.at[val, "name"]:
                    return check_kinase_domain(
                        pos, upDF.at[val, "chromStart"], upDF.at[val, "chromEnd"]
                    )
        if egene2 == gene:
            for val in domainIdx:
                if "Protein kinase" in upDF.at[val, "name"]:
                    return check_kinase_domain(
                        pos, upDF.at[val, "chromStart"], upDF.at[val, "chromEnd"]
                    )

    return kanno


def getValueOrDefault(value, index, default=None):
    """
    The function `getValueOrDefault` returns the value at a specified index in a list or string, or a
    default value if an error occurs.

    :param value: The `value` parameter is the list or string from which you want to retrieve a value
    :param index: The `index` parameter in the `getValueOrDefault` function is used to specify the
    position of the element in the `value` parameter that you want to retrieve. It is the index of the
    element you want to access in the `value` list or string
    :param default: The `default` parameter in the `getValueOrDefault` function is a value that will be
    returned if an exception occurs while trying to access the element at the specified `index` in the
    `value` list. If no `default` value is provided when calling the function, it will default to `
    :return: The `getValueOrDefault` function is designed to return the value at the specified index in
    the `value` list. If the index is out of range or if an exception occurs during the retrieval, the
    function will return the `default` value provided (which is `None` by default).
    """
    try:
        return value[index]
    except Exception:
        return default
