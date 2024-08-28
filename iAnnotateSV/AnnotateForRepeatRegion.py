"""
Created on 12/23/2015
@Ronak Shah

"""

from collections import defaultdict
import pandas as pd
import logging
import coloredlogs

coloredlogs.install(level="DEBUG")


def ReadRepeatFile(filename, verbose):
    """
    The function `ReadRepeatFile` reads a TSV file, processes the data, and stores it in a dictionary.
    
    :param filename: The `filename` parameter is a string that represents the path to the file that you
    want to read and process in the `ReadRepeatFile` function
    :param verbose: The `verbose` parameter in the `ReadRepeatFile` function is a boolean flag that
    determines whether to display additional information or messages during the execution of the
    function. If `verbose` is set to `True`, it will log a specific message using the `logging.info`
    function. If `verbose
    :return: The function `ReadRepeatFile` returns a dictionary where the keys are processed chromosome
    data (without the "chr" prefix) and the values are lists of joined data from the input file
    specified by the `filename` parameter.
    """
    if verbose:
        logging.info(
            "iAnnotateSV::AnnotateForRepeatRegion: Reading & Storing Repeat TSV file as dictionary"
        )
    dataDict = defaultdict(list)
    with open(filename, "r") as filecontent:
        header = filecontent.readline()
        for line in filecontent:
            data = line.rstrip("\n").split("\t")
            processedData = data[0].replace("chr", "")
            joinedData = "\t".join(data[1:])
            dataDict[processedData].append(joinedData)
    return dataDict


def check_overlap(sv_chr, sv_pos, rrDict, verbose):
    """
    The function `check_overlap` takes a chromosome, position, repeat region dictionary, and a verbosity
    flag as input, and returns a list of repeat regions that overlap with the given position on the
    chromosome.
    
    :param sv_chr: sv_chr is a string representing the chromosome of a structural variant
    :param sv_pos: The `sv_pos` parameter in the `check_overlap` function represents the position of a
    structural variant (SV) on a specific chromosome (`sv_chr`). The function checks for any overlap
    between this SV position and the positions of repeat regions (RR) stored in the `rrDict` dictionary
    for the
    :param rrDict: `rrDict` is a dictionary containing information about repeat regions. The keys are
    chromosome names (e.g., 'chr1', 'chr2') and the values are lists of strings representing the start
    and end positions of repeat regions on that chromosome, followed by additional data separated by
    tabs
    :param verbose: The `verbose` parameter in the `check_overlap` function is a boolean flag that
    controls whether additional information or logging messages should be displayed during the execution
    of the function. If `verbose` is set to `True`, the function will log a message using the
    `logging.info` function when a specific
    :return: The function `check_overlap` returns a list of strings containing data from the repeat
    regions that overlap with the given structural variant (SV) position on the specified chromosome.
    """
    list_svloc = []
    list_loc = rrDict.get(sv_chr)
    if list_loc:
        for loc in list_loc:
            data = loc.split("\t")
            rr_pos1, rr_pos2 = int(data[0]), int(data[1])
            if rr_pos1 <= sv_pos <= rr_pos2:
                joinedData = "-".join(data[2:])
                list_svloc.append(joinedData)
    else:
        if verbose:
            logging.info(
                "iAnnotateSV::AnnotateForRepeatRegion: Chromosome %s is not there in the repeat dictionary",
                sv_chr,
            )
    return list_svloc


def AnnotateRepeatRegion(verbose, count, sv, rrDict):
    """
    The function `AnnotateRepeatRegion` checks for overlap between two genomic positions and a
    dictionary of repeat regions.
    
    :param verbose: The `verbose` parameter is a boolean flag that determines whether to print detailed
    information or not during the execution of the function. If `verbose` is `True`, additional
    information will be logged using the `logging.info` function
    :param count: The `count` parameter in the `AnnotateRepeatRegion` function is used to keep track of
    the entry being processed in the repeat data. It is an integer value that represents the current
    entry number being checked in the repeat data
    :param sv: The `sv` parameter seems to be a DataFrame with columns "chr1" and "chr2" representing
    chromosome names, and columns "pos1" and "pos2" representing positions on those chromosomes
    :param rrDict: The `rrDict` parameter is likely a dictionary containing information about repeat
    regions. This dictionary is used in the `AnnotateRepeatRegion` function to check for overlap with
    the specified structural variant (sv) at positions sv_pos1 and sv_pos2 on chromosomes sv_chr1 and
    sv_chr2
    :return: The function `AnnotateRepeatRegion` is returning a tuple containing two lists:
    `list_svloc1` and `list_svloc2`.
    """
    if verbose:
        logging.info(
            "iAnnotateSV::AnnotateForRepeatRegion: Checking Entry %d in Repeat data",
            count,
        )
    sv_chr1, sv_pos1 = str(sv.loc["chr1"]), int(sv.loc["pos1"])
    sv_chr2, sv_pos2 = str(sv.loc["chr2"]), int(sv.loc["pos2"])
    list_svloc1 = check_overlap(sv_chr1, sv_pos1, rrDict, verbose)
    list_svloc2 = check_overlap(sv_chr2, sv_pos2, rrDict, verbose)
    return (list_svloc1, list_svloc2)
