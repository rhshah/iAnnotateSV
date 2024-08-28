"""
Created on 12/23/2015
@Ronak Shah

"""

from collections import defaultdict
import pandas as pd
import logging
import coloredlogs

# Initialize logging
coloredlogs.install(level="DEBUG")


def ReadDGvFile(filename, verbose):
    """
    The function `ReadDGvFile` reads a TSV file, processes the data, and stores it in a dictionary.
    
    :param filename: The `filename` parameter in the `ReadDGvFile` function is the name of the file that
    you want to read and process. It should be a string representing the path to the file you want to
    open and read
    :param verbose: The `verbose` parameter in the `ReadDGvFile` function is a boolean flag that
    determines whether to display additional information or log messages during the execution of the
    function. If `verbose` is set to `True`, the function will log a message indicating that it is
    reading and storing a D
    :return: The function `ReadDGvFile` returns a dictionary where the keys are chromosome numbers
    (without the "chr" prefix) and the values are lists of data entries corresponding to each
    chromosome.
    """
    if verbose:
        logging.info(
            "iAnnotateSV::AnnotateForDGv: Reading & Storing DGV TSV file as dictionary"
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


def annotate_chromosome(sv_chr, sv_pos, dgvDict, verbose):
    """
    This Python function annotates a given chromosome and position with data from a dictionary based on
    specific conditions.
    
    :param sv_chr: The `sv_chr` parameter in the `annotate_chromosome` function represents the
    chromosome of the structural variant (SV) for which you want to annotate the position. It is a
    string that specifies the chromosome where the SV is located
    :param sv_pos: The `sv_pos` parameter represents the position of a structural variant on a specific
    chromosome. It is used in the `annotate_chromosome` function to find and annotate any relevant
    information from the DGV dictionary based on the chromosome and position provided
    :param dgvDict: `dgvDict` is a dictionary containing genomic locations from the Database of Genomic
    Variants (DGV). The keys in this dictionary are chromosome names (e.g., 'chr1', 'chr2') and the
    values are lists of genomic locations on that chromosome in the format 'start_position
    :param verbose: The `verbose` parameter in the `annotate_chromosome` function is a boolean flag that
    controls whether warning messages should be displayed during the annotation process. If `verbose` is
    set to `True`, warning messages will be logged using the `logging.warning` function. These messages
    provide additional information about the
    :return: The function `annotate_chromosome` returns a list of annotations for a given structural
    variant (SV) chromosome and position based on the provided DGV dictionary.
    """
    list_svloc = []
    list_loc = dgvDict.get(sv_chr, None)
    if list_loc:
        for loc in list_loc:
            data = loc.split("\t")
            dgv_pos1, dgv_pos2 = int(data[0]), int(data[1])
            if dgv_pos1 <= sv_pos <= dgv_pos2:
                joinedData = "-".join([data[2], data[8]])
                list_svloc.append(joinedData)
    else:
        if verbose:
            logging.warning(
                f"iAnnotateSV::AnnotateForDGv: Chromosome {sv_chr} is not there in the DGv dictionary"
            )
    return list_svloc


def AnnotateDGv(verbose, count, sv, dgvDict):
    """
    The function `AnnotateDGv` takes in verbose mode, count, structural variant information, and a
    dictionary of genomic variants, then annotates the genomic locations of the structural variant and
    returns the annotated locations.
    
    :param verbose: The `verbose` parameter is a boolean flag that determines whether to print detailed
    information or not during the execution of the `AnnotateDGv` function. If `verbose` is `True`,
    additional information will be logged using the `logging.info` function
    :param count: Count is the index of the entry being checked in the DGv data. It is used for logging
    purposes to track the progress of the annotation process
    :param sv: The `sv` parameter is a dictionary containing information about a structural variant
    (SV). It likely includes keys such as "chr1" for the chromosome of the first breakpoint, "pos1" for
    the position of the first breakpoint, "chr2" for the chromosome of the second breakpoint, and
    :param dgvDict: The function `AnnotateDGv` takes four parameters: `verbose`, `count`, `sv`, and
    `dgvDict`. Here is a brief description of each parameter:
    :return: The function `AnnotateDGv` is returning two lists, `list_svloc1` and `list_svloc2`, which
    contain the annotations for the genomic locations of the input structural variant `sv` based on the
    data in `dgvDict`.
    """
    if verbose:
        logging.info(f"iAnnotateSV::AnnotateForDGv: Checking Entry {count} in DGv data")

    sv_chr1, sv_pos1 = str(sv["chr1"]), int(sv["pos1"])
    sv_chr2, sv_pos2 = str(sv["chr2"]), int(sv["pos2"])

    list_svloc1 = annotate_chromosome(sv_chr1, sv_pos1, dgvDict, verbose)
    list_svloc2 = annotate_chromosome(sv_chr2, sv_pos2, dgvDict, verbose)

    return list_svloc1, list_svloc2
