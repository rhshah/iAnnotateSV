"""
Created on 12/23/2015
@Ronak Shah

"""
import pandas as pd
import numpy as np
import logging
import coloredlogs

# Initialize logging
coloredlogs.install(level="DEBUG")


def getVar(searchList, indices):
    """
    The function `getVar` takes a list `searchList` and a list of indices `indices`, and returns a
    sublist of `searchList` based on the specified indices.
    
    :param searchList: The `searchList` parameter is a list or array containing the elements you want to
    search through
    :param indices: The `indices` parameter is a list of indices that you want to use to extract
    specific elements from the `searchList`. By providing a list of indices, you can retrieve the
    elements at those specific positions in the `searchList`
    :return: The `getVar` function takes in two parameters: `searchList` and `indices`. It converts
    `searchList` into a NumPy array, then uses the `indices` to select specific elements from the array.
    Finally, it converts the selected elements back into a Python list and returns them.
    """
    return np.array(searchList)[indices].tolist()


def AnnotateFromCosmicCensusFile(filename, verbose, count, sv):
    """
    The function `AnnotateFromCosmicCensusFile` reads data from a file, filters it based on gene names,
    processes the data, and returns a list of annotated information.
    
    :param filename: The `filename` parameter is the name of the file from which the data will be read
    for annotating the cosmic census
    :param verbose: The `verbose` parameter is a boolean flag that controls whether logging information
    will be displayed during the execution of the function. If `verbose` is set to `True`, informational
    messages will be logged, such as the entry being checked in the Cosmic file. If `verbose` is set to
    `False
    :param count: The `count` parameter in the `AnnotateFromCosmicCensusFile` function is used to keep
    track of the entry being processed in the Cosmic file. It is used for logging purposes when
    `verbose` mode is enabled
    :param sv: The `sv` parameter is a dictionary containing information about a structural variant
    (SV). It likely includes keys such as "gene1" and "gene2" which represent the genes involved in the
    SV. The function `AnnotateFromCosmicCensusFile` reads a file (specified by
    :return: The function `AnnotateFromCosmicCensusFile` returns a list of annotated data related to the
    genes specified in the input structural variant (sv) data. The data is extracted from a Cosmic
    Census file (specified by the input filename) and processed based on specific columns. The function
    iterates over the filtered data for the genes provided in the sv data, extracts relevant
    information, prefixes the
    """
    if verbose:
        logging.info(
            f"iAnnotateSV::AnnotateForCosmic: Checking Entry {count} in Cosmic"
        )

    sv_gene1 = str(sv["gene1"])
    sv_gene2 = str(sv["gene2"])

    df = pd.read_csv(filename, sep="\t")
    filtered_df = df[df.iloc[:, 0].isin([sv_gene1, sv_gene2])]

    list_ccData = []
    for _, row in filtered_df.iterrows():
        prefix = "site1:" if row.iloc[0] == sv_gene1 else "site2:"
        slicedData = getVar(row.tolist(), [4, 7, 9, 12, 13])
        slicedProcessedData = [
            f"{prefix}{sData}" if sData else " " for sData in slicedData
        ]
        joinedData = "\t".join(slicedProcessedData)
        list_ccData.append(joinedData)

    return list_ccData


def AnnotateFromCosmicFusionCountsFile(filename, verbose, count, sv):
    """
    The function `AnnotateFromCosmicFusionCountsFile` reads a file containing Cosmic fusion counts data,
    matches gene combinations with a given structural variant, and returns the corresponding count if a
    match is found.
    
    :param filename: The `filename` parameter is the name of the file containing the Cosmic Fusion
    Counts data that you want to read and process in the `AnnotateFromCosmicFusionCountsFile` function
    :param verbose: The `verbose` parameter is a boolean flag that determines whether to display
    detailed logging information during the execution of the function. If `verbose` is set to `True`,
    the function will log information about the progress of checking entries in the Cosmic Counts data.
    If `verbose` is set to `False
    :param count: The `count` parameter in the `AnnotateFromCosmicFusionCountsFile` function is used to
    keep track of the entry being processed in the Cosmic Counts data. It is used for logging purposes
    when `verbose` is set to True. The function logs a message indicating which entry in
    :param sv: The `sv` parameter seems to be a dictionary containing information about a structural
    variant (SV). It likely includes keys such as "gene1" and "gene2" that represent genes involved in
    the SV
    :return: The function `AnnotateFromCosmicFusionCountsFile` returns an integer value representing the
    counts from the Cosmic Counts data for a specific gene combination in the input file. If there is no
    match found for the gene combination in the file, it returns `None`.
    """
    if verbose:
        logging.info(
            f"iAnnotateSV::AnnotateForCosmic: Checking Entry {count} in Cosmic Counts data"
        )

    sv_gene1 = str(sv["gene1"])
    sv_gene2 = str(sv["gene2"])
    sv_combos = {f"{sv_gene1}-{sv_gene2}", f"{sv_gene2}-{sv_gene1}"}

    countDF = pd.read_csv(filename, sep="\t", keep_default_na=True)
    countDF['Combo'] = countDF.apply(lambda row: {f"{row['Gene1']}-{row['Gene2']}", f"{row['Gene2']}-{row['Gene1']}"}, axis=1)
    match = countDF[countDF['Combo'].apply(lambda x: bool(sv_combos & x))]

    return None if match.empty else int(match.iloc[0]["Counts"])
