"""
Created on  Mar 4, 2015

@author: Ronak H Shah
"""
from operator import itemgetter
import logging
import coloredlogs

"""
Get the mininmumIndex for a list and return the variables that match canonical transcript
Preference:# zone: 1=exon, 2=intron, 3=3'-UTR, 4=5'-UTR, 5=promoter
"""
coloredlogs.install(level="DEBUG")


def FindCT(
    geneList,
    transcriptList,
    siteList,
    zoneList,
    strandList,
    intronnumList,
    intronframeList,
    ctDict,
):
    """
    The function `FindCT` retrieves specific data elements based on input lists and a dictionary.

    :param geneList: The `FindCT` function takes in several lists (`geneList`, `transcriptList`,
    `siteList`, `zoneList`, `strandList`, `intronnumList`, `intronframeList`) and a dictionary `ctDict`
    as input parameters. The function then processes these inputs
    :param transcriptList: The `transcriptList` parameter in the `FindCT` function is a list containing
    transcript information. This list likely includes identifiers or names of transcripts associated
    with a gene. The function appears to iterate through the provided `geneList` to find the canonical
    transcript associated with a gene and returns various details
    :param siteList: The `siteList` parameter in the `FindCT` function likely refers to a list
    containing information about specific sites related to genetic sequences. This list could include
    details such as the location or position of certain genetic elements, mutations, or other relevant
    features within a gene or transcript. The function appears to
    :param zoneList: `zoneList` appears to be a list of zones. The function `FindCT` seems to be
    designed to find specific values corresponding to a given gene from the provided lists. If you have
    any specific questions or need further assistance with this function or any other part of the code,
    feel free to
    :param strandList: The `strandList` parameter in the `FindCT` function seems to be a list containing
    information about the strands of genes or transcripts. The function appears to iterate through the
    `geneList` to find the canonical transcript and then retrieves corresponding information from other
    lists based on the minimum index found in the
    :param intronnumList: The `intronnumList` parameter in the `FindCT` function seems to be a list
    containing information about intron numbers. The function is designed to find and return specific
    values based on input lists and dictionaries provided as arguments
    :param intronframeList: The `intronframeList` parameter in the `FindCT` function seems to be a list
    that contains information about the intron frame for each entry in the list. The function is
    designed to find and return specific details (gene, transcript, site, zone, strand, intron number,
    :param ctDict: The function `FindCT` takes in several lists as input parameters along with a
    dictionary `ctDict`. The purpose of this function is to find and return specific values based on the
    input lists and the dictionary
    :return: The function `FindCT` returns a tuple containing the values of `gene`, `transcript`,
    `site`, `zone`, `strand`, `intronnum`, and `intronframe`.
    """
    gene = None
    transcript = None
    site = None
    zone = None
    strand = None
    intronnum = None
    intronframe = None
    chkval = isinstance(geneList, list)
    # print geneList,transcriptList
    if chkval:
        cts = None
        minIndex = None
        for gene in geneList:
            # print gene
            if gene in ctDict:
                cts = ctDict.get(gene)
                break
        if cts:
            if len(cts) > 1:
                minIndex = min(enumerate(zoneList), key=itemgetter(1))[0]
            else:
                # print "CTS",cts[0]
                try:
                    # print "I am here", transcriptList.index(cts[0])
                    minIndex = transcriptList.index(cts[0])
                except ValueError:
                    logging.warn(
                        "iAnnotateSV::FindCanonicalTranscript: The given canonical transcript does not cover the coordinates."
                    )
                    minIndex = min(enumerate(zoneList), key=itemgetter(1))[0]
        else:
            minIndex = min(enumerate(zoneList), key=itemgetter(1))[0]
        # print minIndex
        gene = geneList[minIndex]
        transcript = transcriptList[minIndex]
        site = siteList[minIndex]
        zone = zoneList[minIndex]
        strand = strandList[minIndex]
        intronnum = intronnumList[minIndex]
        if intronnum == "Null":
            intronnum = None
        intronframe = intronframeList[minIndex]
        if intronframe == "Null":
            intronframe = None
    else:
        gene = geneList
        transcript = transcriptList
        site = siteList
        zone = zoneList
        strand = strandList
        intronnum = intronnumList
        intronframe = intronframeList
    return (gene, transcript, site, zone, strand, intronnum, intronframe)
