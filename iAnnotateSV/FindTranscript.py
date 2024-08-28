"""

Created on Mar 4, 2015

@author: Ronak H Shah
"""

import helper as hp
import pandas as pd

"""
This function will return a single highest priority transcript

The Python script defines functions to find and return information about transcripts based on
query data and reference data.
    
:param queryDF: The `queryDF` parameter in the functions `FindATranscript` and `FindAllTranscripts`
is a DataFrame containing query data that is used to find transcripts. It likely includes
information about genomic locations, exon numbers, intron positions, and other relevant details for
transcripts
:param refDF: The `refDF` parameter in the functions `FindATranscript` and `FindAllTranscripts`
seems to be a reference DataFrame containing information about transcripts. This DataFrame likely
includes columns such as `#name`, `name2`, `strand`, `exonCount`, and other relevant information
:return: The `FindATranscript` function returns a tuple containing information about a single
highest priority transcript, including gene name, transcript name, description, zone, strand
direction, intron number, and intron frame.

"""


def FindATranscript(queryDF, refDF):
    
    """
    The function `FindATranscript` takes in two DataFrames, extracts relevant information based on
    specific conditions, and returns a tuple containing gene name, transcript name, description, zone,
    strand direction, intron number, and intron frame.
    
    :param queryDF: The function `FindATranscript` takes two parameters `queryDF` and `refDF`. The
    `queryDF` parameter is a DataFrame containing query data, and the `refDF` parameter is a reference
    DataFrame. The function processes the query data to find information about a specific transcript
    based on
    :param refDF: The function `FindATranscript` takes two parameters, `queryDF` and `refDF`. The
    `refDF` parameter is expected to be a DataFrame containing reference data for transcripts. This data
    likely includes information such as transcript names, gene names, strand directions, exon counts,
    and other
    :return: The function `FindATranscript` returns a tuple containing the following elements:
    - `geneName`: The name of the gene associated with the transcript
    - `transcript`: The name of the transcript
    - `desc`: A description of the location within the gene/transcript (e.g., Exon, Intron, UTR,
    Promoter)
    - `zone`: The zone/category of
    """
    tmpDF = queryDF.copy()
    desc = None
    intronnum = None
    intronframe = None
    if len(queryDF.index) > 1:
        try:
            tmpDF["d"] = pd.to_numeric(queryDF["d"], downcast="signed")
            idxMin = int(tmpDF["d"].idxmin(axis=1))
        except ValueError:
            idxMin = queryDF.index.values.tolist()[0]
        # print idxMin
    else:
        idxMin = queryDF.index.values.tolist()[0]
        # print "1",idxMin
    c = int(queryDF.loc[idxMin]["c"])

    if queryDF.loc[idxMin]["d"]:
        d = int(queryDF.loc[idxMin]["d"])
    else:
        d = queryDF.loc[idxMin]["d"]
    if queryDF.loc[idxMin]["d1"]:
        d1 = int(queryDF.loc[idxMin]["d1"])
    else:
        d1 = queryDF.loc[idxMin]["d1"]
    if queryDF.loc[idxMin]["d2"]:
        d2 = int(queryDF.loc[idxMin]["d2"])
    else:
        d2 = queryDF.loc[idxMin]["d2"]
    if queryDF.loc[idxMin]["e"]:
        e = int(queryDF.loc[idxMin]["e"])
    else:
        e = queryDF.loc[idxMin]["e"]
    if queryDF.loc[idxMin]["e1"]:
        e1 = int(queryDF.loc[idxMin]["e1"])
    else:
        e1 = queryDF.loc[idxMin]["e1"]
    if queryDF.loc[idxMin]["e2"]:
        e2 = int(queryDF.loc[idxMin]["e2"])
    else:
        e2 = queryDF.loc[idxMin]["e2"]
    f = queryDF.loc[idxMin]["c"]
    zone = c
    transcript = refDF.iloc[idxMin]["#name"]
    geneName = refDF.iloc[idxMin]["name2"]
    strandDirection = refDF.iloc[idxMin]["strand"]
    if strandDirection == "-":
        if e:
            e = int(refDF.iloc[idxMin]["exonCount"]) - e + 1
        else:
            e = int(refDF.iloc[idxMin]["exonCount"]) + 1
        if e1:
            e1 = int(refDF.iloc[idxMin]["exonCount"]) - e1 + 1
        else:
            e1 = int(refDF.iloc[idxMin]["exonCount"]) + 1
        if e2:
            e2 = int(refDF.iloc[idxMin]["exonCount"]) - e2 + 1
        else:
            e2 = int(refDF.iloc[idxMin]["exonCount"]) + 1
    # in Exon:
    if c == 1:
        desc = "Exon " + str(e) + " of " + geneName + "(" + strandDirection + ")"
        # print desc

    # In Intron
    elif c == 2:
        if strandDirection == "+":
            if d1 < d2:
                desc = (
                    "Intron of "
                    + geneName
                    + "("
                    + strandDirection
                    + "):"
                    + hp.bp2str(d1, 2)
                    + " after exon "
                    + str(e1)
                )
            else:
                desc = (
                    "Intron of "
                    + geneName
                    + "("
                    + strandDirection
                    + "):"
                    + hp.bp2str(d2, 2)
                    + " before exon "
                    + str(e2)
                )
            # print desc

        else:
            if d1 < d2:
                desc = (
                    "Intron of "
                    + geneName
                    + "("
                    + strandDirection
                    + "):"
                    + hp.bp2str(d1, 2)
                    + " before exon "
                    + str(e1)
                )
            else:
                desc = (
                    "Intron of "
                    + geneName
                    + "("
                    + strandDirection
                    + "):"
                    + hp.bp2str(d2, 2)
                    + " after exon "
                    + str(e2)
                )
            # print desc

        intronnum = e1
        intronframe = f
    # In 3'-UTR
    elif c == 3:
        desc = (
            "3'-UTR of "
            + geneName
            + "("
            + strandDirection
            + "):"
            + hp.bp2str(d, 2)
            + " after coding stop"
        )
        # print desc

    # In 5'-UTR
    elif c == 4:
        desc = (
            "5'-UTR of "
            + geneName
            + "("
            + strandDirection
            + "):"
            + hp.bp2str(d, 2)
            + " before coding start"
        )
        # print desc

    # In Promoter
    elif c == 5:
        desc = (
            "Promoter of "
            + geneName
            + "("
            + strandDirection
            + "):"
            + hp.bp2str(d, 2)
            + " from tx start"
        )
        # print desc

    else:
        desc = "Unexpected Error"
        # print desc
    return (geneName, transcript, desc, zone, strandDirection, intronnum, intronframe)


"""
This function will return all the transcripts
"""


def FindAllTranscripts(queryDF, refDF):
    """
    The function `FindAllTranscripts` processes data from two DataFrames to extract information about
    gene transcripts, including gene names, descriptions, and genomic locations.
    
    :param queryDF: The function `FindAllTranscripts` takes two parameters, `queryDF` and `refDF`. The
    `queryDF` parameter seems to be a DataFrame containing information related to gene transcripts,
    while the `refDF` parameter appears to be a reference DataFrame containing additional details about
    the gene transcripts
    :param refDF: The function `FindAllTranscripts` takes two parameters, `queryDF` and `refDF`. The
    `refDF` parameter seems to be a DataFrame containing information about gene transcripts, such as
    gene names, transcript names, strand directions, exon counts, etc
    :return: The function `FindAllTranscripts` returns multiple lists containing information about
    transcripts based on the input query dataframe and reference dataframe. The lists returned are:
    - `geneNameList`: List of gene names
    - `transcriptList`: List of transcript names
    - `descList`: List of descriptions for each transcript
    - `zoneList`: List of zones for each transcript
    - `strandDirection
    """
    geneNameList = []
    transcriptList = []
    strandDirectionList = []
    descList = []
    zoneList = []
    intronnumList = []
    intronframeList = []

    for count, row in queryDF.iterrows():
        desc = "Null"
        intronnum = "Null"
        intronframe = "Null"
        transcript = refDF.iloc[count]["#name"]
        transcriptList.append(transcript)
        geneName = refDF.iloc[count]["name2"]
        geneNameList.append(geneName)
        strandDirection = refDF.iloc[count]["strand"]
        strandDirectionList.append(strandDirection)
        # print transcript,geneName,strandDirection
        c = int(row.loc["c"])
        if row.loc["d"]:
            d = int(row.loc["d"])
        else:
            d1 = row.loc["d"]
        e = row.loc["e"]
        if row.loc["d1"]:
            d1 = int(row.loc["d1"])
        else:
            d1 = row.loc["d1"]
        if row.loc["d2"]:
            d2 = int(row.loc["d2"])
        else:
            d1 = row.loc["d2"]
        e1 = row.loc["e1"]
        e2 = row.loc["e2"]
        f = row.loc["f"]
        zone = c
        zoneList.append(zone)
        if strandDirection == "-":
            if e:
                e = int(refDF.iloc[count]["exonCount"]) - int(row.loc["e"]) + 1
            else:
                e = int(refDF.iloc[count]["exonCount"]) + 1
            if e1:
                e1 = int(refDF.iloc[count]["exonCount"]) - int(row.loc["e1"]) + 1
            else:
                e1 = int(refDF.iloc[count]["exonCount"]) + 1
            if e2:
                e2 = int(refDF.iloc[count]["exonCount"]) - int(row.loc["e2"]) + 1
            else:
                e2 = int(refDF.iloc[count]["exonCount"]) + 1
        # in Exon:
        if zone == 1:
            desc = "Exon " + str(e) + " of " + geneName + "(" + strandDirection + ")"
            descList.append(desc)
            intronnumList.append(intronnum)
            intronframeList.append(intronframe)
            # print desc
            continue
        # In Intron
        elif zone == 2:
            if strandDirection == "+":
                if d1 < d2:
                    desc = (
                        "Intron of "
                        + geneName
                        + "("
                        + strandDirection
                        + "):"
                        + hp.bp2str(d1, 2)
                        + " after exon "
                        + str(e1)
                    )
                else:
                    desc = (
                        "Intron of "
                        + geneName
                        + "("
                        + strandDirection
                        + "):"
                        + hp.bp2str(d2, 2)
                        + " before exon "
                        + str(e2)
                    )
                descList.append(desc)
                intronnum = e1
                if intronnum:
                    intronnumList.append(intronnum)
                else:
                    intronnumList.append("Null")
                intronframe = f
                if intronframe:
                    intronframeList.append(intronframe)
                else:
                    intronframeList.append("Null")
                # print desc
                continue
            else:
                if d1 < d2:
                    desc = (
                        "Intron of "
                        + geneName
                        + "("
                        + strandDirection
                        + "):"
                        + hp.bp2str(d1, 2)
                        + " before exon "
                        + str(e1)
                    )
                else:
                    desc = (
                        "Intron of "
                        + geneName
                        + "("
                        + strandDirection
                        + "):"
                        + hp.bp2str(d2, 2)
                        + " after exon "
                        + str(e2)
                    )
                descList.append(desc)
                intronnum = e1
                if intronnum:
                    intronnumList.append(intronnum)
                else:
                    intronnumList.append("Null")
                intronframe = f
                if intronframe:
                    intronframeList.append(intronframe)
                else:
                    intronframeList.append("Null")
                # print desc
                continue

        # In 3'-UTR
        elif zone == 3:
            desc = (
                "3'-UTR of "
                + geneName
                + "("
                + strandDirection
                + "):"
                + hp.bp2str(d, 2)
                + " after coding stop"
            )
            descList.append(desc)
            intronnumList.append(intronnum)
            intronframeList.append(intronframe)
            # print desc
            continue
        # In 5'-UTR
        elif zone == 4:
            desc = (
                "5'-UTR of "
                + geneName
                + "("
                + strandDirection
                + "):"
                + hp.bp2str(d, 2)
                + " before coding start"
            )
            descList.append(desc)
            intronnumList.append(intronnum)
            intronframeList.append(intronframe)
            # print desc
            continue
        # In Promoter
        elif zone == 5:
            desc = (
                "Promoter of "
                + geneName
                + "("
                + strandDirection
                + "):"
                + hp.bp2str(d, 2)
                + " from tx start"
            )
            descList.append(desc)
            intronnumList.append(intronnum)
            intronframeList.append(intronframe)
            # print desc
            continue
        else:
            desc = "Unexpected Error"
            descList.append(desc)
            intronnumList.append(intronnum)
            intronframeList.append(intronframe)
            # print desc
            continue
    return (
        geneNameList,
        transcriptList,
        descList,
        zoneList,
        strandDirectionList,
        intronnumList,
        intronframeList,
    )
