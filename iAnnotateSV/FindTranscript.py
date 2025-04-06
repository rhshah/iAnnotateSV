'''
Created on Mar 4, 2015

@author: Ronak H Shah
'''

import helper as hp
import polars as pl
from rich import print

'''
This function will return a single highest priority transcript
'''


def FindATranscript(queryDF, refDF):
    '''
    This function identifies the highest priority transcript from a set of transcripts
    overlapping a given genomic region.

    It selects the transcript based on criteria such as minimum distance to a feature
    and returns details about the selected transcript, including its gene name,
    transcript ID, genomic location description, zone, strand direction, intron
    number (if applicable), and intron frame (if applicable).

    Args:
        queryDF (pl.DataFrame): DataFrame containing information about the genomic region.
        refDF (pl.DataFrame): DataFrame containing reference transcript annotations.

    Returns:
        tuple: A tuple containing the following information about the selected transcript:
            - geneName (str): Name of the gene.
            - transcript (str): Transcript ID.
            - desc (str): Description of the transcript's genomic location.
            - zone (int): Genomic zone of the transcript.
            - strandDirection (str): Strand direction ('+' or '-').
            - intronnum (int or None): Intron number if the transcript is intronic, otherwise None.
            - intronframe (int or None): Intron frame if the transcript is intronic, otherwise None.
    '''
    tmpDF = queryDF.clone()  # Use clone instead of copy for Polars

    desc = None
    intronnum = None
    intronframe = None

    if tmpDF.height > 1:  # Use height instead of len(queryDF.index)
        try:
            tmpDF = tmpDF.with_columns(
                pl.col("d").cast(pl.Int64).alias("d")
            )  # Convert 'd' to numeric
            idxMin = tmpDF.select(pl.col("d").min().alias("min_d"))[0, "min_d"]
            idxMin = tmpDF.filter(pl.col("d") == idxMin).row(0)[0]  # Get index of min 'd'
        except ValueError:
            idxMin = tmpDF[0, 0]  # First index if ValueError
    else:
        idxMin = 0  # First index if only one row

    row = tmpDF[idxMin]  # Get the row as a Series

    c = row["c"]

    d = row["d"]
    d1 = row["d1"]
    d2 = row["d2"]
    e = row["e"]
    e1 = row["e1"]
    e2 = row["e2"]
    f = row["c"]
    zone = c
    transcript = refDF[idxMin, "#name"]
    geneName = refDF[idxMin, "name2"]
    strandDirection = refDF[idxMin, "strand"]

    if (strandDirection == '-'):
        if(e):
            e = int(refDF[idxMin, 'exonCount']) - e + 1
        else:
            e = int(refDF[idxMin, 'exonCount']) + 1
        if(e1):
            e1 = int(refDF[idxMin, 'exonCount']) - e1 + 1
        else:
            e1 = int(refDF[idxMin, 'exonCount']) + 1
        if(e2):
            e2 = int(refDF[idxMin, 'exonCount']) - e2 + 1
        else:
            e2 = int(refDF[idxMin, 'exonCount']) + 1
    # in Exon:
    if(c == 1):
        desc = 'Exon ' + str(e) + " of " + geneName + \
            '(' + strandDirection + ')'
        # print(f"[green]{desc}[/green]")

    # In Intron
    elif(c == 2):
        if(strandDirection == "+"):
            if(d1 < d2):
                desc = 'Intron of ' + geneName + \
                    '(' + strandDirection + '):' + \
                    hp.bp2str(d1, 2) + ' after exon ' + str(e1)
            else:
                desc = 'Intron of ' + geneName + \
                    '(' + strandDirection + '):' + \
                    hp.bp2str(d2, 2) + ' before exon ' + str(e2)
            # print(f"[green]{desc}[/green]")

        else:
            if(d1 < d2):
                desc = 'Intron of ' + geneName + \
                    '(' + strandDirection + '):' + \
                    hp.bp2str(d1, 2) + ' before exon ' + str(e1)
            else:
                desc = 'Intron of ' + geneName + \
                    '(' + strandDirection + '):' + \
                    hp.bp2str(d2, 2) + ' after exon ' + str(e2)
            # print(f"[green]{desc}[/green]")

        intronnum = e1
        intronframe = f
    # In 3'-UTR
    elif(c == 3):
        desc = '3\'-UTR of ' + geneName + \
            '(' + strandDirection + '):' + \
            hp.bp2str(d, 2) + ' after coding stop'
        # print(f"[green]{desc}[/green]")

    # In 5'-UTR
    elif(c == 4):
        desc = '5\'-UTR of ' + geneName + \
            '(' + strandDirection + '):' + \
            hp.bp2str(d, 2) + ' before coding start'
        # print(f"[green]{desc}[/green]")

    # In Promoter
    elif(c == 5):
        desc = 'Promoter of ' + geneName + \
            '(' + strandDirection + '):' + hp.bp2str(d, 2) + ' from tx start'
        # print(f"[green]{desc}[/green]")

    else:
        desc = 'Unexpected Error'
        # print(f"[red]{desc}[/red]")
    return(geneName, transcript, desc, zone, strandDirection, intronnum, intronframe)


'''
This function will return all the transcripts
'''


def FindAllTranscripts(queryDF, refDF):
    '''
    This function identifies all transcripts overlapping a given genomic region.

    It iterates through each transcript and determines its genomic location
    description, zone, strand direction, intron number (if applicable), and intron
    frame (if applicable).

    Args:
        queryDF (pl.DataFrame): DataFrame containing information about the genomic region.
        refDF (pl.DataFrame): DataFrame containing reference transcript annotations.

    Returns:
        tuple: A tuple containing lists of information about all selected transcripts:
            - geneNameList (list): List of gene names.
            - transcriptList (list): List of transcript IDs.
            - descList (list): List of descriptions of transcript genomic locations.
            - zoneList (list): List of genomic zones.
            - strandDirectionList (list): List of strand directions.
            - intronnumList (list): List of intron numbers (or "Null").
            - intronframeList (list): List of intron frames (or "Null").
    '''
    geneNameList = []
    transcriptList = []
    strandDirectionList = []
    descList = []
    zoneList = []
    intronnumList = []
    intronframeList = []

    for row in queryDF.iter_rows(named=True):  # Iterate over rows as dictionaries
        desc = "Null"
        intronnum = "Null"
        intronframe = "Null"
        transcript = refDF[row, '#name']
        transcriptList.append(transcript)
        geneName = refDF[row, 'name2']
        geneNameList.append(geneName)
        strandDirection = refDF[row, 'strand']
        strandDirectionList.append(strandDirection)
        # print(f"[blue]Transcript: {transcript}, Gene: {geneName}, Strand: {strandDirection}[/blue]")
        c = int(row['c'])
        if(row['d']):
            d = int(row['d'])
        else:
            d1 = row['d']
        e = row['e']
        if(row['d1']):
            d1 = int(row['d1'])
        else:
            d1 = row['d1']
        if(row['d2']):
            d2 = int(row['d2'])
        else:
            d1 = row['d2']
        e1 = row['e1']
        e2 = row['e2']
        f = row['f']
        zone = c
        zoneList.append(zone)
        if (strandDirection == '-'):
            if(e):
                e = int(refDF[row]['exonCount']) - int(row['e']) + 1
            else:
                e = int(refDF[row]['exonCount']) + 1
            if(e1):
                e1 = int(refDF[row]['exonCount']) - \
                    int(row['e1']) + 1
            else:
                e1 = int(refDF[row]['exonCount']) + 1
            if(e2):
                e2 = int(refDF[row]['exonCount']) - \
                    int(row['e2']) + 1
            else:
                e2 = int(refDF[row]['exonCount']) + 1
        # in Exon:
        if(zone == 1):
            desc = 'Exon ' + str(e) + " of " + geneName + \
                '(' + strandDirection + ')'
            descList.append(desc)
            intronnumList.append(intronnum)
            intronframeList.append(intronframe)
            # print(f"[green]{desc}[/green]")
            continue
        # In Intron
        elif(zone == 2):
            if(strandDirection == "+"):
                if(d1 < d2):
                    desc = 'Intron of ' + geneName + \
                        '(' + strandDirection + '):' + \
                        hp.bp2str(d1, 2) + ' after exon ' + str(e1)
                else:
                    desc = 'Intron of ' + geneName + \
                        '(' + strandDirection + '):' + \
                        hp.bp2str(d2, 2) + ' before exon ' + str(e2)
                descList.append(desc)
                intronnum = e1
                if(intronnum):
                    intronnumList.append(intronnum)
                else:
                    intronnumList.append("Null")
                intronframe = f
                if(intronframe):
                    intronframeList.append(intronframe)
                else:
                    intronframeList.append("Null")
                # print(f"[green]{desc}[/green]")
                continue
            else:
                if(d1 < d2):
                    desc = 'Intron of ' + geneName + \
                        '(' + strandDirection + '):' + \
                        hp.bp2str(d1, 2) + ' before exon ' + str(e1)
                else:
                    desc = 'Intron of ' + geneName + \
                        '(' + strandDirection + '):' + \
                        hp.bp2str(d2, 2) + ' after exon ' + str(e2)
                descList.append(desc)
                intronnum = e1
                if(intronnum):
                    intronnumList.append(intronnum)
                else:
                    intronnumList.append("Null")
                intronframe = f
                if(intronframe):
                    intronframeList.append(intronframe)
                else:
                    intronframeList.append("Null")
                # print(f"[green]{desc}[/green]")
                continue

        # In 3'-UTR
        elif(zone == 3):
            desc = '3\'-UTR of ' + geneName + \
                '(' + strandDirection + '):' + \
                hp.bp2str(d, 2) + ' after coding stop'
            descList.append(desc)
            intronnumList.append(intronnum)
            intronframeList.append(intronframe)
            # print(f"[green]{desc}[/green]")
            continue
        # In 5'-UTR
        elif(zone == 4):
            desc = '5\'-UTR of ' + geneName + \
                '(' + strandDirection + '):' + \
                hp.bp2str(d, 2) + ' before coding start'
            descList.append(desc)
            intronnumList.append(intronnum)
            intronframeList.append(intronframe)
            # print(f"[green]{desc}[/green]")
            continue
        # In Promoter
        elif(zone == 5):
            desc = 'Promoter of ' + geneName + \
                '(' + strandDirection + '):' + \
                hp.bp2str(d, 2) + ' from tx start'
            descList.append(desc)
            intronnumList.append(intronnum)
            intronframeList.append(intronframe)
            # print(f"[green]{desc}[/green]")
            continue
        else:
            desc = 'Unexpected Error'
            descList.append(desc)
            intronnumList.append(intronnum)
            intronframeList.append(intronframe)
            # print(f"[red]{desc}[/red]")
            continue
    return(geneNameList, transcriptList, descList, zoneList, strandDirectionList, intronnumList, intronframeList)