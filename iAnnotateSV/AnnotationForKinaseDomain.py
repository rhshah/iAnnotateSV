'''
Created on 01/09/2018
@Ronak Shah

'''
import os
import sys
import pandas as pd
import logging
import coloredlogs
import re
import helper as hp
coloredlogs.install(level='DEBUG')


def run(svDFA, refPath, ctPath, allctPath, upPath, verbose):
    if(os.path.isfile(upPath)):
        upDF = hp.ReadFile(upPath)
    else:
        if(verbose):
            logging.critical(
                "iAnnotateSV::AnnotationForKinaseDomain: Location of Uniprot Annoation file is incorrect!!!")
        sys.exit(1)
    if(os.path.isfile(ctPath)):
        ctDF = hp.ReadFile(ctPath)
    else:
        if(verbose):
            logging.warn(
                "iAnnotateSV::AnnotationForKinaseDomain: Location of assay specific canonical transcript file is incorrect!!!")
        ctDF = pd.DataFrame()
    if(os.path.isfile(allctPath)):
        allctDF = hp.ReadFile(allctPath)
    else:
        if(verbose):
            logging.critical(
                "iAnnotateSV::AnnotationForKinaseDomain: Location of all canonical transcript file is incorrect!!!")
        sys.exit(1)
    if(os.path.isfile(refPath)):
        refDF = hp.ReadFile(refPath)
        refDF.columns = refDF.columns.str.replace('#', '')
    else:
        if(verbose):
            logging.critical(
                "iAnnotateSV::AnnotationForKinaseDomain: Location of reference based annotation file is incorrect!!!")
        sys.exit(1)
    svDF = svDFA.copy()
    svDF.insert(loc=9, column='kinase_domain1', value=None)
    svDF.insert(loc=13, column='kinase_domain2', value=None)
    #svDF["kinase_domain1"] = None
    #svDF["kinase_domain2"] = None
    for count, row in svDFA.iterrows():
        # print row
        if(verbose):
            logging.info(
                "iAnnotateSV::AnnotateForKinaseDomain: Checking Entry %d in Uniprot data", count)
        chr1 = str(row.loc['chr1'])
        chr2 = str(row.loc['chr2'])
        if(chr1.startswith('chr')):
            chr1 = chr1
        else:
            chr1 = "chr" + chr1
        if(chr2.startswith('chr')):
            chr2 = chr2
        else:
            chr2 = "chr" + chr2
        pos1 = int(row.loc['pos1'])
        pos2 = int(row.loc['pos2'])
        gene1 = str(row.loc['gene1'])
        gene2 = str(row.loc['gene2'])
        site1 = str(row.loc['site1'])
        site2 = str(row.loc['site2'])

        try:
            transcript1 = ctDF.Transcripts[ctDF.Gene[ctDF.Gene == gene1].index.tolist()[
                0]]
        except IndexError:
            try:
                transcript1 = allctDF.Transcripts[allctDF.Gene[allctDF.Gene == gene1].index.tolist()[
                    0]]
            except IndexError:
                transcript1 = None
        try:
            transcript2 = ctDF.Transcripts[ctDF.Gene[ctDF.Gene == gene2].index.tolist()[
                0]]
        except IndexError:
            try:
                transcript2 = allctDF.Transcripts[allctDF.Gene[allctDF.Gene == gene2].index.tolist()[
                    0]]
            except IndexError:
                transcript2 = None
        fusion = str(row.loc['fusion'])

        kanno1 = None
        kanno2 = None

        if(fusion != "-"):
            # First Gene +, Second Gene -
            fusionevent = re.search(r'\{(.*)\}', fusion)
            if(fusionevent):
                eventType = fusionevent.group(1)
                if(":" in eventType):
                    # print fusion, fusionevent, eventType
                    (egene1, egene2) = (str(eventType)).split(":")

                    if(transcript1):
                        kanno1 = getKinaseInfo(
                            chr1, pos1, gene1, egene1, egene2, transcript1, refDF, upDF)
                    else:
                        kanno1 = None

                    if(transcript2):
                        kanno2 = getKinaseInfo(
                            chr2, pos2, gene2, egene1, egene2, transcript2, refDF, upDF)
                    else:
                        kanno2 = None
                else:
                    kanno1 = None
                    kanno2 = None
            else:
                kanno1 = None
                kanno2 = None
            svDF.loc[count, 'kinase_domain1'] = kanno1
            svDF.loc[count, 'kinase_domain2'] = kanno2

    return(svDF)


def processData(chrom, transcript, refDF, upDF):
    transcripts = (refDF[refDF['name'] == transcript])
    if (len(transcripts) > 1):
        transcriptIdx, = (transcripts[transcripts['chrom'] == chrom].index)
    else:
        try:
            transcriptIdx, = (refDF[refDF['name'] == transcript].index)
        except ValueError:
            return (None, None, None)

    refTxSt = int(refDF.iloc[transcriptIdx]['txStart'])
    refTxEn = int(refDF.iloc[transcriptIdx]['txEnd'])
    # print "1:",transcriptIdx,"\n",refTxSt,"\n", refTxEn, "\n"
    up_idxList = upDF[upDF['#chrom'] == chrom].index.tolist()
    # Find all overlapping transcripts
    up_recordIndex = []
    for index in (up_idxList):
        # print upDF.iloc[index],"\n"
        chromStart = upDF.iloc[index]['chromStart']
        chromEnd = upDF.iloc[index]['chromEnd']
        if ((chromStart >= refTxSt) and (chromEnd <= refTxEn)):
            # print "Chr" , chromStart,chromEnd, refTxSt, refTxEn,"\n"
            if (upDF.iloc[index]['annotationType'] == 'domain'):
                up_recordIndex.append(index)
    allMaxVal = []
    allMinVal = []
    for index, val in enumerate(up_recordIndex):
        chromStart = upDF.iloc[val]['chromStart']
        chromEnd = upDF.iloc[val]['chromEnd']
        maxVal = max(refTxEn, chromEnd)
        allMaxVal.append(maxVal)
        minVal = min(refTxSt, chromStart)
        allMinVal.append(minVal)
    if (allMaxVal):
        max_len = max(allMaxVal)
    else:
        max_len = refTxEn
    if (allMinVal):
        min_len = max(allMinVal)
    else:
        min_len = refTxSt
    return (up_recordIndex, max_len, min_len)


def getKinaseInfo(chrom, pos, gene, egene1, egene2, transcript, refDF, upDF):
    (domainIdx, maxLen, minLen) = processData(chrom, transcript, refDF, upDF)
    if(domainIdx is None):
        return None
    strand = refDF.strand[refDF.name[refDF.name == transcript].index.tolist()[
        0]]
    #kanno = None
    if(strand == "+"):
        if(egene1 == gene):
            # print "Here1"
            # See if Kinase occurs after the breakpoint or within the breakpoint
            for index, val in enumerate(domainIdx):
                chromStart = upDF.iloc[val]['chromStart']
                chromEnd = upDF.iloc[val]['chromEnd']
                fname = upDF.iloc[val]['name']
                if("Protein kinase" in fname):
                    if (pos > chromEnd):
                        kanno = "Kinase Domain Included"
                    else:
                        if(chromStart <= pos):
                            if(pos <= chromEnd):
                                kanno = "Partial Kinase Domain Included"
                            else:
                                kanno = "Kinase Domain Not Included"
                        else:
                            if(chromEnd <= pos):
                                if(pos <= chromStart):
                                    kanno = "Partial Kinase Domain Included"
                                else:
                                    kanno = "Kinase Domain Not Included"
                            else:
                                kanno = "Kinase Domain Not Included"
                    # print gene, pos, chromStart, chromEnd, transcript, strand, kanno
                    return(kanno)

        if(egene2 == gene):
            # print "Here2"
            # See if Kinase occurs after the breakpoint or within the breakpoint
            for index, val in enumerate(domainIdx):
                chromStart = upDF.iloc[val]['chromStart']
                chromEnd = upDF.iloc[val]['chromEnd']
                fname = upDF.iloc[val]['name']
                if("Protein kinase" in fname):
                    if(pos < chromStart):
                        kanno = "Kinase Domain Included"
                    else:
                        if(chromStart <= pos):
                            if(pos <= chromEnd):
                                kanno = "Partial Kinase Domain Included"
                            else:
                                kanno = "Kinase Domain Not Included"
                        else:
                            if(chromEnd <= pos):
                                if(pos <= chromStart):
                                    kanno = "Partial Kinase Domain Included"
                                else:
                                    kanno = "Kinase Domain Not Included"
                            else:
                                kanno = "Kinase Domain Not Included"
                    # print gene, pos, chromStart, chromEnd, transcript, strand, kanno
                    return(kanno)
    else:
        if(egene1 == gene):
            # print "Here3"
            # See if Kinase occurs after the breakpoint or within the breakpoint
            for index, val in enumerate(domainIdx):
                chromStart = upDF.iloc[val]['chromStart']
                chromEnd = upDF.iloc[val]['chromEnd']
                fname = upDF.iloc[val]['name']
                if ("Protein kinase" in fname):
                    if(pos < chromStart):
                        kanno = "Kinase Domain Included"
                    else:
                        if(chromStart <= pos):
                            if(pos <= chromEnd):
                                kanno = "Partial Kinase Domain Included"
                            else:
                                kanno = "Kinase Domain Not Included"
                        else:
                            if(chromEnd <= pos):
                                if(pos <= chromStart):
                                    kanno = "Partial Kinase Domain Included"
                                else:
                                    kanno = "Kinase Domain Not Included"
                            else:
                                kanno = "Kinase Domain Not Included"
                    # print gene, pos, chromStart, chromEnd, transcript, strand, kanno
                    return(kanno)

        if(egene2 == gene):
            # print "Here4"
            # See if Kinase occurs after the breakpoint or within the breakpoint
            for index, val in enumerate(domainIdx):
                chromStart = upDF.iloc[val]['chromStart']
                chromEnd = upDF.iloc[val]['chromEnd']
                fname = upDF.iloc[val]['name']
                if("Protein kinase" in fname):
                    if(pos > chromEnd):
                        kanno = "Kinase Domain Included"
                    else:
                        if(chromStart <= pos):
                            if(pos <= chromEnd):
                                kanno = "Partial Kinase Domain Included"
                            else:
                                kanno = "Kinase Domain Not Included"
                        else:
                            if(chromEnd <= pos):
                                if(pos <= chromStart):
                                    kanno = "Partial Kinase Domain Included"
                                else:
                                    kanno = "Kinase Domain Not Included"
                            else:
                                kanno = "Kinase Domain Not Included"
                    # print gene, pos, chromStart, chromEnd, transcript, strand, kanno
                    return(kanno)
