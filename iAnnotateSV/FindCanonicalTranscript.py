'''
Created on  Mar 4, 2015

@author: Ronak H Shah
'''
from operator import itemgetter
import logging
import coloredlogs
'''
Get the mininmumIndex for a list and return the variables that match canonical transcript
Preference:# zone: 1=exon, 2=intron, 3=3'-UTR, 4=5'-UTR, 5=promoter
'''  
coloredlogs.install(level='DEBUG')
def FindCT(geneList,transcriptList,siteList,zoneList,strandList,intronnumList,intronframeList,ctDict):
    gene = None
    transcript = None
    site = None
    zone = None
    strand = None
    intronnum = None
    intronframe = None
    chkval = isinstance(geneList,list)
    # print geneList,transcriptList
    if(chkval):
        cts = None
        minIndex = None
        for gene in geneList:
            #print gene
            if gene in ctDict:
                cts = ctDict.get(gene)
                break
        if(cts):
            if(len(cts) > 1 ):
                minIndex = min(enumerate(zoneList), key=itemgetter(1))[0] 
            else:
                #print "CTS",cts[0]
                try:
                    #print "I am here", transcriptList.index(cts[0])
                    minIndex = transcriptList.index(cts[0])
                except ValueError:
                    logging.warn("iAnnotateSV::FindCanonicalTranscript: The given canonical transcript does not cover the coordinates.")
                    minIndex = min(enumerate(zoneList), key=itemgetter(1))[0]
        else:
            minIndex = min(enumerate(zoneList), key=itemgetter(1))[0] 
        #print minIndex
        gene = geneList[minIndex]
        transcript = transcriptList[minIndex]
        site = siteList[minIndex] 
        zone = zoneList[minIndex]
        strand = strandList[minIndex]
        intronnum = intronnumList[minIndex] 
        if(intronnum == "Null"):
            intronnum = None
        intronframe = intronframeList[minIndex]  
        if(intronframe == "Null"):
            intronframe = None
    else:
        gene = geneList
        transcript = transcriptList
        site = siteList
        zone = zoneList
        strand = strandList
        intronnum = intronnumList
        intronframe = intronframeList
    return(gene,transcript,site,zone,strand,intronnum,intronframe)