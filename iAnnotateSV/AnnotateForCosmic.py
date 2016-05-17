'''
Created on 12/23/2015
@Ronak Shah

'''
import pandas as pd
import logging
# Gives elements at particular index in list
getVar = lambda searchList, ind: [searchList[i] for i in ind]

def AnnotateFromCosmicCensusFile (filename, verbose, count, sv):
    if(verbose):
        logging.info("iAnnotateSV::AnnotateForCosmic: Checking Entry in Cosmic data for entry %d", count)
    # Initialize List to store comic annotations
    list_ccData = []
    sv_gene1 = str(sv.loc['gene1'])
    sv_gene2 = str(sv.loc['gene2'])
 
    with open(filename, 'r') as filecontent:
        header = filecontent.readline()
        for line in filecontent:
            data = line.rstrip('\n').split('\t')
            if(str(data[0]) == sv_gene1): 
                slicedData = getVar(data,[4,7,9,12,13])
                slicedProcessedData = []
                for sData in slicedData:
                    if(sData):
                        sData = "site1:" + sData
                        slicedProcessedData.append(sData)
                    else:
                        slicedProcessedData.append(" ")
                joinedData = '\t'.join(slicedProcessedData)
                list_ccData.append(joinedData)  
            if(str(data[0]) == sv_gene2):
                slicedData = getVar(data,[4,7,9,12,13]) 
                slicedProcessedData = []
                for sData in slicedData:
                    if(sData):
                        sData = "site2:" + sData   
                        slicedProcessedData.append(sData)
                    else:
                        slicedProcessedData.append(" ")
                joinedData = '\t'.join(slicedProcessedData)
                list_ccData.append(joinedData)             
    return list_ccData