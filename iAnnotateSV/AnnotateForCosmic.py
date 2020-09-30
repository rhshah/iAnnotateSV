'''
Created on 12/23/2015
@Ronak Shah

'''
import pandas as pd
import logging
import coloredlogs
# Gives elements at particular index in list
getVar = lambda searchList, ind: [searchList[i] for i in ind]
coloredlogs.install(level='DEBUG')
def AnnotateFromCosmicCensusFile (filename, verbose, count, sv):
    if(verbose):
        logging.info("iAnnotateSV::AnnotateForCosmic: Checking Entry %d in Cosmic", count)
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

def AnnotateFromCosmicFusionCountsFile (filename, verbose, count, sv):
    if(verbose):
        logging.info("iAnnotateSV::AnnotateForCosmic: Checking Entry %d in Cosmic Counts data", count)
    # Initialize List to store comic annotations
    sv_gene1 = str(sv.loc['gene1'])
    sv_gene2 = str(sv.loc['gene2'])
    sv_combo1 = sv_gene1 + "-" + sv_gene2
    sv_combo2 = sv_gene2 + "-" + sv_gene1
 
    countDF = pd.read_csv(filename, sep='\t', header=0, keep_default_na='True')
    for index,row in countDF.iterrows():
        gene1 = str(row.loc["Gene1"])
        gene2 = str(row.loc["Gene2"])
        combo1 = gene1 + "-" + gene2
        combo2 = gene2 + "-" + gene1
        counts = None
        if(sv_combo1 == combo1 or sv_combo1 == combo2 or sv_combo2 == combo1 or sv_combo2 == combo2 ):
            counts = int(row.loc['Counts'])
            break
        else:
            continue
    return counts