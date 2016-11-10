'''
Created on 12/23/2015
@Ronak Shah

'''

from collections import defaultdict
import pandas as pd
import logging
import coloredlogs
# Gives elements at particular index in list
getVar = lambda searchList, ind: [searchList[i] for i in ind]
coloredlogs.install(level='DEBUG')
def ReadRepeatFile(filename, verbose):
    if(verbose):
        logging.info("iAnnotateSV::AnnotateForRepeatRegion: Reading & Storing Repeat TSV file as dictionary")
    # Initialize dictionary of lists 
    dataDict = defaultdict(list)
    with open(filename, 'r') as filecontent:
        header = filecontent.readline()
        for line in filecontent:
            data = line.rstrip('\n').split('\t')
            processedData = (data[0].replace('chr', ''))
            slicedData = data[1:]
            joinedData = '\t'.join(slicedData)
            dataDict[processedData].append(joinedData)    
    return dataDict        

def AnnotateRepeatRegion (verbose, count, sv, rrDict):
    if(verbose):
        logging.info("iAnnotateSV::AnnotateForRepeatRegion: Checking Entry %d in Repeat data", count)
    # Initialize List to store repeat annotation
    list_svloc1 = []
    list_svloc2 = []
    # Read SV Data
    sv_chr1 = str(sv.loc['chr1'])
    sv_pos1 = int(sv.loc['pos1'])
    sv_chr2 = str(sv.loc['chr2'])
    sv_pos2 = int(sv.loc['pos2'])
    # Traverse through Repeat Data Dict
    list_loc1 = rrDict.get(sv_chr1, "None")  # Get the values for the chromosome
    if(list_loc1 != "None"):  # Check if there are no keys with a particular chromosome
        for loc in list_loc1:  # For each location in all values check the overlap
            data = loc.split('\t') 
            rr_pos1 = int(data[0])
            rr_pos2 = int(data[1])
            if (rr_pos1 <= sv_pos1 <= rr_pos2):
                slicedData = data[2:]
                joinedData = '-'.join(slicedData)
                list_svloc1.append(joinedData)
    else:
        if(verbose):
            logging.info("iAnnotateSV::AnnotateForRepeatRegion: Chromosome %s is not there in the repeat dictionary", sv_chr1)
    list_loc2 = rrDict.get(sv_chr2, "None")
    if(list_loc2 != "None"):
        for loc in list_loc2:
            data = loc.split('\t') 
            rr_pos1 = int(data[0])
            rr_pos2 = int(data[1])
            if (rr_pos1 <= sv_pos2 <= rr_pos2):
                slicedData = data[2:]
                joinedData = '-'.join(slicedData)
                list_svloc2.append(joinedData)
    else:
        if(verbose):
           logging.info("iAnnotateSV::AnnotateForRepeatRegion: Chromosome %s is not there in the repeat dictionary", sv_chr2) 
    return (list_svloc1, list_svloc2)  