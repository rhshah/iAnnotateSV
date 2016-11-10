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
def ReadDGvFile(filename, verbose):
    if(verbose):
        logging.info("iAnnotateSV::AnnotateForDGv: Reading & Storing DGV TSV file as dictionary")
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
         
def AnnotateDGv (verbose, count, sv, dgvDict):
    if(verbose):
        logging.info("iAnnotateSV::AnnotateForDGv: Checking Entry %d in DGv data", count)
    # Initialize List to store repeat annotation
    list_svloc1 = []
    list_svloc2 = []
    # Read SV Data
    sv_chr1 = str(sv.loc['chr1'])
    sv_pos1 = int(sv.loc['pos1'])
    sv_chr2 = str(sv.loc['chr2'])
    sv_pos2 = int(sv.loc['pos2'])
    # Traverse through DGv Data Dict
    list_loc1 = dgvDict.get(sv_chr1, "None")  # Get the values for the chromosome
    if(list_loc1 != "None"):  # Check if there are no keys with a particular chromosome
        for loc in list_loc1:  # For each location in all values check the overlap
            data = loc.split('\t') 
            dgv_pos1 = int(data[0])
            dgv_pos2 = int(data[1])
            if (dgv_pos1 <= sv_pos1 <= dgv_pos2):
                slicedData = getVar(data, [2, 8])
                joinedData = '-'.join(slicedData)
                list_svloc1.append(joinedData)
    else:
        if(verbose):
            logging.warn("iAnnotateSV::AnnotateForDGv: Chromosome %s is not there in the DGv dictionary", sv_chr1)        
    list_loc2 = dgvDict.get(sv_chr2, "None")
    if(list_loc2 != "None"):
        for loc in list_loc2:
            data = loc.split('\t') 
            dgv_pos1 = int(data[0])
            dgv_pos2 = int(data[1])
            if (dgv_pos1 <= sv_pos2 <= dgv_pos2):
                slicedData = getVar(data, [2, 8])
                joinedData = '-'.join(slicedData)
                list_svloc2.append(joinedData)
    else:
        if(verbose):
            logging.warn("iAnnotateSV::AnnotateForDGv: Chromosome %s is not there in the DGv dictionary", sv_chr2)    
    return (list_svloc1, list_svloc2)