'''
Created on 12/23/2015
@Ronak Shah

'''
import polars as pl
import logging
from rich import print
from rich.logging import RichHandler

# Gives elements at particular index in list
getVar = lambda searchList, ind: [searchList[i] for i in ind]

FORMAT = "%(message)s"
logging.basicConfig(
    level="INFO", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)

log = logging.getLogger("rich")


def AnnotateFromCosmicCensusFile(filename, verbose, count, sv):
    if(verbose):
        log.info("iAnnotateSV::AnnotateForCosmic: Checking Entry %d in Cosmic", count)
    # Initialize List to store comic annotations
    list_ccData = []
    sv_gene1 = str(sv['gene1'])
    sv_gene2 = str(sv['gene2'])

    with open(filename, 'r') as filecontent:
        header = filecontent.readline()
        for line in filecontent:
            data = line.rstrip('\n').split('\t')
            if (str(data[0]) == sv_gene1):
                slicedData = getVar(data, [4, 7, 9, 12, 13])
                slicedProcessedData = []
                for sData in slicedData:
                    if sData:
                        sData = f"site1:{sData}"
                        slicedProcessedData.append(sData)
                    else:
                        slicedProcessedData.append(" ")
                joinedData = '\t'.join(slicedProcessedData)
                list_ccData.append(joinedData)
            if (str(data[0]) == sv_gene2):
                slicedData = getVar(data, [4, 7, 9, 12, 13])
                slicedProcessedData = []
                for sData in slicedData:
                    if sData:
                        sData = f"site2:{sData}"
                        slicedProcessedData.append(sData)
                    else:
                        slicedProcessedData.append(" ")
                joinedData = '\t'.join(slicedProcessedData)
                list_ccData.append(joinedData)
    return list_ccData


def AnnotateFromCosmicFusionCountsFile(filename, verbose, count, sv):
    if(verbose):
        log.info("iAnnotateSV::AnnotateForCosmic: Checking Entry %d in Cosmic Counts data", count)
    # Initialize List to store comic annotations
    sv_gene1 = str(sv['gene1'])
    sv_gene2 = str(sv['gene2'])
    sv_combo1 = f"{sv_gene1}-{sv_gene2}"
    sv_combo2 = f"{sv_gene2}-{sv_gene1}"

    countDF = pl.read_csv(filename, sep='\t', infer_schema_length=1000)
    counts = None
    for row in countDF.iter_rows(named=True):
        gene1 = str(row["Gene1"])
        gene2 = str(row["Gene2"])
        combo1 = f"{gene1}-{gene2}"
        combo2 = f"{gene2}-{gene1}"

        if(sv_combo1 == combo1 or sv_combo1 == combo2 or sv_combo2 == combo1 or sv_combo2 == combo2):
            counts = int(row['Counts'])
            break
        else:
            continue
    return counts