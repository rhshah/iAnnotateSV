"""
Created on 12/23/2015.

@Ronak Shah

"""

import sys
import argparse
import pandas as pd
import time
import logging
import AnnotateForRepeatRegion as afr
import AnnotateForCosmic as afc
import AnnotateForDGv as afd
import coloredlogs

def main(command=None):
    parser = argparse.ArgumentParser(
        prog='AddExternalAnnotations.py',
        description='Add External Annotation to the Structural Variants',
        usage='%(prog)s [options]')
    parser.add_argument(
        "-r",
        "--repeatFile",
        action="store",
        dest="rrFilename",
        required=True,
        metavar='RepeatRegionFile.tsv',
        help="Location of the Repeat Region Bed File")
    parser.add_argument(
        "-d",
        "--dgvFile",
        action="store",
        dest="dgvFilename",
        required=True,
        metavar='DGvFile.tsv',
        help="Location of the Database of Genomic Variants Bed File")
    parser.add_argument(
        "-c",
        "--cosmicConsensusFile",
        action="store",
        dest="ccFilename",
        required=True,
        metavar='CosmicConsensus.tsv',
        help="Location of the Cosmic Consensus TSV file")
    parser.add_argument(
        "-cct",
        "--cosmicCountsFile",
        action="store",
        dest="cctFilename",
        required=True,
        metavar='cosmic_fusion_counts.tsv',
        help="Location of the Cosmic Counts TSV file")
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        default=True,
        help="make lots of noise [default]")
    parser.add_argument(
        "-s",
        "--svFile",
        action="store",
        dest="svFilename",
        required=True,
        metavar='SVfile.txt',
        help="Location of the structural variant file to be annotated")
    parser.add_argument(
        "-ofp",
        "--outputFilePrefix",
        action="store",
        dest="outFilePrefix",
        required=True,
        metavar='AnnotatedSV',
        help="Full path with prefix name for the output file")
    parser.add_argument(
        "-o",
        "--outputDir",
        action="store",
        dest="outDir",
        required=True,
        metavar='/somedir',
        help="Full Path to the output dir")
    args = ''
    if(command is None):
        args = parser.parse_args()
    else:
        args = parser.parse_args(command.split())
    coloredlogs.install(level='DEBUG')
    outFileTxt = args.outDir + "/" + args.outFilePrefix + ".txt"
    outFileExl = args.outDir + "/" + args.outFilePrefix + ".xlsx"
    outFileJson = args.outDir + "/" + args.outFilePrefix + ".json"
    if args.verbose:
        logging.info("iAnnotateSV::AddExternalAnnotations: Reading %s...", args.svFilename)
        data = ReadSVFile(args)
        logging.info("iAnnotateSV::AddExternalAnnotations: Finished Reading %s", args.svFilename)
        logging.info("iAnnotateSV::AddExternalAnnotations: Reading %s...", args.rrFilename)
        repeatRegionDict = afr.ReadRepeatFile(args.rrFilename, args.verbose)
        logging.info("iAnnotateSV::AddExternalAnnotations: Finished Reading %s", args.rrFilename)
        logging.info("iAnnotateSV::AddExternalAnnotations: Reading %s...", args.dgvFilename)
        dgvDict = afd.ReadDGvFile(args.dgvFilename, args.verbose)
        logging.info("Finished Reading %s", args.dgvFilename)
        data['Cosmic_Fusion_Counts'] = "-"
        data['repName-repClass-repFamily:-site1'] = "-"
        data['repName-repClass-repFamily:-site2'] = "-"
        data['CC_Chr_Band'] = "-"
        data['CC_Tumour_Types(Somatic)'] = "-"
        data['CC_Cancer_Syndrome'] = "-"
        data['CC_Mutation_Type'] = "-"
        data['CC_Translocation_Partner'] = "-"
        data['DGv_Name-DGv_VarType-site1'] = "-"
        data['DGv_Name-DGv_VarType-site2'] = "-"
        for count, row in data.iterrows():
            sv_chr1 = row.loc['chr1']
            sv_pos1 = row.loc['pos1']
            sv_chr2 = row.loc['chr2']
            sv_pos2 = row.loc['pos2']
            sv_gene1 = row.loc['gene1']
            sv_gene2 = row.loc['gene2']
            logging.info("iAnnotateSV::AddExternalAnnotations: Processing Record: %s\t%s\t%s\t%s\t%s\t%s", sv_chr1, sv_pos1, sv_chr2, sv_pos2, sv_gene1, sv_gene2)
            # Repeat Region Data
            (rr_loc1, rr_loc2) = afr.AnnotateRepeatRegion(
                args.verbose, count, row.copy(), repeatRegionDict)
            data.ix[count, 'repName-repClass-repFamily:-site1'] = "<=>".join(rr_loc1)
            data.ix[count, 'repName-repClass-repFamily:-site2'] = "<=>".join(rr_loc2)
            # Cosmic Consensus Data
            cc_SV = afc.AnnotateFromCosmicCensusFile(args.ccFilename, args.verbose, count, row.copy())
            cct_SV = afc.AnnotateFromCosmicFusionCountsFile(args.cctFilename, args.verbose, count, row.copy())
            ccA, ccB, ccC, ccD, ccE = ([] for i in range(5))
            for cc in cc_SV:
                ccData = cc.split('\t')
                ccA.append(ccData[0])
                ccB.append(ccData[1])
                ccC.append(ccData[2])
                ccD.append(ccData[3])
                ccE.append(ccData[4])
            data.ix[count, 'Cosmic_Fusion_Counts'] = cct_SV
            data.ix[count, 'CC_Chr_Band'] = "<=>".join(ccA)
            data.ix[count, 'CC_Tumour_Types(Somatic)'] = "<=>".join(ccB)
            data.ix[count, 'CC_Cancer_Syndrome'] = "<=>".join(ccC)
            data.ix[count, 'CC_Mutation_Type'] = "<=>".join(ccD)
            data.ix[count, 'CC_Translocation_Partner'] = "<=>".join(ccE)
            # DGvData
            (dgv_loc1, dgv_loc2) = afd.AnnotateDGv(args.verbose, count, row.copy(), dgvDict)
            data.ix[count, 'DGv_Name-DGv_VarType-site1'] = "<=>".join(dgv_loc1)
            data.ix[count, 'DGv_Name-DGv_VarType-site2'] = "<=>".join(dgv_loc2)
    else:
        data = ReadSVFile(args)
        repeatRegionDict = afr.ReadRepeatFile(args.rrFilename, args.verbose)
        dgvDict = afr.ReadRepeatFile(args.dgvFilename, args.verbose)
        for count, row in data.iterrows():
            (rr_loc1, rr_loc2) = afr.AnnotateRepeatRegion(
                args.verbose, count, row.copy(), repeatRegionDict)
            cc_SV = afc.AnnotateFromCosmicCensusFile(args.ccFilename, args.verbose, count, row.copy())
            cct_SV = afc.AnnotateFromCosmicFusionCountsFile(args.cctFilename, args.verbose, count, row.copy())
            (dgv_loc1, dgv_loc2) = afd.AnnotateDGv(args.verbose, count, row.copy(), dgvDict)
            data.ix[count, 'repName-repClass-repFamily:-site1'] = "<=>".join(rr_loc1)
            data.ix[count, 'repName-repClass-repFamily:-site2'] = "<=>".join(rr_loc2)
            ccA, ccB, ccC, ccD = ([] for i in range(4))
            for cc in cc_SV:
                ccData = cc.split('\t')
                ccA.append(ccData[0])
                ccB.append(ccData[1])
                ccC.append(ccData[2])
                ccD.append(ccData[3])
            data.ix[count, 'Cosmic_Fusion_Counts'] = cct_SV
            data.ix[count, 'CC_Chr_Band'] = "<=>".join(ccA)
            data.ix[count, 'CC_Tumour_Types(Somatic)'] = "<=>".join(ccB)
            data.ix[count, 'CC_Cancer_Syndrome'] = "<=>".join(ccC)
            data.ix[count, 'CC_Mutation_Type'] = "<=>".join(ccD)
            data.ix[count, 'CC_Translocation_Partner'] = "<=>".join(ccE)
            data.ix[count, 'DGv_Name-DGv_VarType-site1'] = "<=>".join(dgv_loc1)
            data.ix[count, 'DGv_Name-DGv_VarType-site2'] = "<=>".join(dgv_loc2)

    # Print to TSV file
    data.to_csv(outFileTxt, sep='\t', index=False)
    # Print to Json
    data.to_json(outFileJson)
    # Print to Excel
    data.to_excel(outFileExl, sheet_name='Annotated_SVs', index=False)


def ReadSVFile(args):
    if(args.verbose):
        logging.info("Reading Structural Variant File")
        count = len(open(args.svFilename).readlines())
        if(count > 1):
            data = pd.read_csv(args.svFilename, sep='\t', header=0, keep_default_na='True')
            return (data)
        else:
            if(args.verbose):
                logging.warn("File %s doesnot have any structural variants to annotate.", (args.svFilename))
            data = pd.read_csv(args.svFilename, sep='\t', header=0, keep_default_na='True')
            outFileTxt = args.outDir + "/" + args.outFilePrefix + ".txt"
            outFileExl = args.outDir + "/" + args.outFilePrefix + ".xlsx"
            outFileJson = args.outDir + "/" + args.outFilePrefix + ".json"
            # Print to TSV file
            data.to_csv(outFileTxt, sep='\t', index=False)
            # Print to Excel
            data.to_excel(outFileExl, sheet_name='Annotated_SVs', index=False)
            # Print to Json
            data.to_json(outFileJson)
            sys.exit()
   


if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    print("Elapsed time was %g seconds" % (end_time - start_time))
