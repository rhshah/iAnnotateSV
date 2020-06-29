"""
Created on 25/11/2014.

@author: Ronak H Shah

"""

import argparse
import time
import pandas as pd
import helper as hp
import AnnotateEachBreakpoint as aeb
import PredictFunction as pf
import FindCanonicalTranscript as fct
import AddExternalAnnotations as aea
import AnnotationForKinaseDomain as kda
import VisualizeSV as vsv
import os
import sys
import logging
import coloredlogs

'''
Driver function to drive the whole process
'''


def main(command=None):

    parser = argparse.ArgumentParser(
        prog='iAnnotateSV.py',
        description='Annotate SV based on a specific human reference',
        usage='%(prog)s [options]')
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        help="make lots of noise [default]")
    parser.add_argument(
        "-r",
        "--refFileVersion",
        action="store",
        dest="refVersion",
        required=True,
        metavar='hg19',
        help="Which human reference file to be used, hg18,hg19 or hg38")
    parser.add_argument(
        "-rf",
        "--refFile",
        action="store",
        dest="refFile",
        required=False,
        metavar='hg19.sv.table.txt',
        help="Human reference file location to be used")
    parser.add_argument(
        "-ofp",
        "--outputFilePrefix",
        action="store",
        dest="outFilePrefix",
        required=True,
        metavar='test',
        help="Prefix for the output file")
    parser.add_argument(
        "-o",
        "--outputDir",
        action="store",
        dest="outDir",
        required=True,
        metavar='/somedir',
        help="Full Path to the output dir")
    parser.add_argument(
        "-i",
        "--svFile",
        action="store",
        dest="svFilename",
        required=True,
        metavar='svfile.txt',
        help="Location of the structural variants file to annotate")
    parser.add_argument(
        "-d",
        "--distance",
        action="store",
        dest="distance",
        default=3000,
        required=False,
        metavar='3000',
        help="Distance used to extend the promoter region")
    parser.add_argument(
        "-a",
        "--autoSelect",
        action="store_true",
        dest="autoSelect",
        default=True,
        help="Auto Select which transcript to be used[default]")
    parser.add_argument(
        "-c",
        "--canonicalTranscripts",
        action="store",
        dest="canonicalTranscripts",
        required=False,
        metavar='canonicalExons.txt',
        help="Location of canonical transcript list for each gene. Use only if you want the output for specific transcripts for each gene.")
    parser.add_argument(
        "-p",
        "--plotSV",
        action="store_true",
        dest="plotSV",
        help="Plot the structural variant in question")
    parser.add_argument(
        "-u",
        "--uniprotFile",
        action="store",
        dest="uniprot",
        required=False,
        metavar='uniprot.txt',
        help="Location of UniProt list contain information for protein domains. Use only if you want to plot the structural variant")
    parser.add_argument(
        "-rr",
        "--repeatFile",
        action="store",
        dest="rrFilename",
        required=False,
        metavar='RepeatRegionFile.tsv',
        help="Location of the Repeat Region Bed File")
    parser.add_argument(
        "-dgv",
        "--dgvFile",
        action="store",
        dest="dgvFilename",
        required=False,
        metavar='DGvFile.tsv',
        help="Location of the Database of Genomic Variants Bed File")
    parser.add_argument(
        "-cc",
        "--cosmicConsensusFile",
        action="store",
        dest="ccFilename",
        required=False,
        metavar='CosmicConsensus.tsv',
        help="Location of the Cosmic Consensus TSV file")
    parser.add_argument(
        "-cct",
        "--cosmicCountsFile",
        action="store",
        dest="cctFilename",
        required=False,
        metavar='cosmic_fusion_counts.tsv',
        help="Location of the Cosmic Counts TSV file")
    args = ""
    if(command is None):
        args = parser.parse_args()
    else:
        args = parser.parse_args(command.split())

    # Create Logger if verbose
    loggeroutput = args.outDir + "/" + args.outFilePrefix + "_iAnnotateSV.log"
    logging.basicConfig(
        filename=loggeroutput,
        filemode='w',
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logging.DEBUG)

    coloredlogs.install(level='DEBUG')
    # Get current location
    this_dir, this_filename = os.path.split(__file__)
    
    # Check if file for canonical transcript is given or not
    if(args.canonicalTranscripts):
        args.autoSelect = False
    
    if(args.refVersion == 'hg18' or args.refVersion == 'hg19' or args.refVersion == 'hg38'):
        if(args.refFile):
            pass
        else:
            refFile = args.refVersion + ".sv.table.txt"
            refFile = os.path.join(this_dir, "data/references", refFile)
            args.refFile = refFile
        if(args.rrFilename):
            rrPath = args.rrFilename
        else:
            rrFilename = args.refVersion + "_repeatRegion.tsv"
            rrPath = os.path.join(this_dir, "data/repeat_region", rrFilename)
            args.rrFilename = rrPath
        if(args.dgvFilename):
            dgvPath = args.dgvFilename
        else:
            dgvFilename = args.refVersion + "_DGv_Annotation.tsv"
            dgvPath = os.path.join(
                this_dir, "data/database_of_genomic_variants", dgvFilename)
            args.dgvFilename = dgvPath
        if(args.ccFilename):
            ccPath = args.ccFilename
        else:
            ccFilename = "cancer_gene_census.tsv"
            ccPath = os.path.join(this_dir, "data/cosmic", ccFilename)
            args.ccFilename = ccPath
        if(args.cctFilename):
            cctPath = args.cctFilename
        else:
            cctFilename = "cosmic_fusion_counts.tsv"
            cctPath = os.path.join(this_dir, "data/cosmic", cctFilename)
            args.cctFilename = cctPath
        if(args.uniprot):
            uniprotPath = args.uniprot
        else:
            upFilename = args.refVersion + ".uniprot.spAnnot.table.txt"
            args.uniprot = str(os.path.join(
                this_dir, "data/UcscUniprotdomainInfo", upFilename))
            uniprotPath = args.uniprot
        args.allCanonicalTranscriptsPath = str(os.path.join(
            this_dir, "data/canonicalInfo/canonical_transcripts.txt"))
    else:
        if(args.verbose):
            logging.fatal(
                "iAnnotateSV: Please enter correct reference file version. Values can be: hg18 or hg19 or hg38")
            sys.exit()
    (refDF) = hp.ReadFile(refFile)
    NewRefDF = hp.ExtendPromoterRegion(refDF, args.distance)
    svDF = hp.ReadFile(args.svFilename)
    annDF = processSV(svDF, NewRefDF, args)
    plotDF = annDF.copy()
    # Print to TSV file
    outFilePrefixPath = args.outDir + "/" + args.outFilePrefix + "_functional.txt"
    annDF.to_csv(outFilePrefixPath, sep='\t', index=False)
    # Add External Annotations
    if args.verbose:
        logging.info("iAnnotateSV: Adding External Annotations...")
    makeCommandLineForAEA = "-r " + rrPath + " -d " + dgvPath + " -c " + ccPath + " -cct " + cctPath + " -s " + \
        outFilePrefixPath + " -ofp " + args.outFilePrefix + \
        "_Annotated" + " -o " + args.outDir
    aea.main(makeCommandLineForAEA)
    # Plot if required
    if(args.plotSV):
        if args.verbose:
            logging.info("iAnnotateSV: Plotting Each Structural Variants")
        plotSV(plotDF, NewRefDF, uniprotPath, args)

    if(args.verbose):
        logging.info("iAnnotateSV: Finished Running the Annotation Process!!!")


'''
Process Each Structural Variant
'''


def processSV(svDF, refDF, args):
    if args.verbose:
        logging.info("iAnnotateSV: Processing Each Structural Variants...")
    # Read Canonical Transcript if the file is given in the cmdline
    if(args.canonicalTranscripts):
        ctDict = hp.ReadTranscriptFile(args.canonicalTranscripts)
    annDF = pd.DataFrame(
        columns=[
            'chr1',
            'pos1',
            'str1',
            'chr2',
            'pos2',
            'str2',
            'gene1',
            'transcript1',
            'site1',
            'gene2',
            'transcript2',
            'site2',
            'fusion'])
    for count, row in svDF.iterrows():
        # print row
        chr1 = str(row.loc['chr1'])
        chr2 = str(row.loc['chr2'])
        pos1 = int(row.loc['pos1'])
        pos2 = int(row.loc['pos2'])
        str1 = int(row.loc['str1'])
        str2 = int(row.loc['str2'])
        if(args.autoSelect):
            (gene1, transcript1, site1, zone1, strand1, intronnum1,
             intronframe1) = aeb.AnnotateEachBreakpoint(chr1, pos1, str1, refDF, args.autoSelect)
            (gene2, transcript2, site2, zone2, strand2, intronnum2,
             intronframe2) = aeb.AnnotateEachBreakpoint(chr2, pos2, str2, refDF, args.autoSelect)
            ann1S = pd.Series([gene1, transcript1, site1, zone1, strand1, str1, intronnum1, intronframe1], index=[
                              'gene1', 'transcript1', 'site1', 'zone1', 'txstrand1', 'readstrand1', 'intronnum1', 'intronframe1'])
            ann2S = pd.Series([gene2, transcript2, site2, zone2, strand2, str2, intronnum2, intronframe2], index=[
                              'gene2', 'transcript2', 'site2', 'zone2', 'txstrand2', 'readstrand2', 'intronnum2', 'intronframe2'])
            fusionFunction = pf.PredictFunctionForSV(ann1S, ann2S)
            #annS = pd.Series([chr1,pos1,str1,chr2,pos2,str2,gene1,transcript1,site1,gene2,transcript2,site2,fusionFunction],index=['chr1','pos1','str1','chr2','pos2','str2','gene1','transcript1','site1','gene2','transcript2','site2','fusion'])
            annDF.loc[
                count,
                ['chr1', 'pos1', 'str1', 'chr2', 'pos2', 'str2', 'gene1', 'transcript1', 'site1',
                 'gene2', 'transcript2', 'site2', 'fusion']] = [
                chr1, pos1, str1, chr2, pos2, str2, gene1, transcript1, site1, gene2, transcript2,
                site2, fusionFunction]
        else:
            (gene1List, transcript1List, site1List, zone1List, strand1List, intronnum1List,
             intronframe1List) = aeb.AnnotateEachBreakpoint(chr1, pos1, str1, refDF, args.autoSelect)
            # print "1:\n", gene1List, transcript1List, site1List, zone1List, strand1List, intronnum1List, intronframe1List
            (gene2List, transcript2List, site2List, zone2List, strand2List, intronnum2List,
             intronframe2List) = aeb.AnnotateEachBreakpoint(chr2, pos2, str2, refDF, args.autoSelect)
            # print "\n2:\n", gene2List, transcript2List, site2List, zone2List, strand2List, intronnum2List, intronframe2List
            (gene1, transcript1, site1, zone1, strand1, intronnum1, intronframe1) = fct.FindCT(
                gene1List, transcript1List, site1List, zone1List, strand1List, intronnum1List, intronframe1List, ctDict)
            # print "\n3:\n", gene1, transcript1, site1, zone1, strand1, intronnum1, intronframe1
            (gene2, transcript2, site2, zone2, strand2, intronnum2, intronframe2) = fct.FindCT(
                gene2List, transcript2List, site2List, zone2List, strand2List, intronnum2List, intronframe2List, ctDict)
            # print "\n4:\n", gene2, transcript2, site2, zone2, strand2, intronnum2, intronframe2
            ann1S = pd.Series([gene1, transcript1, site1, zone1, strand1, str1, intronnum1, intronframe1], index=[
                              'gene1', 'transcript1', 'site1', 'zone1', 'txstrand1', 'readstrand1', 'intronnum1', 'intronframe1'])
            ann2S = pd.Series([gene2, transcript2, site2, zone2, strand2, str2, intronnum2, intronframe2], index=[
                              'gene2', 'transcript2', 'site2', 'zone2', 'txstrand2', 'readstrand2', 'intronnum2', 'intronframe2'])
            fusionFunction = pf.PredictFunctionForSV(ann1S, ann2S)
            # print "\n5:\n", fusionFunction
            annDF.loc[
                count,
                ['chr1', 'pos1', 'str1', 'chr2', 'pos2', 'str2', 'gene1', 'transcript1', 'site1',
                 'gene2', 'transcript2', 'site2', 'fusion']] = [
                chr1, pos1, str1, chr2, pos2, str2, gene1, transcript1, site1, gene2, transcript2,
                site2, fusionFunction]
    if(args.canonicalTranscripts):
        (svDF) = kda.run(annDF, args.refFile, args.canonicalTranscripts,
                        args.allCanonicalTranscriptsPath, args.uniprot, args.verbose)
        return(svDF)
    else:
        return(annDF)


'''
Plot Annotated Structural Variants
'''


def plotSV(svDF, refDF, uniprotPath, args):
    if args.verbose:
        logging.info(
            "iAnnotateSV: Will now try to plot Each Structural Variants")
    upDF = None
    if(os.path.isfile(uniprotPath)):
        upDF = hp.ReadFile(uniprotPath)
    else:
        if args.verbose:
            logging.fatal(
                "iAnnotateSV: %s file does not exist!!, Please use it to plot structural variants",
                uniprotPath)
            sys.exit()

    vsv.VisualizeSV(svDF, refDF, upDF, args)


'''
Initializing the Driver function
'''

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    logging.info("iAnnotateSV: Elapsed time was %g seconds",
                 (end_time - start_time))
