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
import VisualizeSV as vsv
from models import *
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
        required=True,
        metavar='canonicalExons.txt',
        help="Location of canonical transcript list for each gene.")
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
            refFile = args.refFile
            pass
        else:
            refFile = args.refVersion + ".sv.table.txt"
            refFile = os.path.join(this_dir, "data/references", refFile)
            args.refFile = refFile
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
        ctDict = pd.read_csv(args.canonicalTranscripts, sep='\t', header=0)
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
        id1 = str(row.loc['id1'])
        id2 = str(row.loc['id2'])
        b1, b2 = (None,)*2
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
            try:
                (gene1List, transcript1List, site1List, zone1List, strand1List, intronnum1List, intronframe1List) = aeb.AnnotateEachBreakpoint(chr1, pos1, str1, refDF, args.autoSelect)
                (gene1, transcript1, site1, zone1, strand1, intronnum1, intronframe1) = fct.FindCT(gene1List, transcript1List, site1List, zone1List, strand1List, intronnum1List, intronframe1List, ctDict, id1)
            except (IntergenicError, ChrError) as b1:
                logging.info("iAnnotateSV: " + str(b1))
                (gene1, transcript1, site1, zone1, strand1, intronnum1, intronframe1) = ("-",)*7
            try:
                (gene2List, transcript2List, site2List, zone2List, strand2List, intronnum2List, intronframe2List) = aeb.AnnotateEachBreakpoint(chr2, pos2, str2, refDF, args.autoSelect)
                (gene2, transcript2, site2, zone2, strand2, intronnum2, intronframe2) = fct.FindCT(gene2List, transcript2List, site2List, zone2List, strand2List, intronnum2List, intronframe2List, ctDict, id2)
            except (IntergenicError, ChrError) as b2:
                logging.info("iAnnotateSV: " + str(b2))
                (gene2, transcript2, site2, zone2, strand2, intronnum2, intronframe2) = ("-",)*7
            ann1S = pd.Series([gene1, transcript1, site1, zone1, strand1, str1, intronnum1, intronframe1], index=['gene1', 'transcript1', 'site1', 'zone1', 'txstrand1', 'readstrand1', 'intronnum1', 'intronframe1'])
            ann2S = pd.Series([gene2, transcript2, site2, zone2, strand2, str2, intronnum2, intronframe2], index=['gene2', 'transcript2', 'site2', 'zone2', 'txstrand2', 'readstrand2', 'intronnum2', 'intronframe2'])
            if not any([b1, b2]):
                fusionFunction = pf.PredictFunctionForSV(ann1S, ann2S)
            else:
                fusionFunction = "-"
            annDF.loc[
                count,
                ['chr1', 'pos1', 'str1', 'chr2', 'pos2', 'str2', 'gene1', 'transcript1', 'site1',
                 'gene2', 'transcript2', 'site2', 'fusion']] = [
                chr1, pos1, str1, chr2, pos2, str2, gene1, transcript1, site1, gene2, transcript2,
                site2, fusionFunction]
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
