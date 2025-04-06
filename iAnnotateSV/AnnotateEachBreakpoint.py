'''
Created on Dec 29, 2014

@author: Ronak H Shah
'''
import polars as pl
import numpy as np
import helper as hp
import re
import FindTranscript as ft
from models import *
import logging
from rich import print
from rich.logging import RichHandler

FORMAT = "%(message)s"
logging.basicConfig(
    level="INFO", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)

log = logging.getLogger("rich")


def AnnotateEachBreakpoint(chromosome, position, strand, df, autoSelect):
    """
    Annotates a genomic breakpoint with transcript information.

    Args:
        chromosome (str): Chromosome of the breakpoint.
        position (int): Genomic position of the breakpoint.
        strand (str): Strand of the breakpoint.
        df (pl.DataFrame): DataFrame containing transcript annotations.
        autoSelect (bool): Whether to automatically select the best transcript.

    Returns:
        tuple: A tuple containing gene name, transcript ID, description, zone,
               strand direction, intron number, and intron frame.

    Raises:
        ChrError: If the chromosome is invalid.
        IntergenicError: If the breakpoint is in an intergenic region and no
                         nearby genes are found.
    """
    log.debug(f"Annotating breakpoint: {chromosome}:{position}")

    chromosome = chromosome if chromosome.startswith('chr') else f"chr{chromosome}"

    if not re.match(r"(chr[1-9]$|chr1[0-9]$|chr2[0-2]$|chr[X|Y]$)", chromosome):
        log.error(f"Invalid chromosome: {chromosome}")
        raise ChrError(":".join([str(chromosome), str(position)]))

    # Find all transcripts on the chromosome
    idxList = df.filter(pl.col('chrom') == chromosome).select(pl.col('chrom')).to_series().to_list()
    log.debug(f"Found {len(idxList)} transcripts on chromosome {chromosome}")

    transcriptIndex = []
    # Find all overlapping transcripts
    for index in (idxList):
        geneStart = df[index, 'geneStart']
        geneEnd = df[index, 'geneEnd']
        if geneStart <= position and geneEnd >= position:
            transcriptIndex.append(index)
            log.debug(f"Transcript {df[index, '#name']} overlaps breakpoint")

    geneName, transcript, desc, zone, strandDirection, intronnum, intronframe = (
        None,) * 7

    if transcriptIndex:
        log.debug(f"Found {len(transcriptIndex)} overlapping transcripts")

        coordData = pl.DataFrame({
            'c': [None] * len(transcriptIndex),
            'd': [None] * len(transcriptIndex),
            'e': [None] * len(transcriptIndex),
            'd1': [None] * len(transcriptIndex),
            'd2': [None] * len(transcriptIndex),
            'e1': [None] * len(transcriptIndex),
            'e2': [None] * len(transcriptIndex),
            'f': [None] * len(transcriptIndex)
        })

        for tindex in transcriptIndex:
            log.debug(f"Processing transcript {df[tindex, '#name']}")

            # Annotate coords with each transcript
            c, d, e, d1, d2, e1, e2, f = (None,) * 8  # Reset variables

            # in promoter region ?
            if position < df[tindex, 'txStart']:
                c = 5
                d = df[tindex, 'txStart'] - position
                coordData = coordData.with_columns([
                    pl.Series(name='c', values=[c]),
                    pl.Series(name='d', values=[d]),
                    pl.Series(name='e', values=[e]),
                    pl.Series(name='d1', values=[d1]),
                    pl.Series(name='d2', values=[d2]),
                    pl.Series(name='e1', values=[e1]),
                    pl.Series(name='e2', values=[e2]),
                    pl.Series(name='f', values=[f])
                ])
                log.debug(f"Breakpoint in promoter region (5'), distance={d}")
                continue

            if position > df[tindex, 'txEnd']:
                c = 5
                d = position - df[tindex, 'txEnd']
                coordData = coordData.with_columns([
                    pl.Series(name='c', values=[c]),
                    pl.Series(name='d', values=[d]),
                    pl.Series(name='e', values=[e]),
                    pl.Series(name='d1', values=[d1]),
                    pl.Series(name='d2', values=[d2]),
                    pl.Series(name='e1', values=[e1]),
                    pl.Series(name='e2', values=[e2]),
                    pl.Series(name='f', values=[f])
                ])
                log.debug(f"Breakpoint in promoter region (3'), distance={d}")
                continue

            # in UTR region ?
            if df[tindex, 'strand'] == '+':
                if df[tindex, 'cdsStart'] > position:
                    c = 4
                    d = df[tindex, 'cdsStart'] - position
                    coordData = coordData.with_columns([
                        pl.Series(name='c', values=[c]),
                        pl.Series(name='d', values=[d]),
                        pl.Series(name='e', values=[e]),
                        pl.Series(name='d1', values=[d1]),
                        pl.Series(name='d2', values=[d2]),
                        pl.Series(name='e1', values=[e1]),
                        pl.Series(name='e2', values=[e2]),
                        pl.Series(name='f', values=[f])
                    ])
                    log.debug(f"Breakpoint in 5' UTR, distance={d}")
                    continue
                if position > df[tindex, 'cdsEnd']:
                    c = 3
                    d = position - df[tindex, 'cdsStart']
                    coordData = coordData.with_columns([
                        pl.Series(name='c', values=[c]),
                        pl.Series(name='d', values=[d]),
                        pl.Series(name='e', values=[e]),
                        pl.Series(name='d1', values=[d1]),
                        pl.Series(name='d2', values=[d2]),
                        pl.Series(name='e1', values=[e1]),
                        pl.Series(name='e2', values=[e2]),
                        pl.Series(name='f', values=[f])
                    ])
                    log.debug(f"Breakpoint in 3' UTR, distance={d}")
                    continue
            else:  # negative strand
                if df[tindex, 'cdsStart'] > position:
                    c = 3
                    d = (df[tindex, 'cdsStart'] - position)
                    coordData = coordData.with_columns([
                        pl.Series(name='c', values=[c]),
                        pl.Series(name='d', values=[d]),
                        pl.Series(name='e', values=[e]),
                        pl.Series(name='d1', values=[d1]),
                        pl.Series(name='d2', values=[d2]),
                        pl.Series(name='e1', values=[e1]),
                        pl.Series(name='e2', values=[e2]),
                        pl.Series(name='f', values=[f])
                    ])
                    log.debug(f"Breakpoint in 3' UTR (negative strand), distance={d}")
                    continue
                if position > df[tindex, 'cdsEnd']:
                    c = 4
                    d = position - df[tindex, 'cdsStart']
                    coordData = coordData.with_columns([
                        pl.Series(name='c', values=[c]),
                        pl.Series(name='d', values=[d]),
                        pl.Series(name='e', values=[e]),
                        pl.Series(name='d1', values=[d1]),
                        pl.Series(name='d2', values=[d2]),
                        pl.Series(name='e1', values=[e1]),
                        pl.Series(name='e2', values=[e2]),
                        pl.Series(name='f', values=[f])
                    ])
                    log.debug(f"Breakpoint in 5' UTR (negative strand), distance={d}")
                    continue

            # In exonic region
            exonStarts = list(
                filter(None, df[tindex, 'exonStarts'].split(",")))
            exonEnds = list(filter(None, df[tindex, 'exonEnds'].split(",")))
            in_exon = None
            for k in range(len(exonStarts)):
                if int(exonStarts[k]) <= int(position) and int(exonEnds[k]) >= int(position):
                    in_exon = k + 1
                    break
            if in_exon:
                c = 1
                e = in_exon
                coordData = coordData.with_columns([
                    pl.Series(name='c', values=[c]),
                    pl.Series(name='d', values=[d]),
                    pl.Series(name='e', values=[e]),
                    pl.Series(name='d1', values=[d1]),
                    pl.Series(name='d2', values=[d2]),
                    pl.Series(name='e1', values=[e1]),
                    pl.Series(name='e2', values=[e2]),
                    pl.Series(name='f', values=[f])
                ])
                log.debug(f"Breakpoint in exon {in_exon}")
                continue

            # In Intronic Region
            c = 2
            exonCount = df[tindex, 'exonCount']
            exonFrames = list(
                filter(None, df[tindex, 'exonFrames'].split(",")))
            for k in range(exonCount - 1):  # Up to exonCount - 1 to avoid index out of bounds
                if int(exonEnds[k]) < int(position) and int(exonStarts[k + 1]) > position:
                    f = exonFrames[k + 1] if df[tindex, 'strand'] == '+' else exonFrames[k]
                    e1 = k + 1
                    e2 = k + 2
                    d1 = int(position) - int(exonEnds[k])
                    d2 = int(exonStarts[k + 1]) - int(position)
                    d = min(d1, d2)
                    coordData = coordData.with_columns([
                        pl.Series(name='c', values=[c]),
                        pl.Series(name='d', values=[d]),
                        pl.Series(name='e', values=[e]),
                        pl.Series(name='d1', values=[d1]),
                        pl.Series(name='d2', values=[d2]),
                        pl.Series(name='e1', values=[e1]),
                        pl.Series(name='e2', values=[e2]),
                        pl.Series(name='f', values=[f])
                    ])
                    log.debug(f"Breakpoint in intron between exons {e1} and {e2}, distances d1={d1}, d2={d2}")
                    break

        if autoSelect:
            geneName, transcript, desc, zone, strandDirection, intronnum, intronframe = ft.FindATranscript(
                coordData, df)
            log.info(f"Auto-selected transcript: {transcript}, gene: {geneName}")
        else:
            geneNameList, transcriptList, descList, zoneList, strandDirectionList, intronnumList, intronframeList = ft.FindAllTranscripts(
                coordData, df)
            log.info(f"Found multiple transcripts: {transcriptList}, genes: {geneNameList}")
            geneName, transcript, desc, zone, strandDirection, intronnum, intronframe = geneNameList, transcriptList, descList, zoneList, strandDirectionList, intronnumList, intronframeList
    else:
        log.warn("No overlapping transcripts found. Searching for nearest genes.")
        distBefore = abs(df[idxList, 'txStart'] - position)
        distAfter = abs(df[idxList, 'txEnd'] - position)
        geneName, transcript, desc, zone, strandDirection = None, None, None, None, None
        for y in 1000.0 ** (np.arange(1, 4, 0.3)):
            cmpDB = distBefore[distBefore <= y]
            if not cmpDB.is_empty():
                beforeIdx = cmpDB.arg_min()
                strandDirection = df[beforeIdx, 'strand']
                transcript = df[beforeIdx, '#name']
                zone = 0
                geneName = df[beforeIdx, 'name2']
                desc = f'IGR: {hp.bp2str(distBefore[beforeIdx], 2)} before {geneName}({strandDirection})'
                log.info(f"Found gene before breakpoint: {geneName}, distance={distBefore[beforeIdx]}")
                break

            cmpDA = distAfter[distAfter <= y]
            if not cmpDA.is_empty():
                afterIdx = cmpDA.arg_min()
                geneName = df[afterIdx, 'name2']
                strandDirection = df[afterIdx, 'strand']
                transcript = df[afterIdx, '#name']
                zone = 0
                desc = f'IGR: {hp.bp2str(distAfter[afterIdx], 2)} after {geneName}({strandDirection})'
                log.info(f"Found gene after breakpoint: {geneName}, distance={distAfter[afterIdx]}")
                break

        if not all([geneName, transcript, desc]):
            log.error("Breakpoint is in intergenic region with no nearby genes.")
            raise IntergenicError(":".join([str(chromosome), str(position)]))

    return geneName, transcript, desc, zone, strandDirection, intronnum, intronframe