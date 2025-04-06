'''
Created on Apr 15, 2015
Description: This will help to plot SV using Genome Diagram from bio python
@author: Ronak H Shah
'''
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from reportlab.lib.colors import red, grey, orange, green, brown, blue, lightblue, purple
import sys
import os
import logging
from PIL import Image
import polars as pl
from rich import print
from rich.logging import RichHandler

FORMAT = "%(message)s"
logging.basicConfig(
    level="INFO", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)

log = logging.getLogger("rich")


def VisualizeSV(svDF, refDF, upDF, args):
    """
    Visualizes structural variants using GenomeDiagram.

    Args:
        svDF (pl.DataFrame): DataFrame containing structural variant information.
        refDF (pl.DataFrame): DataFrame containing reference annotation information.
        upDF (pl.DataFrame): DataFrame containing UniProt annotation information.
        args (argparse.Namespace): Command-line arguments.
    """
    staticDir = f"{args.outFilePrefix}_iAnnotateSVplots"
    AnalysisDir = os.path.join(args.outDir, staticDir)
    try:
        os.mkdir(AnalysisDir)
        log.info(f"Created directory: {AnalysisDir}")
    except OSError as e:
        if args.verbose:
            log.warn(f"iAnnotateSV::VisualizeSV: Directory {AnalysisDir} exists. Results will be overwritten. Error: {e}")

    for row in svDF.iter_rows(named=True):
        chr1, chr2 = str(row['chr1']), str(row['chr2'])
        chr1 = chr1 if chr1.startswith('chr') else f"chr{chr1}"
        chr2 = chr2 if chr2.startswith('chr') else f"chr{chr2}"
        fusion = str(row['fusion'])

        if fusion != "-":
            transcript1, transcript2 = str(row['transcript1']), str(row['transcript2'])
            domain1Idx, maxLen1, minLen1 = processData(chr1, transcript1, refDF, upDF)
            domain2Idx, maxLen2, minLen2 = processData(chr2, transcript2, refDF, upDF)

            eventtype = None
            str1, str2 = int(row['str1']), int(row['str2'])

            if chr1 == chr2:
                if str1 == 0 and str2 == 0:
                    eventtype = "Inversion"
                elif str1 == 1 and str2 == 1:
                    eventtype = "Inversion"
                elif str1 == 1 and str2 == 0:
                    eventtype = "Duplication"
                elif str1 == 0 and str2 == 1:
                    eventtype = "Deletion"
            else:
                eventtype = "Translocation"

            pos1, pos2 = int(row['pos1']), int(row['pos2'])
            gene1, gene2 = str(row['gene1']), str(row['gene2'])
            site1, site2 = str(row['site1']), str(row['site2'])

            if gene1 != gene2:
                baseFileName = f"{AnalysisDir}/{gene1}-{chr1}_{pos1}_{gene2}-{chr2}_{pos2}_{eventtype}"
                outFile1Name = f"{baseFileName}-part1.jpg"
                outFile2Name = f"{baseFileName}-part2.jpg"
                outFileName = f"{baseFileName}.jpg"

                d1Name, d2Name = f"{eventtype}-{gene1}", f"{eventtype}-{gene2}"

                gdd1, gdd2 = GenomeDiagram.Diagram(d1Name), GenomeDiagram.Diagram(d2Name)

                feature1Name, feature2Name, feature3Name = f"GeneTrack:{gene1}:{eventtype}", f"AlignmentTrack:{gene1}:{eventtype}", f"DomainTrack:{gene1}:{eventtype}"
                feature4Name, feature5Name, feature6Name = f"GeneTrack:{gene2}:{eventtype}", f"AlignmentTrack:{gene2}:{eventtype}", f"DomainTrack:{gene2}:{eventtype}"

                gdt1_features = gdd1.new_track(1, greytrack=True, name=feature1Name)
                gdt2_features = gdd1.new_track(1, greytrack=True, name=feature2Name)
                gdt3_features = gdd1.new_track(1, greytrack=True, name=feature3Name)
                gdt4_features = gdd2.new_track(1, greytrack=True, name=feature4Name)
                gdt5_features = gdd2.new_track(1, greytrack=True, name=feature5Name)
                gdt6_features = gdd2.new_track(1, greytrack=True, name=feature6Name)

                gds_features = gdt1_features.new_set()
                gds_features = makeReferenceFeatures(transcript1, site1, chr1, pos1, refDF, gds_features)
                gds_features = gdt2_features.new_set()
                gds_features = makeReadFeatures(chr1, pos1, str1, gds_features)
                gds_features = gdt3_features.new_set()
                if domain1Idx:
                    gds_features = makeUniProtFeatures(domain1Idx, upDF, gds_features)

                gdd1.draw(format='linear', fragments=1, start=minLen1 - 1000, end=maxLen1 + 1000)
                gdd1.write(outFile1Name, "JPG", dpi=300)

                gds_features = gdt4_features.new_set()
                gds_features = makeReferenceFeatures(transcript2, site2, chr2, pos2, refDF, gds_features)
                gds_features = gdt5_features.new_set()
                gds_features = makeReadFeatures(chr2, pos2, str2, gds_features)
                gds_features = gdt6_features.new_set()
                if domain2Idx:
                    gds_features = makeUniProtFeatures(domain2Idx, upDF, gds_features)

                gdd2.draw(format='linear', fragments=1, start=minLen2 - 1000, end=maxLen2 + 1000)
                gdd2.write(outFile2Name, "JPG", dpi=300)

                try:
                    img1, img2 = Image.open(outFile1Name), Image.open(outFile2Name)
                    images = [img1, img2]
                    w = max(i.size[0] for i in images)
                    mh = sum(i.size[1] for i in images)
                    result = Image.new("RGB", (w, mh), (255, 255, 255))
                    x = 0
                    for i in images:
                        result.paste(i, (0, x))
                        x += i.size[1]
                    result.save(outFileName)
                    log.info(f"Successfully merged images to {outFileName}")

                    if os.path.isfile(outFileName):
                        os.remove(outFile1Name)
                        os.remove(outFile2Name)
                        log.debug(f"Removed temporary files {outFile1Name} and {outFile2Name}")

                except Exception as e:
                    log.error(f"Error processing or merging images: {e}")

            else:
                outFileName = f"{AnalysisDir}/{gene1}-{chr1}_{pos1}_{gene2}-{chr2}_{pos2}_{eventtype}.jpg"
                gdd = GenomeDiagram.Diagram('Test Diagram')

                feature1Name, feature2Name, feature3Name = f"GeneTrack:{gene1}:{eventtype}", f"AlignmentTrack:{gene1}:{eventtype}", f"ProteinDomainTrack:{gene1}:{eventtype}"
                feature4Name, feature5Name, feature6Name = f"GeneTrack:{gene2}:{eventtype}", f"AlignmentTrack:{gene2}:{eventtype}", f"ProteinDomainTrack:{gene2}:{eventtype}"

                gdt1_features = gdd.new_track(1, greytrack=True, name=feature1Name)
                gdt2_features = gdd.new_track(1, greytrack=True, name=feature2Name)
                gdt3_features = gdd.new_track(1, greytrack=True, name=feature3Name)
                gdt4_features = gdd.new_track(1, greytrack=True, name=feature4Name)
                gdt5_features = gdd.new_track(1, greytrack=True, name=feature5Name)
                gdt6_features = gdd.new_track(1, greytrack=True, name=feature6Name)

                gds_features = gdt1_features.new_set()
                gds_features = makeReferenceFeatures(transcript1, site1, chr1, pos1, refDF, gds_features)
                gds_features = gdt2_features.new_set()
                gds_features = makeReadFeatures(chr1, pos1, str1, gds_features)
                gds_features = gdt3_features.new_set()
                if domain1Idx:
                    gds_features = makeUniProtFeatures(domain1Idx, upDF, gds_features)

                gds_features = gdt4_features.new_set()
                gds_features = makeReferenceFeatures(transcript2, site2, chr2, pos2, refDF, gds_features)
                gds_features = gdt5_features.new_set()
                gds_features = makeReadFeatures(chr2, pos2, str2, gds_features)
                gds_features = gdt6_features.new_set()
                if domain2Idx:
                    gds_features = makeUniProtFeatures(domain2Idx, upDF, gds_features)

                max_len = max(maxLen1, maxLen2)
                min_len = min(minLen1, minLen2)

                gdd.draw(format='linear', pagesize='A4', fragments=1, start=min_len - 1000, end=max_len + 1000)
                gdd.write(outFileName, "JPG", dpi=300)
                log.info(f"Successfully created combined image: {outFileName}")


def processData(chrom, transcript, refDF, upDF):
    """
    Processes transcript and UniProt data to find overlapping domain information.

    Args:
        chrom (str): Chromosome.
        transcript (str): Transcript ID.
        refDF (pl.DataFrame): DataFrame containing reference transcript annotations.
        upDF (pl.DataFrame): DataFrame containing UniProt annotations.

    Returns:
        tuple: A tuple containing UniProt record indices, max length, and min length.
    """
    log.debug(f"Processing data for chrom: {chrom}, transcript: {transcript}")
    transcripts = refDF.filter((pl.col('#name') == transcript) & (pl.col('chrom') == chrom))
    if transcripts.height > 0:
        transcriptIdx = 0
    else:
        transcripts = refDF.filter(pl.col('#name') == transcript)
        if transcripts.height > 0:
            transcriptIdx = 0
        else:
            log.warn(f"No transcript found for {transcript} on chromosome {chrom}")
            return ([], 0, 0)

    refTxSt = int(transcripts[transcriptIdx, 'txStart'])
    refTxEn = int(transcripts[transcriptIdx, 'txEnd'])
    log.debug(f"Transcript start: {refTxSt}, end: {refTxEn}")

    up_idxList = upDF.filter(pl.col('#chrom') == chrom).select(pl.col('#chrom')).to_series().to_list()
    up_recordIndex = []
    for index in up_idxList:
        chromStart, chromEnd = upDF[index, 'chromStart'], upDF[index, 'chromEnd']
        if (chromStart >= refTxSt) and (chromEnd <= refTxEn) and upDF[index, 'annotationType'] == 'domain':
            up_recordIndex.append(index)
            log.debug(f"Found overlapping domain at {chromStart}-{chromEnd}")

    allMaxVal = [upDF[index, 'chromEnd'] for index in up_recordIndex]
    allMinVal = [upDF[index, 'chromStart'] for index in up_recordIndex]

    max_len = max(allMaxVal, default=refTxEn)
    min_len = max(allMinVal, default=refTxSt)
    log.debug(f"Max length: {max_len}, Min length: {min_len}")

    return (up_recordIndex, max_len, min_len)


def makeReferenceFeatures(transcript, site, chrom, pos, refDF, gds_features):
    """
    Creates reference features for GenomeDiagram.

    Args:
        transcript (str): Transcript ID.
        site (str): Site information.
        chrom (str): Chromosome.
        pos (int): Position.
        refDF (pl.DataFrame): DataFrame containing reference annotation information.
        gds_features (GenomeDiagram.FeatureSet): Feature set to add features to.

    Returns:
        GenomeDiagram.FeatureSet: Updated feature set.
    """
    log.debug(f"Making reference features for transcript: {transcript}, site: {site}, chrom: {chrom}, pos: {pos}")
    transcripts = refDF.filter((pl.col('#name') == transcript) & (pl.col('chrom') == chrom))
    if transcripts.height > 0:
        transcriptIdx = 0
    else:
        transcripts = refDF.filter(pl.col('#name') == transcript)
        if transcripts.height > 0:
            transcriptIdx = 0
        else:
            log.warn(f"No transcript found for {transcript} on chromosome {chrom}")
            return gds_features

    refTxSt = int(transcripts[transcriptIdx, 'txStart'])
    refTxEn = int(transcripts[transcriptIdx, 'txEnd'])
    log.debug(f"Transcript start: {refTxSt}, end: {refTxEn}")

    ExonSts = list(filter(None, transcripts[transcriptIdx, 'exonStarts'].split(",")))
    ExonEnds = list(filter(None, transcripts[transcriptIdx, 'exonEnds'].split(",")))
    transcriptStrand = str(transcripts[transcriptIdx, 'strand'])
    transcriptStrand = +1 if transcriptStrand == "+" else -1 if transcriptStrand == "-" else None

    for idx, val in enumerate(ExonSts):
        feature = SeqFeature(FeatureLocation(int(val), int(ExonEnds[idx]), strand=transcriptStrand))
        fname = f"exon{len(ExonSts) - idx if transcriptStrand == -1 else idx + 1}"
        gds_features.add_feature(
            feature,
            sigil="ARROW",
            color=brown,
            arrowshaft_height=1.0,
            name=fname,
            label=True, label_position="middle", label_size=5, label_angle=90 if transcriptStrand == -1 else 0)
        log.debug(f"Added exon feature: {fname} at {val}-{ExonEnds[idx]}")

    feature = SeqFeature(FeatureLocation(pos - 5, pos + 5))
    gds_features.add_feature(
        feature,
        color=orange,
        name=site,
        label=True,
        label_size=6, label_color=orange)
    log.debug(f"Added site feature: {site} at {pos - 5}-{pos + 5}")

    return gds_features


def makeUniProtFeatures(domainIdx, upDF, gds_features):
    """
    Creates UniProt domain features for GenomeDiagram.

    Args:
        domainIdx (list): List of UniProt record indices.
        upDF (pl.DataFrame): DataFrame containing UniProt annotation information.
        gds_features (GenomeDiagram.FeatureSet): Feature set to add features to.

    Returns:
        GenomeDiagram.FeatureSet: Updated feature set.
    """
    log.debug(f"Making UniProt features for {len(domainIdx)} domains")
    for index in domainIdx:
        chromStart, chromEnd = upDF[index, 'chromStart'], upDF[index, 'chromEnd']
        fname = upDF[index, 'name']
        feature = SeqFeature(FeatureLocation(chromStart, chromEnd), strand=None)
        gds_features.add_feature(
            feature,
            name=fname,
            label=True,
            color=green,
            label_position="middle",
            label_size=6,
            label_color=green)
        log.debug(f"Added UniProt feature: {fname} at {chromStart}-{chromEnd}")
    return gds_features


def makeReadFeatures(chrom, pos, strand, gds_features):
    """
    Creates read alignment features for GenomeDiagram.

    Args:
        chrom (str): Chromosome.
        pos (int): Position.
        strand (int): Strand (0 or 1).
        gds_features (GenomeDiagram.FeatureSet): Feature set to add features to.

    Returns:
        GenomeDiagram.FeatureSet: Updated feature set.
    """
    log.debug(f"Making read features for chrom: {chrom}, pos: {pos}, strand: {strand}")
    start, end = int(pos) - 1000, int(pos) + 1000
    bname = f"{chrom}:{pos}"
    strandDirection = +1 if strand == 0 else -1
    color = blue if strand == 0 else red

    feature = SeqFeature(FeatureLocation(start, end, strand=strandDirection))
    gds_features.add_feature(
        feature,
        sigil="ARROW",
        arrowshaft_height=0.1,
        color=color,
        name=bname, label_size=8, label=True, label_angle=0 if strand == 0 else -90, label_color=purple)
    log.debug(f"Added read feature: {bname} at {start}-{end}, strand: {strandDirection}")

    return gds_features


def makePlainImage(
        refDF,
        eventtype,
        transcript1,
        transcript2,
        chr1,
        chr2,
        pos1,
        pos2,
        str1,
        str2,
        site1,
        site2,
        fusion,
        gds_features):
    """
    Creates a plain image for visualization (work in progress).

    Args:
        refDF (pl.DataFrame): DataFrame containing reference annotation information.
        eventtype (str): Type of event.
        transcript1 (str): Transcript ID for the first breakpoint.
        transcript2 (str): Transcript ID for the second breakpoint.
        chr1 (str): Chromosome for the first breakpoint.
        chr2 (str): Chromosome for the second breakpoint.
        pos1 (int): Position for the first breakpoint.
        pos2 (int): Position for the second breakpoint.
        str1 (int): Strand for the first breakpoint.
        str2 (int): Strand for the second breakpoint.
        site1 (str): Site information for the first breakpoint.
        site2 (str): Site information for the second breakpoint.
        fusion (str): Fusion information.
        gds_features (GenomeDiagram.FeatureSet): Feature set to add features to.
    """
    log.debug("Making plain image (work in progress)")
    ExonCounts1 = _extracted_from_makePlainImage_16(refDF, transcript1)
    ExonCounts2 = _extracted_from_makePlainImage_16(refDF, transcript2)
    before_exonNum1, after_exonNum1, before_exonNum2, after_exonNum2 = None, None, None, None

    if "before" in site1:
        before_exonNum1 = site1[-1:]
    if "after" in site1:
        after_exonNum1 = site1[-1:]
    if "before" in site2:
        before_exonNum2 = site2[-1:]
    if "after" in site2:
        after_exonNum2 = site2[-1:]

    beforeExons1, afterExons1, beforeExons2, afterExons2 = [], [], [], []
    for i in range(1, ExonCounts1):
        if before_exonNum1 and i <= int(before_exonNum1):
            beforeExons1.append(f"exon{i}")
        if after_exonNum1 and i >= int(after_exonNum1):
            afterExons1.append(f"exon{i}")
    for i in range(1, ExonCounts2):
        if before_exonNum2 and i <= int(before_exonNum2):
            beforeExons2.append(f"exon{i}")
        if after_exonNum2 and i >= int(after_exonNum2):
            afterExons2.append(f"exon{i}")


def _extracted_from_makePlainImage_16(refDF, arg1):
    """
    Extracts exon-related information from the reference DataFrame.

    Args:
        refDF (pl.DataFrame): DataFrame containing reference annotation information.
        arg1 (str): Transcript ID.

    Returns:
        int: Exon count.
    """
    transcript1Idx, = refDF[refDF['#name'] == arg1].index
    ExonSts1 = filter(None, refDF[transcript1Idx, 'exonStarts'].split(","))
    ExonEnds1 = filter(None, refDF[transcript1Idx, 'exonEnds'].split(","))
    transcript1Strand = str(refDF[transcript1Idx, 'strand'])
    return int(refDF[transcript1Idx, 'exonCount'])


def get_concat_v_resize(im1, im2, resample=Image.BICUBIC, resize_big_image=True):
    """
    Concatenates two images vertically with resizing.

    Args:
        im1 (PIL.Image.Image): First image.
        im2 (PIL.Image.Image): Second image.
        resample (int): Resampling filter.
        resize_big_image (bool): Whether to resize the larger image.

    Returns:
        PIL.Image.Image: Concatenated image.
    """
    if im1.width == im2.width:
        _im1, _im2 = im1, im2
    elif ((im1.width > im2.width and resize_big_image) or (im1.width < im2.width and not resize_big_image)):
        _im1 = im1.resize((im2.width, int(im1.height * im2.width / im1.width)), resample=resample)
        _im2 = im2
    else:
        _im1, _im2 = im1, im2.resize((im1.width, int(im2.height * im1.width / im2.width)), resample=resample)

    dst = Image.new('RGB', (_im1.width, _im1.height + _im2.height))
    dst.paste(_im1, (0, 0))
    dst.paste(_im2, (0, _im1.height))
    return dst