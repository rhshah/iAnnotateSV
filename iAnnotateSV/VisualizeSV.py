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
import coloredlogs

coloredlogs.install(level='DEBUG')

def VisualizeSV(svDF, refDF, upDF, args):
    staticDir = args.outFilePrefix + "_iAnnotateSVplots"
    AnalysisDir = os.path.join(args.outDir, staticDir)
    try:
        os.mkdir(AnalysisDir)
    except OSError:
        if(args.verbose):
            logging.warn("iAnnotateSV::VisualizeSV: Dir: %s exists thus we wont be making it. Thus Results would be over-written", AnalysisDir)
    for count, row in svDF.iterrows():
        # print row
        chr1 = str(row.loc['chr1'])
        chr2 = str(row.loc['chr2'])
        if(chr1.startswith('chr')):
            chr1 = chr1
        else:
            chr1 = "chr" + chr1
        if(chr2.startswith('chr')):
            chr2 = chr2
        else:
            chr2 = "chr" + chr2
        pos1 = int(row.loc['pos1'])
        pos2 = int(row.loc['pos2'])
        str1 = int(row.loc['str1'])
        str2 = int(row.loc['str2'])
        gene1 = str(row.loc['gene1'])
        gene2 = str(row.loc['gene2'])
        site1 = str(row.loc['site1'])
        site2 = str(row.loc['site2'])
        transcript1 = str(row.loc['transcript1'])
        transcript2 = str(row.loc['transcript2'])
        fusion = str(row.loc['fusion'])
        if(fusion != "-"):
            (domain1Idx, maxLen1, minLen1) = processData(chr1, transcript1, refDF, upDF)
            (domain2Idx, maxLen2, minLen2) = processData(chr2, transcript2, refDF, upDF)
            # print fusion, "\n", "1", domain1Idx, "\n", "2", domain2Idx, "\n\n",
            eventtype = None
            if(str1 == 0 and str2 == 0 and chr1 == chr2):
                eventtype = "Inversion"
            elif(str1 == 1 and str2 == 1 and chr1 == chr2):
                eventtype = "Inversion"
            elif(str1 == 1 and str2 == 0 and chr1 == chr2):
                eventtype = "Duplication"
            elif(str1 == 0 and str2 == 1 and chr1 == chr2):
                eventtype = "Deletion"
            else:
                eventtype = "Translocation"
            #imageMSG = eventtype + " causing " + fusion
            if(gene1 != gene2):
                outFile1Name = AnalysisDir + "/" + gene1 + "-" + str(chr1) + "_" + str(
                    pos1) + "_" + gene2 + "-" + str(chr2) + "_" + str(pos2) + "_" + str(eventtype) + "-part1.jpg"
                outFile2Name = AnalysisDir + "/" + gene1 + "-" + str(chr1) + "_" + str(
                    pos1) + "_" + gene2 + "-" + str(chr2) + "_" + str(pos2) + "_" + str(eventtype) + "-part2.jpg"
                outFileName = AnalysisDir + "/" + gene1 + "-" + str(chr1) + "_" + str(
                    pos1) + "_" + gene2 + "-" + str(chr2) + "_" + str(pos2) + "_" + str(eventtype) + ".jpg"
                d1Name = eventtype + "-" + gene1
                d2Name = eventtype + "-" + gene2
                # Make an instace of class diagram
                gdd1 = GenomeDiagram.Diagram(d1Name)
                gdd2 = GenomeDiagram.Diagram(d2Name)
                # Make name of the tracks
                feature1Name = "GeneTrack:" + gene1 + ":" + eventtype
                feature2Name = "AlignmentTrack:" + gene1 + ":" + eventtype
                feature3Name = "DomainTrack:" + gene1 + ":" + eventtype
                feature4Name = "GeneTrack:" + gene2 + ":" + eventtype
                feature5Name = "AlignmentTrack:" + gene2 + ":" + eventtype
                feature6Name = "DomainTrack:" + gene2 + ":" + eventtype
                # Make track for each feature
                gdt1_features = gdd1.new_track(1, greytrack=True, name=feature1Name)
                gdt2_features = gdd1.new_track(1, greytrack=True, name=feature2Name)
                gdt3_features = gdd1.new_track(1, greytrack=True, name=feature3Name)
                gdt4_features = gdd2.new_track(1, greytrack=True, name=feature4Name)
                gdt5_features = gdd2.new_track(1, greytrack=True, name=feature5Name)
                gdt6_features = gdd2.new_track(1, greytrack=True, name=feature6Name)
                # Write features to a track
                gds_features = gdt1_features.new_set()
                (gds_features) = makeReferenceFeatures(
                    transcript1,
                    site1,
                    chr1,
                    pos1,
                    refDF,
                    gds_features)
                gds_features = gdt2_features.new_set()
                (gds_features) = makeReadFeatures(chr1, pos1, str1, gds_features)
                gds_features = gdt3_features.new_set()
                if(domain1Idx):
                    (gds_features) = makeUniProtFeatures(domain1Idx, upDF, gds_features)
                gdd1.draw(
                    format='linear',
                    #pagesize='A4',
                    fragments=1,
                    start=minLen1 -
                    1000,
                    end=maxLen1 +
                    1000)
                gdd1.write(outFile1Name, "JPG", dpi=300)
                gds_features = gdt4_features.new_set()
                (gds_features) = makeReferenceFeatures(
                    transcript2,
                    site2,
                    chr2,
                    pos2,
                    refDF,
                    gds_features)
                gds_features = gdt5_features.new_set()
                (gds_features) = makeReadFeatures(chr2, pos2, str2, gds_features)
                gds_features = gdt6_features.new_set()
                if(domain2Idx):
                    (gds_features) = makeUniProtFeatures(domain2Idx, upDF, gds_features)
                # draw the object and store in memory
                gdd2.draw(
                    format='linear',
                    #pagesize='A4',
                    fragments=1,
                    start=minLen2 -
                    1000,
                    end=maxLen2 +
                    1000)
                # Write the object to a file
                gdd2.write(outFile2Name, "JPG", dpi=300)
                # merge Images
                img1 = Image.open(outFile1Name)
                img2 = Image.open(outFile2Name)
                images = map(Image.open, [outFile1Name, outFile2Name])
                w = max(i.size[0] for i in images)
                mh = sum(i.size[1] for i in images)
                result = Image.new("RGB", (w, mh), (255, 255, 255))
                x = 0
                for i in images:
                    result.paste(i, (0, x))
                    x += i.size[1]
                result.save(outFileName)
                if(os.path.isfile(outFileName)):
                    os.remove(outFile1Name)
                    os.remove(outFile2Name)

            else:
                outFileName = AnalysisDir + "/" + gene1 + "-" + str(chr1) + "_" + str(
                    pos1) + "_" + gene2 + "-" + str(chr2) + "_" + str(pos2) + "_" + str(eventtype) + ".jpg"
                gdd = GenomeDiagram.Diagram('Test Diagram')
                feature1Name = "GeneTrack:" + gene1 + ":" + eventtype
                feature2Name = "AlignmentTrack:" + gene1 + ":" + eventtype
                feature3Name = "ProteinDomainTrack:" + gene1 + ":" + eventtype
                feature4Name = "GeneTrack:" + gene2 + ":" + eventtype
                feature5Name = "AlignmentTrack:" + gene2 + ":" + eventtype
                feature6Name = "ProteinDomainTrack:" + gene2 + ":" + eventtype
                gdt1_features = gdd.new_track(1, greytrack=True, name=feature1Name)
                gdt2_features = gdd.new_track(1, greytrack=True, name=feature2Name)
                gdt3_features = gdd.new_track(1, greytrack=True, name=feature3Name)
                gdt4_features = gdd.new_track(1, greytrack=True, name=feature4Name)
                gdt5_features = gdd.new_track(1, greytrack=True, name=feature5Name)
                gdt6_features = gdd.new_track(1, greytrack=True, name=feature6Name)
                gds_features = gdt1_features.new_set()
                (gds_features) = makeReferenceFeatures(
                    transcript1,
                    site1,
                    chr1,
                    pos1,
                    refDF,
                    gds_features)
                gds_features = gdt2_features.new_set()
                (gds_features) = makeReadFeatures(chr1, pos1, str1, gds_features)
                gds_features = gdt3_features.new_set()
                if(domain1Idx):
                    (gds_features) = makeUniProtFeatures(domain1Idx, upDF, gds_features)
                gds_features = gdt4_features.new_set()
                (gds_features) = makeReferenceFeatures(
                    transcript2,
                    site2,
                    chr2,
                    pos2,
                    refDF,
                    gds_features)
                gds_features = gdt5_features.new_set()
                (gds_features) = makeReadFeatures(chr2, pos2, str2, gds_features)
                gds_features = gdt6_features.new_set()
                if(domain2Idx):
                    (gds_features) = makeUniProtFeatures(domain2Idx, upDF, gds_features)
                max_len = max(maxLen1, maxLen2)
                min_len = min(minLen1, minLen2)
                gdd.draw(
                    format='linear',
                    pagesize='A4',
                    fragments=1,
                    start=min_len -
                    1000,
                    end=max_len +
                    1000)
                gdd.write(outFileName, "JPG", dpi=300)
            #gdd = GenomeDiagram.Diagram("Fusion Image")
            #featureName = gene1 + ":" + gene2 + ":" + eventtype
            #gdt_features = gdd.new_track(1, greytrack=True, name=featureName)
            #gds_features = gdt_features.new_set()
            # makePlainImage(refDF,eventtype,transcript1,transcript2,chr1,chr2,pos1,pos2,str1,str2,site1,site2,fusion,gds_features)


def processData(chrom, transcript, refDF, upDF):
    transcripts = (refDF[refDF['#name'] == transcript])
    if(len(transcripts) > 1):
        transcriptIdx, = (transcripts[transcripts['chrom'] == chrom].index)
    else:
        transcriptIdx, = (refDF[refDF['#name'] == transcript].index)
    refTxSt = int(refDF.iloc[transcriptIdx]['txStart'])
    refTxEn = int(refDF.iloc[transcriptIdx]['txEnd'])
    # print "1:",transcriptIdx,"\n",refTxSt,"\n", refTxEn, "\n"
    up_idxList = upDF[upDF['#chrom'] == chrom].index.tolist()
    # Find all overlapping transcripts
    up_recordIndex = []
    for index in (up_idxList):
        # print upDF.iloc[index],"\n"
        chromStart = upDF.iloc[index]['chromStart']
        chromEnd = upDF.iloc[index]['chromEnd']
        if((chromStart >= refTxSt) and (chromEnd <= refTxEn)):
            # print "Chr" , chromStart,chromEnd, refTxSt, refTxEn,"\n"
            if(upDF.iloc[index]['annotationType'] == 'domain'):
                up_recordIndex.append(index)
    allMaxVal = []
    allMinVal = []
    for index, val in enumerate(up_recordIndex):
        chromStart = upDF.iloc[val]['chromStart']
        chromEnd = upDF.iloc[val]['chromEnd']
        maxVal = max(refTxEn, chromEnd)
        allMaxVal.append(maxVal)
        minVal = min(refTxSt, chromStart)
        allMinVal.append(minVal)
    if(allMaxVal):
        max_len = max(allMaxVal)
    else:
        max_len = refTxEn
    if(allMinVal):
        min_len = max(allMinVal)
    else:
        min_len = refTxSt
    return(up_recordIndex, max_len, min_len)


def makeReferenceFeatures(transcript, site, chrom, pos, refDF, gds_features):
    transcripts = (refDF[refDF['#name'] == transcript])
    if(len(transcripts) > 1):
        transcriptIdx, = (transcripts[transcripts['chrom'] == chrom].index)
    else:
        transcriptIdx, = (refDF[refDF['#name'] == transcript].index)
    refTxSt = int(refDF.iloc[transcriptIdx]['txStart'])
    refTxEn = int(refDF.iloc[transcriptIdx]['txEnd'])
    # print "2:",transcriptIdx,"\n",refTxSt,"\n", refTxEn, "\n"
    ExonSts = filter(None, refDF.iloc[transcriptIdx]['exonStarts'].split(","))
    ExonEnds = filter(None, refDF.iloc[transcriptIdx]['exonEnds'].split(","))
    #ExonCounts = int(refDF.iloc[transcriptIdx]['exonCount'])
    transcriptStrand = str(refDF.iloc[transcriptIdx]['strand'])
    if(transcriptStrand == "+"):
        transcriptStrand = +1
    if(transcriptStrand == "-"):
        transcriptStrand = -1
    for idx, val in enumerate(ExonSts):
        feature = SeqFeature(
            FeatureLocation(
                int(val),
                int(ExonEnds[idx]),
                strand=transcriptStrand))
        if(transcriptStrand == -1):
            fname = "exon" + str(len(ExonSts)-idx)
            gds_features.add_feature(
                feature,
                sigil="ARROW",
                color=brown,
                arrowshaft_height=1.0,
                name=fname,
                label=True, label_position="middle", label_size=5, label_angle=90)
        else:
            fname = "exon" + str(idx + 1)
            gds_features.add_feature(
                feature,
                sigil="ARROW",
                color=brown,
                arrowshaft_height=1.0,
                name=fname,
                label=True, label_position="middle", label_size=5)
    feature = SeqFeature(FeatureLocation(pos - 5, pos + 5))
    bname = site
    gds_features.add_feature(
        feature,
        color=orange,
        name=bname,
        label=True,
        label_size=6, label_color=orange)
    return(gds_features)


def makeUniProtFeatures(domainIdx, upDF, gds_features):
    #colors = ["green", "purple", "blue", "brown", "teal", "red", "yellow"]
    for index, val in enumerate(domainIdx):
        chromStart = upDF.iloc[val]['chromStart']
        chromEnd = upDF.iloc[val]['chromEnd']
        fname = upDF.iloc[val]['name']
        feature = SeqFeature(FeatureLocation(chromStart, chromEnd), strand=None)
        gds_features.add_feature(
            feature,
            name=fname,
            label=True,
            color=green,
            label_position="middle",
            label_size=6,
            label_color=green)
    return(gds_features)


def makeReadFeatures(chrom, pos, strand, gds_features):
    start = int(pos) - 1000
    end = int(pos) + 1000
    bname = str(chrom) + ":" + str(pos)
    if(strand == 0):
        strandDirection = +1
        feature = SeqFeature(FeatureLocation(start, end, strand=strandDirection))
        gds_features.add_feature(
            feature,
            sigil="ARROW",
            arrowshaft_height=0.1,
            color=blue,
            name=bname, label_size=8, label=True, label_angle=0, label_color=purple)
    if(strand == 1):
        strandDirection = -1
        feature = SeqFeature(FeatureLocation(start, end, strand=strandDirection))
        gds_features.add_feature(
            feature,
            sigil="ARROW",
            arrowshaft_height=0.1,
            color=red,
            name=bname,
            label_size=8,
            label=True,
            label_position="middle",
            label_angle=-90,
            label_color=purple)

    return(gds_features)
# Work In Progress


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
    transcript1Idx, = (refDF[refDF['#name'] == transcript1].index)
    ExonSts1 = filter(None, refDF.iloc[transcript1Idx]['exonStarts'].split(","))
    ExonEnds1 = filter(None, refDF.iloc[transcript1Idx]['exonEnds'].split(","))
    ExonCounts1 = int(refDF.iloc[transcript1Idx]['exonCount'])
    transcript1Strand = str(refDF.iloc[transcript1Idx]['strand'])
    transcript2Idx, = (refDF[refDF['#name'] == transcript2].index)
    Exon1Sts2 = filter(None, refDF.iloc[transcript2Idx]['exonStarts'].split(","))
    Exon1Ends2 = filter(None, refDF.iloc[transcript2Idx]['exonEnds'].split(","))
    ExonCounts2 = int(refDF.iloc[transcript2Idx]['exonCount'])
    transcript2Strand = str(refDF.iloc[transcript2Idx]['strand'])
    if("before" in site1):
        before_exonNum1 = site1[-1:]
    if("after" in site1):
        after_exonNum1 = site1[-1:]
    if("before" in site2):
        before_exonNum2 = site2[-1:]
    if("after" in site2):
        after_exonNum2 = site2[-1:]
    beforeExons1 = []
    afterExons1 = []
    beforeExons2 = []
    afterExons2 = []
    for i in range(1, ExonCounts1):
        if(before_exonNum1):
            if(i <= int(before_exonNum1)):
                beforeExons1.append("exon" + i)
        if(after_exonNum1):
            if(i >= int(after_exonNum1)):
                afterExons1.append("exon" + i)
    for i in range(1, ExonCounts2):
        if(before_exonNum2):
            if(i <= int(before_exonNum2)):
                beforeExons2.append("exon" + i)
        if(after_exonNum2):
            if(i >= int(after_exonNum2)):
                afterExons2.append("exon" + i)

def get_concat_v_resize(im1, im2, resample=Image.BICUBIC, resize_big_image=True):
    if im1.width == im2.width:
        _im1 = im1
        _im2 = im2
    elif (((im1.width > im2.width) and resize_big_image) or
          ((im1.width < im2.width) and not resize_big_image)):
        _im1 = im1.resize((im2.width, int(im1.height * im2.width / im1.width)), resample=resample)
        _im2 = im2
    else:
        _im1 = im1
        _im2 = im2.resize((im1.width, int(im2.height * im1.width / im2.width)), resample=resample)
    dst = Image.new('RGB', (_im1.width, _im1.height + _im2.height))
    dst.paste(_im1, (0, 0))
    dst.paste(_im2, (0, _im1.height))
    return dst