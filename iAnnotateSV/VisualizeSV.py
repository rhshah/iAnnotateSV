"""
Created on Apr 15, 2015
Description: This will help to plot SV using Genome Diagram from bio python
@author: Ronak H Shah
"""
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from reportlab.lib.colors import (
    red,
    grey,
    orange,
    green,
    brown,
    blue,
    lightblue,
    purple,
)
import os
import logging
from PIL import Image
import coloredlogs

coloredlogs.install(level="DEBUG")


def VisualizeSV(svDF, refDF, upDF, args):
    """
    The function VisualizeSV generates visualizations of structural variations based on input data.

    :param svDF: The `svDF` parameter is likely a DataFrame containing structural variant data, such as
    information about chromosomal positions, genes involved, fusion events, and other relevant details.
    The function `VisualizeSV` appears to iterate over each row of this DataFrame to process and
    visualize the structural variants
    :param refDF: The `refDF` parameter in the `VisualizeSV` function likely refers to a DataFrame
    containing reference data for the genomic regions being analyzed. This data could include
    information such as gene annotations, chromosome lengths, and other genomic features that are
    relevant for visualizing the structural variations (SVs) identified
    :param upDF: It seems like the code snippet you provided is incomplete. You mentioned the `upDF`
    parameter, but its definition or usage is missing from the code snippet. If you provide more context
    or details about the `upDF` parameter, I can assist you further with understanding or modifying the
    code
    :param args: The `args` parameter in the `VisualizeSV` function seems to be an object that contains
    various settings or configurations for the visualization process. It likely includes information
    such as the output file prefix, the output directory, verbosity settings (like `verbose`), and
    possibly other parameters needed for the visualization
    """
    staticDir = f"{args.outFilePrefix}_iAnnotateSVplots"
    AnalysisDir = os.path.join(args.outDir, staticDir)
    os.makedirs(AnalysisDir, exist_ok=True)
    if args.verbose:
        logging.warn(
            "iAnnotateSV::VisualizeSV: Dir: %s exists thus we wont be making it. Thus Results would be over-written",
            AnalysisDir,
        )

    for count, row in svDF.iterrows():
        chr1, chr2 = (
            format_chromosome(row.loc["chr1"]),
            format_chromosome(row.loc["chr2"]),
        )
        pos1, pos2 = int(row.loc["pos1"]), int(row.loc["pos2"])
        str1, str2 = int(row.loc["str1"]), int(row.loc["str2"])
        gene1, gene2 = str(row.loc["gene1"]), str(row.loc["gene2"])
        site1, site2 = str(row.loc["site1"]), str(row.loc["site2"])
        transcript1, transcript2 = (
            str(row.loc["transcript1"]),
            str(row.loc["transcript2"]),
        )
        fusion = str(row.loc["fusion"])

        if fusion != "-":
            domain1Idx, maxLen1, minLen1 = processData(chr1, transcript1, refDF, upDF)
            domain2Idx, maxLen2, minLen2 = processData(chr2, transcript2, refDF, upDF)
            eventtype = determine_event_type(str1, str2, chr1, chr2)

            outFileName = os.path.join(
                AnalysisDir,
                f"{gene1}-{chr1}_{pos1}_{gene2}-{chr2}_{pos2}_{eventtype}.jpg",
            )
            if gene1 != gene2:
                outFile1Name = outFileName.replace(".jpg", "-part1.jpg")
                outFile2Name = outFileName.replace(".jpg", "-part2.jpg")
                gdd1, gdd2 = create_diagrams(eventtype, gene1, gene2)
                create_tracks(
                    gdd1,
                    gdd2,
                    gene1,
                    gene2,
                    eventtype,
                    transcript1,
                    site1,
                    chr1,
                    pos1,
                    refDF,
                    str1,
                    domain1Idx,
                    upDF,
                    transcript2,
                    site2,
                    chr2,
                    pos2,
                    str2,
                    domain2Idx,
                )
                draw_and_save(gdd1, outFile1Name, minLen1, maxLen1)
                draw_and_save(gdd2, outFile2Name, minLen2, maxLen2)
                merge_images(outFile1Name, outFile2Name, outFileName)
            else:
                gdd = GenomeDiagram.Diagram("Test Diagram")
                create_tracks(
                    gdd,
                    gdd,
                    gene1,
                    gene2,
                    eventtype,
                    transcript1,
                    site1,
                    chr1,
                    pos1,
                    refDF,
                    str1,
                    domain1Idx,
                    upDF,
                    transcript2,
                    site2,
                    chr2,
                    pos2,
                    str2,
                    domain2Idx,
                )
                draw_and_save(
                    gdd, outFileName, min(minLen1, minLen2), max(maxLen1, maxLen2)
                )


def format_chromosome(chrom):
    """
    The `format_chromosome` function adds "chr" prefix to the chromosome number if it's not already
    present.

    :param chrom: The `format_chromosome` function takes a chromosome identifier as input and returns it
    formatted with the prefix "chr" if it doesn't already start with "chr"
    :return: The `format_chromosome` function takes a chromosome identifier as input and returns it with
    the prefix "chr" if it doesn't already start with "chr".
    """
    return f"chr{chrom}" if not str(chrom).startswith("chr") else str(chrom)


def determine_event_type(str1, str2, chr1, chr2):
    """
    The function `determine_event_type` categorizes genetic events as inversion, duplication, deletion,
    or translocation based on input parameters.

    :param str1: It seems like you were about to provide some information about the parameters `str1`,
    `str2`, `chr1`, and `chr2` for the function `determine_event_type`. Could you please continue with
    the information you wanted to share about `str1`?
    :param str2: It seems like you were about to provide more information about the parameters, but the
    input got cut off. Could you please provide the complete information about the parameters str2,
    chr1, and chr2 so I can assist you further?
    :param chr1: It seems like you were about to provide the explanation for the parameter `chr1` in the
    `determine_event_type` function. Could you please provide more details or complete the information
    you intended to share about `chr1` so that I can assist you further?
    :param chr2: It seems like you have provided the function `determine_event_type` that takes four
    parameters `str1`, `str2`, `chr1`, and `chr2` to determine the type of genetic event based on their
    values. However, you have not provided the value for `chr2`
    :return: The function `determine_event_type` returns the type of genetic event based on the input
    parameters `str1`, `str2`, `chr1`, and `chr2`. The possible return values are "Inversion",
    "Duplication", "Deletion", or "Translocation" depending on the conditions specified in the function.
    """
    if str1 == str2 == 0 and chr1 == chr2:
        return "Inversion"
    elif str1 == str2 == 1 and chr1 == chr2:
        return "Inversion"
    elif str1 == 1 and str2 == 0 and chr1 == chr2:
        return "Duplication"
    elif str1 == 0 and str2 == 1 and chr1 == chr2:
        return "Deletion"
    else:
        return "Translocation"


def create_diagrams(eventtype, gene1, gene2):
    """
    The function create_diagrams takes event type, gene1, and gene2 as input and returns two
    GenomeDiagram.Diagram objects named based on the event type and genes.

    :param eventtype: The `eventtype` parameter is a string that represents the type of event for which
    the diagrams are being created. It could be something like "mutation", "translocation", "deletion",
    etc
    :param gene1: It looks like you were about to provide more information about the parameters, but it
    seems like the message got cut off. Could you please provide more details about the `gene1`
    parameter so that I can assist you further?
    :param gene2: The `gene2` parameter is likely the name or identifier of a second gene that is
    related to the event type specified in the `eventtype` parameter. This function `create_diagrams`
    seems to create two GenomeDiagram objects for visualizing the specified event type with respect to
    two different genes (
    :return: Two GenomeDiagram.Diagram objects named gdd1 and gdd2, each with a specific title based on
    the event type and gene names provided as input.
    """
    gdd1 = GenomeDiagram.Diagram(f"{eventtype}-{gene1}")
    gdd2 = GenomeDiagram.Diagram(f"{eventtype}-{gene2}")
    return gdd1, gdd2


def create_tracks(
    gdd1,
    gdd2,
    gene1,
    gene2,
    eventtype,
    transcript1,
    site1,
    chr1,
    pos1,
    refDF,
    str1,
    domain1Idx,
    upDF,
    transcript2,
    site2,
    chr2,
    pos2,
    str2,
    domain2Idx,
):
    """
    The function `create_tracks` generates tracks for gene features, alignment information, and domain
    information for two genes in a genome browser display.

    :param gdd1: The `gdd1` parameter in the `create_tracks` function seems to be a reference to a
    genomic data display object or a similar data structure used for visualization or analysis of
    genomic data. It is used to create new tracks for displaying various features related to genes,
    alignments, and domains for a
    :param gdd2: It seems like the description of the `gdd2` parameter is missing. Could you please
    provide more information or context about what `gdd2` represents or is used for in the
    `create_tracks` function?
    :param gene1: `gene1` is a variable representing the name of the first gene involved in the event
    for which tracks are being created
    :param gene2: The `gene2` parameter in the `create_tracks` function represents the second gene
    involved in the event for which tracks are being created. It is used to generate track names and
    features related to this gene in the visualization tracks
    :param eventtype: The `eventtype` parameter in the `create_tracks` function represents the type of
    event being analyzed or tracked. It could be a specific biological event or process that is being
    studied in the context of the provided genes (`gene1` and `gene2`). This parameter helps in
    generating track names and
    :param transcript1: The `transcript1` parameter in the `create_tracks` function represents the
    transcript associated with gene1. It is used to specify the transcript information for creating
    features in the genomic data display tracks
    :param site1: It seems like the description of the `site1` parameter is missing in the provided code
    snippet. Could you please provide more information or context about what `site1` represents or how
    it is used within the `create_tracks` function? This will help me provide you with a more accurate
    explanation or
    :param chr1: The parameter `chr1` in the `create_tracks` function likely represents the chromosome
    number or identifier for the genomic location associated with `gene1` and `transcript1`. It is used
    in various parts of the function to specify the chromosome for creating different tracks and
    features related to the genomic data
    :param pos1: The `pos1` parameter likely represents the position of a genomic feature in the genome.
    It could be a numerical value indicating the specific location of a gene, mutation, or other genetic
    element on a chromosome
    :param refDF: The `refDF` parameter in the `create_tracks` function seems to be a DataFrame that is
    used for reference features. It is passed as an argument to the `makeReferenceFeatures` function
    along with other parameters like `transcript1`, `site1`, `chr1`, and `pos
    :param str1: The parameter `str1` in the `create_tracks` function appears to represent the strand
    information for the first gene or transcript. It is likely used to indicate the directionality of
    the gene or transcript on the DNA strand (e.g., whether it is located on the forward (+) or reverse
    (-
    :param domain1Idx: The `domain1Idx` parameter is used to specify the index of a domain for gene 1 in
    the `makeUniProtFeatures` function. If `domain1Idx` has a value (i.e., not None or False), it will
    be used to retrieve information about the domain from the
    :param upDF: The `upDF` parameter in the `create_tracks` function seems to represent a DataFrame
    containing UniProt data. This data is likely used to create UniProt features for the tracks being
    generated in the function. If you need further assistance or have any specific questions regarding
    the `upDF` parameter or
    :param transcript2: The `transcript2` parameter in the `create_tracks` function likely represents
    the transcript associated with the second gene (`gene2`) for which tracks are being created. This
    parameter is used to specify the transcript information needed to generate features for the second
    gene in the visualization tracks
    :param site2: The `site2` parameter likely refers to the site or location associated with
    `transcript2` in the `create_tracks` function. It is used as input to create certain features or
    tracks related to the second gene and event type being processed in the function. The specific
    details of how `site
    :param chr2: The parameter `chr2` in the `create_tracks` function likely represents the chromosome
    number or identifier for the second gene or genomic location being processed. It is used as part of
    the genomic coordinates to locate specific features within the genome
    :param pos2: The `pos2` parameter in the `create_tracks` function likely represents the position of
    a genomic feature on a chromosome for the second gene (`gene2`). This parameter is used in the
    function to create tracks and features related to the second gene's genomic data
    :param str2: The parameter `str2` in the `create_tracks` function appears to represent the strand
    information for the second gene or feature being processed. The strand information typically
    indicates the directionality of the gene or feature on the DNA molecule, with values like '1' or '+'
    representing the forward strand and values
    :param domain2Idx: The `domain2Idx` parameter is likely used to specify the index of a domain in the
    context of the second gene or feature being processed in the `create_tracks` function. This index
    may be used to retrieve specific information or features related to that domain from the `upDF`
    (UniProt
    """
    feature_names = [
        f"GeneTrack:{gene1}:{eventtype}",
        f"AlignmentTrack:{gene1}:{eventtype}",
        f"DomainTrack:{gene1}:{eventtype}",
        f"GeneTrack:{gene2}:{eventtype}",
        f"AlignmentTrack:{gene2}:{eventtype}",
        f"DomainTrack:{gene2}:{eventtype}",
    ]
    tracks1 = [
        gdd1.new_track(1, greytrack=True, name=name) for name in feature_names[:3]
    ]
    tracks2 = [
        gdd2.new_track(1, greytrack=True, name=name) for name in feature_names[3:]
    ]

    gds_features = tracks1[0].new_set()
    gds_features = makeReferenceFeatures(
        transcript1, site1, chr1, pos1, refDF, gds_features
    )
    gds_features = tracks1[1].new_set()
    gds_features = makeReadFeatures(chr1, pos1, str1, gds_features)
    gds_features = tracks1[2].new_set()
    if domain1Idx:
        gds_features = makeUniProtFeatures(domain1Idx, upDF, gds_features)

    gds_features = tracks2[0].new_set()
    gds_features = makeReferenceFeatures(
        transcript2, site2, chr2, pos2, refDF, gds_features
    )
    gds_features = tracks2[1].new_set()
    gds_features = makeReadFeatures(chr2, pos2, str2, gds_features)
    gds_features = tracks2[2].new_set()
    if domain2Idx:
        gds_features = makeUniProtFeatures(domain2Idx, upDF, gds_features)


def draw_and_save(gdd, outFileName, minLen, maxLen):
    """
    The function `draw_and_save` takes a GenomicDiagram object, an output file name, minimum and maximum
    lengths as input, draws a linear format diagram with specified fragments and saves it as a JPG image
    with a DPI of 300.

    :param gdd: The `gdd` parameter seems to be an object that likely represents some kind of graphical
    display or diagram. The function `draw_and_save` takes this object, an output file name
    (`outFileName`), a minimum length (`minLen`), and a maximum length (`maxLen`) as parameters
    :param outFileName: The `outFileName` parameter in the `draw_and_save` function is the name of the
    file where the graph will be saved. You should provide a string value representing the file name
    along with the desired file extension (e.g., "graph.jpg")
    :param minLen: The `minLen` parameter in the `draw_and_save` function represents the minimum length
    of the sequence to be drawn. This value is used to set the start position for drawing the sequence
    on the graphical display
    :param maxLen: The `maxLen` parameter in the `draw_and_save` function likely represents the maximum
    length of the sequence or data that will be visualized and saved as an image. This value is used to
    determine the end point for the visualization
    """
    gdd.draw(format="linear", fragments=1, start=minLen - 1000, end=maxLen + 1000)
    gdd.write(outFileName, "JPG", dpi=300)


def merge_images(outFile1Name, outFile2Name, outFileName):
    """
    The function `merge_images` takes two image files, merges them vertically, saves the result as a new
    image file, and then deletes the original two image files.

    :param outFile1Name: The `outFile1Name` parameter in the `merge_images` function is the file name of
    the first image that you want to merge with another image
    :param outFile2Name: The `outFile2Name` parameter in the `merge_images` function is the file name of
    the second image that you want to merge with another image
    :param outFileName: The `outFileName` parameter in the `merge_images` function represents the file
    name under which the merged image will be saved after combining `outFile1Name` and `outFile2Name`
    """
    img1 = Image.open(outFile1Name)
    img2 = Image.open(outFile2Name)
    images = [img1, img2]
    w = max(i.size[0] for i in images)
    mh = sum(i.size[1] for i in images)
    result = Image.new("RGB", (w, mh), (255, 255, 255))
    x = 0
    for i in images:
        result.paste(i, (0, x))
        x += i.size[1]
    result.save(outFileName)
    os.remove(outFile1Name)
    os.remove(outFile2Name)


def processData(chrom, transcript, refDF, upDF):
    """
    The function processData processes data related to a specific transcript and chromosome, extracting
    relevant information from two input dataframes.

    :param chrom: Chrom is a variable that typically represents the chromosome number or identifier in a
    genetic sequence. It is used to specify the chromosome on which a particular genetic operation or
    analysis is being performed
    :param transcript: Transcript is a variable that represents the name of a specific transcript. It is
    used to filter data from the reference DataFrame (refDF) to find information related to that
    particular transcript
    :param refDF: refDF is a DataFrame containing reference data for transcripts, with columns such as
    "#name", "chrom", "txStart", and "txEnd". It seems to store information about the genomic locations
    of transcripts
    :param upDF: The `upDF` parameter in the `processData` function seems to be a DataFrame containing
    information about genomic regions. The function processes this DataFrame based on the provided
    `chrom` and `transcript` parameters along with reference DataFrames `refDF`
    :return: The function `processData` returns three values:
    1. `up_recordIndex`: a list of indices of records in the `upDF` DataFrame that meet certain
    conditions
    2. `max_len`: the maximum value between the end position of the transcript and the end position of
    the records in `upDF`
    3. `min_len`: the minimum value between the start position of the transcript
    """
    transcripts = refDF[refDF["#name"] == transcript]
    transcriptIdx = (
        transcripts[transcripts["chrom"] == chrom].index[0]
        if len(transcripts) > 1
        else refDF[refDF["#name"] == transcript].index[0]
    )
    refTxSt, refTxEn = (
        int(refDF.iloc[transcriptIdx]["txStart"]),
        int(refDF.iloc[transcriptIdx]["txEnd"]),
    )
    up_idxList = upDF[upDF["#chrom"] == chrom].index.tolist()
    up_recordIndex = [
        index
        for index in up_idxList
        if refTxSt <= upDF.iloc[index]["chromStart"] <= refTxEn
        and upDF.iloc[index]["annotationType"] == "domain"
    ]
    allMaxVal = [max(refTxEn, upDF.iloc[val]["chromEnd"]) for val in up_recordIndex]
    allMinVal = [min(refTxSt, upDF.iloc[val]["chromStart"]) for val in up_recordIndex]
    max_len = max(allMaxVal) if allMaxVal else refTxEn
    min_len = min(allMinVal) if allMinVal else refTxSt
    return up_recordIndex, max_len, min_len


def makeReferenceFeatures(transcript, site, chrom, pos, refDF, gds_features):
    """
    The function `makeReferenceFeatures` generates genomic features for a given transcript and site
    based on reference data and adds them to a graphical display.

    :param transcript: Transcript is a string representing the name of the transcript for which
    reference features are being generated
    :param site: The `site` parameter in the `makeReferenceFeatures` function represents the specific
    genomic site or position for which you want to create reference features. This function is designed
    to generate reference features for a given transcript and genomic site, including exon features and
    a specific site feature indicated by the `site` parameter
    :param chrom: The `chrom` parameter in the `makeReferenceFeatures` function represents the
    chromosome where the genetic variant is located. It is used to filter the reference transcripts
    based on the chromosome information provided
    :param pos: The `pos` parameter in the `makeReferenceFeatures` function represents the genomic
    position where a specific site or feature is located. This position is used to create a sequence
    feature that is highlighted in the visualization generated by the function
    :param refDF: The `refDF` parameter is a DataFrame containing reference information for transcripts.
    It likely includes columns such as `#name`, `chrom`, `txStart`, `txEnd`, `exonStarts`, `exonEnds`,
    and `strand` among others. This DataFrame is used to
    :param gds_features: The `gds_features` parameter is likely an object that represents a graphical
    display of sequence features. The function `makeReferenceFeatures` seems to be adding features to
    this graphical display based on the provided transcript information, genomic site, and reference
    data frame
    :return: The function `makeReferenceFeatures` is returning the `gds_features` object after adding
    features representing exons and a specific site within a transcript.
    """
    transcripts = refDF[refDF["#name"] == transcript]
    transcriptIdx = (
        transcripts[transcripts["chrom"] == chrom].index[0]
        if len(transcripts) > 1
        else refDF[refDF["#name"] == transcript].index[0]
    )
    refTxSt, refTxEn = (
        int(refDF.iloc[transcriptIdx]["txStart"]),
        int(refDF.iloc[transcriptIdx]["txEnd"]),
    )
    ExonSts = list(filter(None, refDF.iloc[transcriptIdx]["exonStarts"].split(",")))
    ExonEnds = list(filter(None, refDF.iloc[transcriptIdx]["exonEnds"].split(",")))
    transcriptStrand = 1 if refDF.iloc[transcriptIdx]["strand"] == "+" else -1
    for idx, val in enumerate(ExonSts):
        feature = SeqFeature(
            FeatureLocation(int(val), int(ExonEnds[idx]), strand=transcriptStrand)
        )
        fname = (
            f"exon{len(ExonSts)-idx}" if transcriptStrand == -1 else f"exon{idx + 1}"
        )
        gds_features.add_feature(
            feature,
            sigil="ARROW",
            color=brown,
            arrowshaft_height=1.0,
            name=fname,
            label=True,
            label_position="middle",
            label_size=5,
            label_angle=90 if transcriptStrand == -1 else 0,
        )
    feature = SeqFeature(FeatureLocation(pos - 5, pos + 5))
    gds_features.add_feature(
        feature, color=orange, name=site, label=True, label_size=6, label_color=orange
    )
    return gds_features


def makeUniProtFeatures(domainIdx, upDF, gds_features):
    """
    The function `makeUniProtFeatures` adds features to a GenBank Feature Table using information from a
    UniProt DataFrame and a list of domain indices.

    :param domainIdx: DomainIdx is a list of indices representing the positions of domains in a dataset
    :param upDF: The `upDF` parameter seems to be a DataFrame containing information about UniProt
    entries. The function `makeUniProtFeatures` takes `domainIdx`, which is a list of indices, `upDF`,
    and `gds_features` as input parameters
    :param gds_features: `gds_features` seems to be an object that likely represents genomic features.
    The function `makeUniProtFeatures` takes `domainIdx`, `upDF`, and `gds_features` as parameters. It
    iterates over the indices in `domainIdx`, extracts relevant information from `upDF
    :return: The function `makeUniProtFeatures` is returning the `gds_features` object after adding
    features based on the `domainIdx` values from the `upDF` DataFrame.
    """
    for val in domainIdx:
        chromStart, chromEnd = upDF.iloc[val]["chromStart"], upDF.iloc[val]["chromEnd"]
        fname = upDF.iloc[val]["name"]
        feature = SeqFeature(FeatureLocation(chromStart, chromEnd), strand=None)
        gds_features.add_feature(
            feature,
            name=fname,
            label=True,
            color=green,
            label_position="middle",
            label_size=6,
            label_color=green,
        )
    return gds_features


def makeReadFeatures(chrom, pos, strand, gds_features):
    """
    The function `makeReadFeatures` creates a SeqFeature object and adds it to a set of genomic features
    with specific visual attributes based on input parameters.

    :param chrom: The `chrom` parameter represents the chromosome on which the feature is located
    :param pos: The `pos` parameter represents the position on the chromosome where a feature will be
    created
    :param strand: The `strand` parameter in the `makeReadFeatures` function represents the orientation
    of the feature on the DNA strand. It can have two possible values: 0 or 1. A value of 0 indicates
    that the feature is on the forward strand, while a value of 1 indicates that
    :param gds_features: The `gds_features` parameter is likely an object or data structure that stores
    genomic features. The `makeReadFeatures` function takes the input parameters `chrom`, `pos`,
    `strand`, and `gds_features`, and creates a new genomic feature based on these inputs. The function
    then adds
    :return: The function `makeReadFeatures` is returning the `gds_features` object after adding a new
    feature to it based on the input parameters `chrom`, `pos`, `strand`, and `gds_features`.
    """
    start, end = int(pos) - 1000, int(pos) + 1000
    bname = f"{chrom}:{pos}"
    strandDirection = +1 if strand == 0 else -1
    feature = SeqFeature(FeatureLocation(start, end, strand=strandDirection))
    gds_features.add_feature(
        feature,
        sigil="ARROW",
        arrowshaft_height=0.1,
        color=blue if strand == 0 else red,
        name=bname,
        label_size=8,
        label=True,
        label_angle=0 if strand == 0 else -90,
        label_color=purple,
    )
    return gds_features
