'''
Created on 01/09/2018
@Ronak Shah

'''

import contextlib
import os
import sys
import polars as pl
import logging
from rich import print
from rich.logging import RichHandler
import re
import helper as hp

FORMAT = "%(message)s"
logging.basicConfig(
    level="INFO", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)

log = logging.getLogger("rich")


def run(svDFA, refPath, ctPath, allctPath, upPath, verbose):
    """
    Annotates structural variants with kinase domain information.

    Args:
        svDFA (pl.DataFrame): DataFrame containing structural variant annotations.
        refPath (str): Path to the reference annotation file.
        ctPath (str): Path to the canonical transcript file.
        allctPath (str): Path to the all canonical transcripts file.
        upPath (str): Path to the UniProt annotation file.
        verbose (bool): Verbosity flag.

    Returns:
        pl.DataFrame: DataFrame with added kinase domain annotations.
    """

    # Load dataframes
    upDF = load_dataframe(upPath, "UniProt annotation", critical=True)
    ctDF = load_dataframe(ctPath, "Assay-specific canonical transcript", critical=False)
    allctDF = load_dataframe(allctPath, "All canonical transcripts", critical=True)
    refDF = load_dataframe(refPath, "Reference annotation", critical=True)

    # Initialize kinase domain columns
    svDF = svDFA.with_columns([
        pl.lit(None).alias('kinase_domain1'),
        pl.lit(None).alias('kinase_domain2')
    ])

    # Annotate each structural variant
    for count, row in enumerate(svDFA.iter_rows(named=True)):
        log.info(f"iAnnotateSV::AnnotateForKinaseDomain: Checking Entry {count} in Uniprot data")

        chr1, chr2 = row['chr1'], row['chr2']
        chr1 = chr1 if chr1.startswith('chr') else f"chr{chr1}"
        chr2 = chr2 if chr2.startswith('chr') else f"chr{chr2}"
        pos1, pos2 = int(row['pos1']), int(row['pos2'])
        gene1, gene2 = row['gene1'], row['gene2']
        fusion = row['fusion']

        transcript1 = get_transcript(ctDF, allctDF, gene1)
        transcript2 = get_transcript(ctDF, allctDF, gene2)

        kanno1, kanno2 = None, None

        if fusion != "-":
            if fusionevent := re.search(r'\{(.*)\}', fusion):
                eventType = fusionevent.group(1)
                if ":" in eventType:
                    egene1, egene2 = eventType.split(":")

                    if transcript1:
                        kanno1 = getKinaseInfo(chr1, pos1, gene1, egene1, egene2, transcript1, refDF, upDF)
                    if transcript2:
                        kanno2 = getKinaseInfo(chr2, pos2, gene2, egene1, egene2, transcript2, refDF, upDF)

            svDF = svDF.with_columns([
                pl.when(pl.lit(True)).then(pl.lit(kanno1)).otherwise(pl.col("kinase_domain1")).alias("kinase_domain1"),
                pl.when(pl.lit(True)).then(pl.lit(kanno2)).otherwise(pl.col("kinase_domain2")).alias("kinase_domain2")
            ])

    return svDF


def load_dataframe(path, description, critical=True):
    """
    Loads a dataframe from a file.

    Args:
        path (str): Path to the file.
        description (str): Description of the data being loaded.
        critical (bool): Whether to exit if the file is not found.

    Returns:
        pl.DataFrame: The loaded dataframe, or an empty dataframe if critical is False and the file is not found.
    """
    if os.path.isfile(path):
        return hp.ReadFile(path)
    message = f"iAnnotateSV::AnnotationForKinaseDomain: Location of {description} file is incorrect!!!"
    if critical:
        log.critical(message)
        sys.exit(1)
    else:
        log.warn(message)
        return pl.DataFrame()


def get_transcript(ctDF, allctDF, gene):
    """
    Retrieves the canonical transcript for a gene.

    Args:
        ctDF (pl.DataFrame): DataFrame containing assay-specific canonical transcripts.
        allctDF (pl.DataFrame): DataFrame containing all canonical transcripts.
        gene (str): Gene symbol.

    Returns:
        str: The canonical transcript, or None if not found.
    """
    transcript = None
    for df in [ctDF, allctDF]:
        if not df.is_empty:
            try:
                transcript = df.filter(pl.col('Gene') == gene)['Transcripts'][0]
                break  # Stop after finding the transcript in the first available dataframe
            except (IndexError, KeyError):
                continue  # Try the next dataframe
    return transcript


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
    transcripts = refDF.filter(pl.col('name') == transcript)
    if transcripts.is_empty():
        return (None, None, None)

    transcriptIdx = 0
    refTxSt = int(transcripts[transcriptIdx, 'txStart'])
    refTxEn = int(transcripts[transcriptIdx, 'txEnd'])

    # Find overlapping UniProt records
    up_idxList = upDF.filter(pl.col('#chrom') == chrom).select(pl.col('#chrom')).to_series().to_list()
    up_recordIndex = [
        index for index in up_idxList
        if (upDF[index, 'chromStart'] >= refTxSt and upDF[index, 'chromEnd'] <= refTxEn and
            upDF[index, 'annotationType'] == 'domain')
    ]

    # Determine max and min lengths
    allMaxVal = [upDF[index, 'chromEnd'] for index in up_recordIndex]
    allMinVal = [upDF[index, 'chromStart'] for index in up_recordIndex]

    max_len = max(allMaxVal + [refTxEn]) if allMaxVal else refTxEn
    min_len = min(allMinVal + [refTxSt]) if allMinVal else refTxSt

    return (up_recordIndex, max_len, min_len)


def getKinaseInfo(chrom, pos, gene, egene1, egene2, transcript, refDF, upDF):
    """
    Determines if a kinase domain is included in a structural variant.

    Args:
        chrom (str): Chromosome.
        pos (int): Breakpoint position.
        gene (str): Gene symbol.
        egene1 (str): First gene in the fusion event.
        egene2 (str): Second gene in the fusion event.
        transcript (str): Transcript ID.
        refDF (pl.DataFrame): DataFrame containing reference transcript annotations.
        upDF (pl.DataFrame): DataFrame containing UniProt annotations.

    Returns:
        str: Kinase domain inclusion status, or None if not found.
    """
    domainIdx, _, _ = processData(chrom, transcript, refDF, upDF)
    if domainIdx is None:
        return None

    strand = refDF.filter(pl.col('name') == transcript)['strand'][0]

    kanno = None
    if egene1 == gene or egene2 == gene:  # Check if either gene matches
        # Determine inclusion status based on breakpoint position and domain location
        for index in domainIdx:
            chromStart, chromEnd, fname = upDF[index, 'chromStart'], upDF[index, 'chromEnd'], upDF[index, 'name']
            if "Protein kinase" in fname:
                if ((strand == "+" and ((egene1 == gene and pos > chromEnd) or (egene2 == gene and pos < chromStart))) or
                        (strand == "-" and ((egene1 == gene and pos < chromStart) or (egene2 == gene and pos > chromEnd)))):
                    kanno = "Kinase Domain Included"
                elif chromStart <= pos <= chromEnd:
                    kanno = "Partial Kinase Domain Included"
                else:
                    kanno = "Kinase Domain Not Included"
                return kanno

    return None


def getValueOrDefault(value, index, default=None):
    returnValue = default
    with contextlib.suppress(Exception):
        returnValue = value[index]
    return returnValue