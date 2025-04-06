"""
Created on Dec 30, 2014.

@author: Ronak H Shah

"""

import polars as pl
import logging
from rich import print
from rich.logging import RichHandler

FORMAT = "%(message)s"
logging.basicConfig(
    level="INFO", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)

log = logging.getLogger("rich")


def PredictFunctionForSV(ann1S, ann2S):
    """
    Predicts the functional consequence of a structural variant (SV) based on
    the annotations of the two breakpoints.

    Args:
        ann1S (pl.Series): Annotation information for the first breakpoint.
        ann2S (pl.Series): Annotation information for the second breakpoint.

    Returns:
        str: A string describing the predicted functional consequence of the SV.
    """
    log.debug("Predicting functional consequence for SV...")
    log.debug(f"Annotation 1: {ann1S}")
    log.debug(f"Annotation 2: {ann2S}")

    strandmatch1 = ((ann1S["txstrand1"] == '+' and ann1S["readstrand1"] == 0)
                    or (ann1S["txstrand1"] == '-' and ann1S["readstrand1"] == 1))
    strandmatch2 = ((ann2S["txstrand2"] == '+' and ann2S["readstrand2"] == 0)
                    or (ann2S["txstrand2"] == '-' and ann2S["readstrand2"] == 1))
    txactive1 = (ann1S["zone1"] > 0 and ann1S["zone1"] != 3)
    txactive2 = (ann2S["zone2"] > 0 and ann2S["zone2"] != 3)

    log.debug(f"Strand match 1: {strandmatch1}")
    log.debug(f"Strand match 2: {strandmatch2}")
    log.debug(f"Transcriptional activity 1: {txactive1}")
    log.debug(f"Transcriptional activity 2: {txactive2}")

    txt = '-'
    if (txactive1 and strandmatch1) and (txactive2 and strandmatch2):
        txt = 'Antisense Fusion'
        log.info(f"Predicted functional consequence: {txt}")
    elif (txactive1 or txactive2) and (strandmatch1 or strandmatch2):
        if txactive1 and txactive2:
            if ann1S["gene1"] == ann2S["gene2"]:
                log.debug("SV within the same gene")
                if ann1S["readstrand1"] == 0 and ann2S["readstrand2"] == 1:
                    typeclass = 'Deletion'
                elif ann1S["readstrand1"] == 1 and ann2S["readstrand2"] == 0:
                    typeclass = 'Duplication'
                else:
                    typeclass = 'Inversion'
                log.debug(f"Type class: {typeclass}")
                if ann1S["zone1"] == 2 and ann2S["zone2"] == 2:
                    if ann1S["intronnum1"] == ann2S["intronnum2"]:
                        txt = f'{typeclass} within intron '
                    else:
                        numExons = abs(ann2S["intronnum2"] - ann1S["intronnum1"])
                        txt = f'{typeclass} of {numExons} exon'
                        if numExons > 1:
                            txt += 's'
                        if typeclass != 'Inversion':
                            if ann1S["intronframe1"] == ann2S["intronframe2"]:
                                txt += ' : in frame'
                            else:
                                txt += ' : out of frame'
                else:
                    txt = f'{typeclass} within transcript'
                    if ann1S["zone1"] == 1 or ann2S["zone2"] == 1:
                        txt += ' : mid-exon'
            else:
                log.debug("SV between different genes")
                fusname = (
                    f'{ann1S["gene1"]}:{ann2S["gene2"]}'
                    if strandmatch1
                    else f'{ann2S["gene2"]}:{ann1S["gene1"]}'
                )
                if ann1S["zone1"] == 2 and ann2S["zone2"] == 2:
                    txt = (
                        'Protein Fusion: in frame '
                        if ann1S["intronframe1"] == ann2S["intronframe2"]
                        else 'Protein Fusion: out of frame '
                    )
                elif ann1S["zone1"] == 1 or ann2S["zone2"] == 1:
                    txt = 'Protein Fusion: mid-exon '
                else:
                    txt = 'Transcript Fusion'
                txt += f' {{{fusname}}}'
            log.info(f"Predicted functional consequence: {txt}")
    log.debug(f"Returning predicted functional consequence: {txt}")
    return txt