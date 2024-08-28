"""
Created on Dec 30, 2014.

@author: Ronak H Shah

"""


def is_strand_match(txstrand, readstrand):
    """
    The function `is_strand_match` checks if the transcription strand matches the read strand based on
    specific conditions.
    
    :param txstrand: The `txstrand` parameter in the `is_strand_match` function represents the strand of
    a transcript, which can be either "+" (plus strand) or "-" (minus strand)
    :param readstrand: It seems like there might be a mistake in the code you provided. The `readstrand`
    parameter is being compared to integers (0 and 1) instead of strings ("0" and "1"). If you want to
    compare `readstrand` as strings, you should update the comparison accordingly
    :return: The function `is_strand_match` is returning a boolean value based on the conditions
    provided. It will return `True` if `txstrand` is "+" and `readstrand` is 0, or if `txstrand` is "-"
    and `readstrand` is 1. Otherwise, it will return `False`.
    """
    return (txstrand == "+" and readstrand == 0) or (
        txstrand == "-" and readstrand == 1
    )


def is_tx_active(zone):
    """
    The function `is_tx_active(zone)` checks if a given zone is active for transactions by ensuring it
    is greater than 0 and not equal to 3.
    
    :param zone: The `is_tx_active` function takes a parameter `zone` and checks if it is greater than 0
    and not equal to 3. If both conditions are met, the function returns `True`, indicating that the
    transaction is active in that zone
    :return: True if the zone is greater than 0 and not equal to 3, otherwise it returns False.
    """
    return zone > 0 and zone != 3


def determine_typeclass(ann1S, ann2S):
    """
    The function `determine_typeclass` determines the type of genetic variation based on the read
    strands of two annotations.
    
    :param ann1S: It seems like you were about to provide more information about the `ann1S` parameter
    but it got cut off. Could you please provide the complete information about the `ann1S` parameter so
    that I can assist you further with the `determine_typeclass` function?
    :param ann2S: It seems like you have not provided the complete information about the `ann2S`
    parameter. Could you please provide more details or specify what information you need help with
    regarding `ann2S`?
    :return: The function `determine_typeclass` is returning a string indicating the type of genetic
    variation based on the values of `readstrand1` and `readstrand2` attributes of the input objects
    `ann1S` and `ann2S`. If `readstrand1` is 0 and `readstrand2` is 1, it returns "Deletion". If
    `readstrand1
    """
    if ann1S.readstrand1 == 0 and ann2S.readstrand2 == 1:
        return "Deletion"
    elif ann1S.readstrand1 == 1 and ann2S.readstrand2 == 0:
        return "Duplication"
    else:
        return "Inversion"


def determine_fusion_name(ann1S, ann2S, strandmatch1):
    """
    The function `determine_fusion_name` returns a fusion name based on input gene annotations and a
    boolean indicating if the strands match.
    
    :param ann1S: It seems like you were about to provide some information about the `ann1S` parameter,
    but it got cut off. Could you please provide more details about the `ann1S` parameter so that I can
    assist you further with the `determine_fusion_name` function?
    :param ann2S: It seems like you were about to provide more information about the `ann2S` parameter
    but it got cut off. Could you please provide more details or context about the `ann2S` parameter so
    that I can assist you further with the `determine_fusion_name` function?
    :param strandmatch1: It seems like you were about to provide a description of the `strandmatch1`
    parameter, but the description is missing. Could you please provide more information about what the
    `strandmatch1` parameter represents or how it is used in the `determine_fusion_name` function?
    :return: The function `determine_fusion_name` returns a fusion name based on the input parameters
    `ann1S`, `ann2S`, and `strandmatch1`. If `strandmatch1` is `True`, it returns the fusion name in the
    format `{ann1S.gene1}:{ann2S.gene2}`. If `strandmatch1` is `False`, it returns
    """
    if strandmatch1:
        return f"{ann1S.gene1}:{ann2S.gene2}"
    else:
        return f"{ann2S.gene2}:{ann1S.gene1}"


def PredictFunctionForSV(ann1S, ann2S):
    """
    The function `PredictFunctionForSV` analyzes structural variations in genetic data to predict
    potential fusion events based on strand matching, transcript activity, gene relationships, and
    structural characteristics.
    
    :param ann1S: It seems like you were about to provide some information about the `ann1S` parameter
    in the `PredictFunctionForSV` function. Could you please provide more details or specify what
    information you need regarding the `ann1S` parameter?
    :param ann2S: It seems like you were about to provide some information about the `ann2S` parameter,
    but the details are missing. Could you please provide the necessary information about `ann2S` so
    that I can assist you further with the `PredictFunctionForSV` function?
    :return: The function `PredictFunctionForSV` returns a string indicating the type of fusion event
    based on the input annotations `ann1S` and `ann2S`. The possible return values are:
    - "Antisense Fusion" if both transcripts are active and have matching strands
    - Various descriptions of fusion events if at least one transcript is active and has a matching
    strand
    - "-" if none of
    """
    strandmatch1 = is_strand_match(ann1S.txstrand1, ann1S.readstrand1)
    strandmatch2 = is_strand_match(ann2S.txstrand2, ann2S.readstrand2)
    txactive1 = is_tx_active(ann1S.zone1)
    txactive2 = is_tx_active(ann2S.zone2)

    if txactive1 and strandmatch1 and txactive2 and strandmatch2:
        return "Antisense Fusion"

    if txactive1 and strandmatch1 or txactive2 and strandmatch2:
        if txactive1 and txactive2:
            if ann1S.gene1 == ann2S.gene2:
                typeclass = determine_typeclass(ann1S, ann2S)
                if ann1S.zone1 == 2 and ann2S.zone2 == 2:
                    if ann1S.intronnum1 == ann2S.intronnum2:
                        return f"{typeclass} within intron "
                    numExons = abs(ann2S.intronnum2 - ann1S.intronnum1)
                    txt = f"{typeclass} of {numExons} exon{'s' if numExons > 1 else ''}"
                    if typeclass != "Inversion":
                        txt += f" : {'in frame' if ann1S.intronframe1 == ann2S.intronframe2 else 'out of frame'}"
                    return txt
                txt = f"{typeclass} within transcript"
                if ann1S.zone1 == 1 or ann2S.zone2 == 1:
                    txt += " : mid-exon"
                return txt
            fusname = determine_fusion_name(ann1S, ann2S, strandmatch1)
            if ann1S.zone1 == 2 and ann2S.zone2 == 2:
                return f"Protein Fusion: {'in frame ' if ann1S.intronframe1 == ann2S.intronframe2 else 'out of frame '}{{{fusname}}}"
            if ann1S.zone1 == 1 or ann2S.zone2 == 1:
                return f"Protein Fusion: mid-exon {{{fusname}}}"
            return f"Transcript Fusion {{{fusname}}}"

    return "-"
