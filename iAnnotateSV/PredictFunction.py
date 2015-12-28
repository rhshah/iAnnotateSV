"""
Created on Dec 30, 2014.

@author: Ronak H Shah

"""


def PredictFunctionForSV(ann1S, ann2S):
    strandmatch1 = ((ann1S.txstrand1 == '+' and ann1S.readstrand1 == 0)
                    or(ann1S.txstrand1 == '-' and ann1S.readstrand1 == 1))
    strandmatch2 = ((ann2S.txstrand2 == '+' and ann2S.readstrand2 == 0)
                    or(ann2S.txstrand2 == '-' and ann2S.readstrand2 == 1))
    txactive1 = (ann1S.zone1 > 0 and ann1S.zone1 != 3)
    txactive2 = (ann2S.zone2 > 0 and ann2S.zone2 != 3)

    txt = '-'
    #typeclass = None
    if(txactive1 and strandmatch1) and (txactive2 and strandmatch2):
        txt = 'Antisense Fusion'
    elif(txactive1 and strandmatch1) or (txactive2 and strandmatch2):
        if(txactive1 and txactive2):
            # Within Gene Even
            if(ann1S.gene1 == ann2S.gene2):
                if(ann1S.readstrand1 == 0 and ann2S.readstrand2 == 1):
                    typeclass = 'Deletion'
                elif(ann1S.readstrand1 == 1 and ann2S.readstrand2 == 0):
                    typeclass = 'Duplication'
                else:
                    typeclass = 'Inversion'
                if (ann1S.zone1 == 2 and ann2S.zone2 == 2):
                    if(ann1S.intronnum1 == ann2S.intronnum2):
                        txt = typeclass + ' within intron '
                    else:
                        numExons = abs(ann2S.intronnum2 - ann1S.intronnum1)
                        txt = typeclass + ' of ' + str(numExons) + ' exon'
                        if(numExons > 1):
                            txt = txt + 's'
                        if(typeclass != 'Inversion'):
                            if(ann1S.intronframe1 == ann2S.intronframe2):
                                txt = txt + ' : in frame'
                            else:
                                txt = txt + ' : out of frame'
                else:
                    txt = typeclass + ' within transcript'
                    if(ann1S.zone1 == 1 or ann2S.zone2 == 1):
                        txt = txt + ' : mid-exon'
            # Inter Gene Event
            else:
                if(strandmatch1):
                    fusname = ann1S.gene1 + ":" + ann2S.gene2
                else:
                    fusname = ann2S.gene2 + ":" + ann1S.gene1
                if(ann1S.zone1 == 2 and ann2S.zone2 == 2):
                    if(ann1S.intronframe1 == ann2S.intronframe2):
                        txt = 'Protein Fusion: in frame '
                    else:
                        txt = 'Protein Fusion: out of frame '
                elif(ann1S.zone1 == 1 or ann2S.zone2 == 1):
                    txt = 'Protein Fusion: mid-exon '
                else:
                    txt = 'Transcript Fusion'
                txt = txt + ' {' + fusname + '}'

    return(txt)
