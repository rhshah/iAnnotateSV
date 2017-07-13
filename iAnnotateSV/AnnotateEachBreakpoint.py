'''
Created on Dec 29, 2014

@author: Ronak H Shah
'''
import pandas as pd
import numpy as np
import helper as hp
import FindTranscript as ft

def AnnotateEachBreakpoint(chromosome,position,strand,df,autoSelect):
    #print "Annotating a coordinate:",position
    if(chromosome.startswith('chr')):
        chromosome = chromosome
    else:
        chromosome = "chr" + chromosome
    #Find all the chromosomes
    idxList = df[df['chrom'] == chromosome].index.tolist()
    transcriptIndex = []
    #Find all overlapping transcripts
    for index in (idxList):
        geneStart = df.iloc[index]['geneStart']
        geneEnd = df.iloc[index]['geneEnd']
        if((geneStart <= position) and (geneEnd >= position)):
            #print position,geneStart,geneEnd
            transcriptIndex.append(index)
    desc = None
    intronnum = None
    intronframe = None
    #print transcriptIndex
    if(transcriptIndex):
        coordData = pd.DataFrame(index=np.asarray(transcriptIndex),columns=['c', 'd', 'e', 'd1', 'd2', 'e1', 'e2','f'])
        
        #print coordData.head()
        #print transcriptIndex,coordData.index
        for tindex in (transcriptIndex):
            #Annotate coords with each transcript
            #For coords within a transcript
            c = None # zone: 1=exon, 2=intron, 3=3'-UTR, 4=5'-UTR, 5=promoter
            d,e = (None for i in range(2)) # for exons: which one, and how far
            d1,d2,e1,e2 = (None for i in range(4)) # for introns: between which exons and how far?
            f = None; #for introns: how many bases in the partially completed codon?
    
            #print df.iloc[tindex]['#name']
            #in promoter region ?
            if position < df.iloc[tindex]['txStart']:
                c = 5
                d = df.iloc[tindex]['txStart'] - position
                apList = [c,d,e,d1,d2,e1,e2,f]
                coordData.loc[tindex,['c', 'd', 'e', 'd1', 'd2', 'e1', 'e2','f']] = apList
                #print "In Promoter 1",position,df.iloc[tindex]['txStart'],"c=",c,"d=",d,"e=",e,"e1=",e1,"e2=",e2,"d1=",d1,"d2=",d2,"f=",f
                continue
            else:
                c = None
                d = None
            if position > df.iloc[tindex]['txEnd']:
                c = 5 
                d = position - df.iloc[tindex]['txEnd']
                apList = [c,d,e,d1,d2,e1,e2,f]
                #coordData.add(apList)
                coordData.loc[tindex,['c', 'd', 'e', 'd1', 'd2', 'e1', 'e2','f']] = apList
                #print "In Promoter 2","c=",c,"d=",d,"e=",e,"e1=",e1,"e2=",e2,"d1=",d1,"d2=",d2,"f=",f
                continue
            else:
                c = None
                d = None 
            
            #in UTR region ?
            if(df.iloc[tindex]['strand'] == '+'):
                if df.iloc[tindex]['cdsStart'] > position:
                    c = 4 
                    d = df.iloc[tindex]['cdsStart'] - position
                    apList = [c,d,e,d1,d2,e1,e2,f]
                    coordData.loc[tindex,['c', 'd', 'e', 'd1', 'd2', 'e1', 'e2','f']] = apList
                    #print "In 5'UTR","c=",c,"d=",d,"e=",e,"e1=",e1,"e2=",e2,"d1=",d1,"d2=",d2,"f=",f
                    continue
                else:
                    c = None
                    d = None
                if position > df.iloc[tindex]['cdsEnd']:
                    c = 3 
                    d = position - df.iloc[tindex]['cdsStart']
                    apList = [c,d,e,d1,d2,e1,e2,f]
                    coordData.loc[tindex,['c', 'd', 'e', 'd1', 'd2', 'e1', 'e2','f']] = apList
                    #print "In 3'UTR","c=",c,"d=",d,"e=",e,"e1=",e1,"e2=",e2,"d1=",d1,"d2=",d2,"f=",f
                    continue
                else:
                    c = None
                    d = None
            else:
                if df.iloc[tindex]['cdsStart'] > position:
                    c = 3
                    d = (df.iloc[tindex]['cdsStart'] - position)
                    apList = [c,d,e,d1,d2,e1,e2,f]
                    coordData.loc[tindex,['c', 'd', 'e', 'd1', 'd2', 'e1', 'e2','f']] = apList
                    #print "In 3'UTR","c=",c,"d=",d,"e=",e,"e1=",e1,"e2=",e2,"d1=",d1,"d2=",d2,"f=",f
                    continue
                else:
                    c = None
                    d = None
                if position > df.iloc[tindex]['cdsEnd']:
                    c = 4
                    d = position - df.iloc[tindex]['cdsStart']
                    apList = [c,d,e,d1,d2,e1,e2,f]
                    coordData.loc[tindex,['c', 'd', 'e', 'd1', 'd2', 'e1', 'e2','f']] = apList
                    #print "In 5'UTR","c=",c,"d=",d,"e=",e,"e1=",e1,"e2=",e2,"d1=",d1,"d2=",d2,"f=",f
                    continue
                else:
                    c = None
                    d = None
            #In exonic region
            exonStarts = filter(None,df.iloc[tindex]['exonStarts'].split(","))
            exonEnds = filter(None,df.iloc[tindex]['exonEnds'].split(","))
            in_exon = None
            for k in range(len(exonStarts)):
                if(int(exonStarts[k])<= int(position) and int(exonEnds[k]) >= int(position)):
                    in_exon = k+1
                    #print in_exon
                    break
            if(in_exon):
                c = 1
                e = in_exon 
                apList = [c,d,e,d1,d2,e1,e2,f]
                coordData.loc[tindex,['c', 'd', 'e', 'd1', 'd2', 'e1', 'e2','f']] = apList
                #print "In Exon","c=",c,"d=",d,"e=",e,"e1=",e1,"e2=",e2,"d1=",d1,"d2=",d2,"f=",f
                continue
            else:
                c = None
                e = None
            #In Intronic Region
            c = 2
            exonCount = df.iloc[tindex]['exonCount']
            exonFrames = filter(None,df.iloc[tindex]['exonFrames'].split(","))
            for k in range(exonCount):
                if(int(exonEnds[k]) < int(position) and int(exonStarts[k+1]) > position ):
                    if(df.iloc[tindex]['strand'] == '+'):
                        f = exonFrames[k+1]
                    else:
                        f = exonFrames[k]
                    e1 = k+1
                    e2 = k+2
                    d1 = int(position) - int(exonEnds[k])
                    d2  = int(exonStarts[k+1]) - int(position)
                    d = (min(d1,d2))
                    apList = [c,d,e,d1,d2,e1,e2,f]
                    coordData.loc[tindex,['c', 'd', 'e', 'd1', 'd2', 'e1', 'e2','f']] = apList
                    break
            #print "In Intron","c=",c,"d=",d,"e=",e,"e1=",e1,"e2=",e2,"d1=",d1,"d2=",d2,"f=",f
        if(autoSelect):
            (geneName,transcript,desc,zone,strandDirection,intronnum,intronframe) = ft.FindATranscript(coordData, df)
        else:
            (geneName,transcript,desc,zone,strandDirection,intronnum,intronframe) = ft.FindAllTranscripts(coordData, df)
    #For Intergenic calls
    else:
        distBefore = abs(df.iloc[idxList]['txStart'] - position)
        distAfter = abs(df.iloc[idxList]['txEnd'] + position)
        for y in 1000.0 ** (np.arange(1,4,0.3)):
            cmpDB = distBefore[distBefore <= y]
            if(not cmpDB.empty):
                beforeIdx = cmpDB.idxmin(axis=1)
                geneName = df.iloc[beforeIdx]['name2']
                strandDirection = df.iloc[beforeIdx]['strand']
                transcript = df.iloc[beforeIdx]['#name']
                zone = 0
                desc = 'IGR: ' + hp.bp2str(distBefore[beforeIdx],2) + ' before ' + geneName + '(' + strandDirection + ')' 
                #print a
                break
            cmpDA = distAfter[distAfter <= y]
            if(not cmpDA.empty):
                afterIdx = cmpDA.idxmin(axis=1)
                geneName = df.iloc[afterIdx]['name2']
                strandDirection = df.iloc[afterIdx]['strand']
                transcript = df.iloc[afterIdx]['#name']
                zone = 0
                desc = 'IGR: ' + hp.bp2str(distBefore[afterIdx],2) + ' after ' + geneName + '(' + strandDirection + ')' 
                #print a
                break
    return(geneName,transcript,desc,zone,strandDirection,intronnum,intronframe)   