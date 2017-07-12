"""
Created on 25/11/2014.

@author: Ronak H Shah

"""
from __future__ import division
import pandas as pd

'''
Read the Human Annotation using pandas
'''


def ReadFile(infile):
    dataDF = pd.read_csv(infile, sep='\t', header=0, keep_default_na='True')
    return(dataDF)

'''
Read the Human Canonical Transcript Files and return a dictionary
'''


def ReadTranscriptFile(infile):
    dataDict = None
    with open(infile) as fin:
        rows = (line.strip("\n").split('\t') for line in fin)
        dataDict = {row[0]: row[1:] for row in rows}
    return(dataDict)

'''
Using the txStart and txEnd extend the Promoter region
and assign its value to geneStart and geneEnd
'''


def ExtendPromoterRegion(df, distance):
    if(distance):
        distance = distance
    else:
        distance = 3000
    txStart = df['txStart']
    txEnd = df['txEnd']
    p_mask = (df['strand'] != "-")
    p_valid = df[p_mask]
    df['geneStart'] = txStart
    df.loc[p_mask,'geneStart'] = p_valid['txStart'] - int(distance)
    n_mask = (df['strand'] != "+")
    n_valid = df[n_mask]
    df['geneEnd'] = txEnd
    df.loc[n_mask,'geneEnd'] = n_valid['txStart'] + int(distance)
#     if(str(df['strand']) == '+'):
#         geneStart = txStart - int(distance)
#     else:
#         geneStart = txStart
#     if(str(df['strand']) == '-'):
#         geneEnd = txEnd + int(distance)
#     else:
#         geneEnd = txEnd
#     df['geneStart'] = geneStart
#     df['geneEnd'] = geneEnd
    return(df)

'''
Convert the total number of bases to a string value like kb,mb
'''


def bp2str(b, decimal_places):
    if not 'decimal_places' in globals():
        decimal_places = 0
    fs = '{0:.' + str(decimal_places) + 'f}'
    x = []
    x.append(fs.format(b) + 'bp')
    if b >= 1000:
        x.append(fs.format(b / 1000) + 'Kb')
    if b >= 1000000:
        x.append(fs.format(b / 1000000) + 'Mb')
    s = min(x, key=len)
    return s
