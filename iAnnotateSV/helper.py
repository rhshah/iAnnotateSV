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
        distance = int(distance)
    else:
        distance = 3000
    df['geneStart'] = df['txStart']
    df['geneEnd'] = df['txEnd']
    #df['geneStart'] = df.apply(lambda row: row['txStart'] if str(row['strand']) == '+' else row['txStart'] + distance, axis=1)
    #df['geneEnd'] = df.apply(lambda row: row['txEnd'] if str(row['strand']) == '-' else row['txEnd'] + distance, axis=1)
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
