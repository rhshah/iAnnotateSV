"""
Created on 25/11/2014.

@author: Ronak H Shah

"""
from __future__ import division
import polars as pl
import logging
from rich import print
from rich.logging import RichHandler

FORMAT = "%(message)s"
logging.basicConfig(
    level="INFO", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)

log = logging.getLogger("rich")

'''
Read the Human Annotation using polars
'''


def ReadFile(infile):
    """
    Reads a tab-separated file into a Polars DataFrame.

    Args:
        infile (str): Path to the input file.

    Returns:
        pl.DataFrame: A Polars DataFrame containing the data, or None if an error occurs.
    """
    log.info(f"Reading file: {infile}")
    try:
        df = pl.read_csv(infile, separator='\t', infer_schema_length=1000)
        log.info(f"Successfully read file: {infile}")
        # Remove '#' from column names
        new_columns = [col.replace('#', '') for col in df.columns]
        df = df.rename(dict(zip(df.columns, new_columns)))
        return df
    except Exception as e:
        log.error(f"Error reading file {infile}: {e}")
        return None


'''
Read the Human Canonical Transcript Files and return a dictionary
'''


def ReadTranscriptFile(infile):
    """
    Reads a transcript file and returns a dictionary.

    Args:
        infile (str): Path to the input file.

    Returns:
        dict: A dictionary where keys are transcript IDs and values are lists of associated data.
    """
    log.info(f"Reading transcript file: {infile}")
    dataDict = {}
    try:
        with open(infile) as fin:
            for line in fin:
                row = line.strip("\n").split('\t')
                if len(row) > 0:
                    dataDict[row[0]] = row[1:]
        log.info(f"Successfully read transcript file: {infile}")
        return dataDict
    except Exception as e:
        log.error(f"Error reading transcript file {infile}: {e}")
        return None

'''
Using the txStart and txEnd extend the Promoter region
and assign its value to geneStart and geneEnd
'''


def ExtendPromoterRegion(df, distance):
    """
    Extends the promoter region of each gene by a specified distance.

    Args:
        df (pl.DataFrame): Input DataFrame containing gene annotations.
        distance (int): Distance to extend the promoter region.

    Returns:
        pl.DataFrame: A new DataFrame with extended promoter regions.
    """
    distance = int(distance) if distance else 3000
    log.info(f"Extending promoter region by {distance} bp")
    try:
        df = df.with_columns([
            pl.col("txStart").alias("geneStart"),
            pl.col("txEnd").alias("geneEnd")
        ])
        log.info("Successfully extended promoter region")
        return(df)
    except Exception as e:
        log.error(f"Error extending promoter region: {e}")
        return None

'''
Convert the total number of bases to a string value like kb,mb
'''


def bp2str(b, decimal_places):
    """
    Converts a number of base pairs to a human-readable string (e.g., "1.23Kb").

    Args:
        b (int): Number of base pairs.
        decimal_places (int): Number of decimal places to use.

    Returns:
        str: A human-readable string representing the number of base pairs.
    """
    if 'decimal_places' not in locals():
        decimal_places = 0
    fs = '{0:.' + str(decimal_places) + 'f}'
    x = [f'{fs.format(b)}bp']
    if b >= 1000:
        x.append(f'{fs.format(b / 1000)}Kb')
    if b >= 1000000:
        x.append(f'{fs.format(b / 1000000)}Mb')
    return min(x, key=len)