# iAnnotateSV: Annotation of structural variants detected from NGS

Author:   [Ronak H Shah](http://github.com/rhshah)

Contributors :   [Gowtham Jayakumaran](https://github.com/andurill) and [Ian Johonson](https://github.com/ionox0)

Contact:   <rons.shah@gmail.com>

Source code:   <http://github.com/rhshah/iAnnotateSV>

License:   [Apache License 2.0](http://www.apache.org/licenses/LICENSE-2.0)

[![image](https://img.shields.io/pypi/v/iAnnotateSV.svg)](https://pypi.python.org/pypi/iAnnotateSV)
[![image](https://zenodo.org/badge/18929/rhshah/iAnnotateSV.svg)](https://zenodo.org/badge/latestdoi/18929/rhshah/iAnnotateSV)

iAnnotateSV is a Python library and command-line software toolkit to
annotate and visualize structural variants detected from Next Generation
DNA sequencing data. This works for majority is just re-writing of a
tool called dRanger_annotate written in matlab by Mike Lawrence at Broad
Institue. But it also has some additional functionality and control over
the annotation w.r.t the what transcripts to be used for annotation. It
is designed for use with hybrid capture, including both whole-exome and
custom target panels, and short-read sequencing platforms such as
Illumina.

## Citation

We are in the process of publishing a manuscript describing iAnnotateSV
as part of the Structural Variant Detection framework. If you use this
software in a publication, for now, please cite our website
[iAnnotateSV](http://github.com/rhshah/iAnnotateSV).

## Acknowledgements

I would like to thanks Mike Lawrence from Broad Institute for sharing his code and Michael Berger for his insights into the dRanger_annotate tool.

## Required Packages

We require that you install:

python: [v3.10](https://www.python.org/downloads/release/)

pandas:  [v2.2.2](http://pandas.pydata.org/)

biopython:  [v1.84](http://biopython.org/wiki/Main_Page)

Pillow:  [v10.4.0](https://pypi.python.org/pypi/Pillow/)

openpyxl:  [v3.1.5](https://pypi.python.org/pypi/openpyxl/)

reportlab:  [v3.6.13](https://pypi.python.org/pypi/reportlab/)

coloredlogs:  [v15.0.1](https://pypi.python.org/pypi/coloredlogs)

## Utilities

[Generate Fusion Peptides](./generate_fusion_peptides/README.md) - Helen
Xie

## Quick Usage

If you know python I have created a small test script in /iAnnotateSV/test directory it runs a test on existing code and compares the result with the output file.

Else To Run:

- If you want to run with default options:

    ```bash
    python /path/to/iAnnotateSV.py -i svFile.txt -ofp outputfilePrefix -o /path/to/output/dir -r hg19 -d 3000
    ```

- If you want to run with your own transcripts:

    ```bash
    python path/to/path/to/iAnnotateSV.py -i svFile.txt -ofp outputfilePrefix -o /path/to/output/dir -r hg19 -d 3000 -c canonicalTranscripts.txt
    ```

- If you want to run with your own transcripts & make plots (making plots is a test module only):

    ```bash
    python path/to/iAnnotateSV.py -i svFile.txt -ofp outputfilePrefix -o /path/to/output/dir -r hg19 -d 3000 -c canonicalTranscripts.txt -u uniprot.txt -p
    ```
