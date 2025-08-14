iAnnotateSV: Annotation of structural variants detected from NGS
================================================================

:Author: `Ronak H Shah <http://github.com/rhshah>`_
:Contributors: `Gowtham Jayakumaran <https://github.com/andurill>`_ and `Ian Johonson <https://github.com/ionox0>`_ and `Kofi Amoah <https://github.com/kofiamoah>`_
:Contact: rons.shah@gmail.com
:Source code: http://github.com/rhshah/iAnnotateSV
:License: `Apache License 2.0 <http://www.apache.org/licenses/LICENSE-2.0>`_

.. image:: https://img.shields.io/pypi/v/iAnnotateSV.svg
        :target: https://pypi.python.org/pypi/iAnnotateSV

.. image:: https://zenodo.org/badge/18929/rhshah/iAnnotateSV.svg
   :target: https://zenodo.org/badge/latestdoi/18929/rhshah/iAnnotateSV

iAnnotateSV is a Python library and command-line software toolkit to annotate and
visualize structural variants detected from Next Generation DNA sequencing data. This works for majority is just re-writing of a tool called dRanger_annotate written in matlab by Mike Lawrence at Broad Institue. 
But it also has some additional functionality and control over the annotation w.r.t the what transcripts to be used for annotation.
It is designed for use with hybrid capture, including both whole-exome and custom target panels, and
short-read sequencing platforms such as Illumina.

Updates
========

In this update, some updates were made so that the iAnnotateSV package has a compatible sv table for hg38 Ensembl annotations. I have also trimmed down the package so that it does not require the external annotations such as DGv, COSMIC, kinase domain and repeats. So this update only produces the site descriptions column for the positions given in the input


Quick Usage
===========

::

    usage: iAnnotateSV.py [options]

    Annotate SV based on a specific human reference

    optional arguments:
    -h, --help            show this help message and exit
    -v, --verbose         make lots of noise [default]
    -r hg19, --refFileVersion hg19
                            Which human reference file to be used, hg18,hg19 or
                            hg38
    -rf hg19.sv.table.txt, --refFile hg19.sv.table.txt
                            Human reference file location to be used
    -ofp test, --outputFilePrefix test
                            Prefix for the output file
    -o /somedir, --outputDir /somedir
                            Full Path to the output dir
    -i svfile.txt, --svFile svfile.txt
                            Location of the structural variants file to annotate
    -d 3000, --distance 3000
                            Distance used to extend the promoter region
    -a, --autoSelect      Auto Select which transcript to be used[default]
    -c canonicalExons.txt, --canonicalTranscripts canonicalExons.txt
                            Location of canonical transcript list for each gene.
                            This file is required now
    -p, --plotSV          Plot the structural variant in question (very primitive)
    -u uniprot.txt, --uniprotFile uniprot.txt
                            Location of UniProt list contain information for
                            protein domains. Use only if you want to plot the
                            structural variant


Input file format is a tab-delimited file containing:

chr1  pos1  str1  chr2  pos2  str2

as the header and where:

* **chr1:** Its the chromosome name for first break point [1,2,3,4,5,6,7 etc..],
* **pos1:** Its the chromosome loaction for first break point [1-based],
* **str1:** Its the read direction for the first break point [0=top/plus/reference, 1=bottom/minus/complement],
* **chr2:** Its the chromosome name for second break point [1,2,3,4,5,6,7 etc..],
* **pos2:** Its the chromosome loaction for second break point [1-based],
* **str2:** Its the read direction for the second break point [0=top/plus/reference, 1=bottom/minus/complement], 

Output file will is a tab-delimited file containing:

chr1  pos1  str1  chr2  pos2  str2  gene1 transcript1 site1 gene2 transcript2 site2 fusion

as the header and where:

* **chr1** : Its the chromosome name for first break point [1,2,3,4,5,6,7 etc..],
* **pos1** : Its the chromosome loaction for first break point [1-based],
* **str1** : Its the read direction for the first break point [0=top/plus/reference, 1=bottom/minus/complement],
* **chr2** : Its the chromosome name for second break point [1,2,3,4,5,6,7 etc..],
* **pos2** : Its the chromosome loaction for second break point [1-based],
* **str2** : Its the read direction for the second break point [0=top/plus/reference, 1=bottom/minus/complement],
* **gene1** : Gene for the first break point,
* **transcript1** : Transcript used for the first breakpoint,
* **site1** : Explanation of the site where the first breakpoint occurs [Example=>Intron of EWSR1(+):126bp after exon 10],
* **gene2** : Gene for the second break point,
* **transcript2** : Transcript used for the second breakpoint,
* **site2** : Explanation of the site where the second breakpoint occurs [Example=>Intron of ERG(-):393bp after exon 4],
* **fusion** : Explanation if the evnet leads to fusion or not. [Example=>Protein Fusion: in frame  {EWSR1:ERG}]


