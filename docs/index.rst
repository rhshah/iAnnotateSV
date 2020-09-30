.. iAnnotateSV documentation master file, created by
   sphinx-quickstart on Fri Jan  2 14:46:47 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

iAnnotateSV: Annotation of structural variants detected from NGS
================================================================

:Author: Ronak H Shah
:Contact: rons.shah@gmail.com
:Source code: http://github.com/rhshah/iAnnotateSV
:License: `Apache License 2.0 <http://www.apache.org/licenses/LICENSE-2.0>`_

iAnnotateSV is a Python library and command-line software toolkit to annotate and
visualize structural variants detected from Next Generation DNA sequencing data. This works for majority is just a re-writing of a tool called dRanger_annotate written in matlab by Mike Lawrence at Broad Institue. 
But it also has some additional functionality and control over the annotation w.r.t the what transcripts to be used for annotation.
It is designed for use with hybrid capture, including both whole-exome and custom target panels, and
short-read sequencing platforms such as Illumina.

.. image:: https://img.shields.io/pypi/v/iAnnotateSV.svg
        :target: https://pypi.python.org/pypi/iAnnotateSV


.. image:: https://landscape.io/github/rhshah/iAnnotateSV/master/landscape.svg?style=flat
   :target: https://landscape.io/github/rhshah/iAnnotateSV/master
   :alt: Code Health


.. image:: https://zenodo.org/badge/18929/rhshah/iAnnotateSV.svg
   :target: https://zenodo.org/badge/latestdoi/18929/rhshah/iAnnotateSV


.. image:: https://travis-ci.org/rhshah/iAnnotateSV.svg?branch=master
    :target: https://travis-ci.org/rhshah/iAnnotateSV


.. image:: https://codecov.io/gh/rhshah/iAnnotateSV/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/rhshah/iAnnotateSV
   

Contents:

.. toctree::
   :maxdepth: 2
   
   modules

Citation
========

We are in the process of publishing a manuscript describing iAnnotateSV as part of the Structural Variant Detection framework.
If you use this software in a publication, for now, please cite our website `iAnnotateSV <http://github.com/rhshah/iAnnotateSV>`_.

Acknowledgements
================

I would like to thanks Mike Lawrence from Broad Institute for sharing his code and Michael Berger for his inshgts into the dRanger_Annoate tool.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

