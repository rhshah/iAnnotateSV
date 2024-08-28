from __future__ import print_function

try:
    from setuptools import setup, find_packages, Extension
except ImportError:
    from distutils.core import setup, find_packages
import io
import codecs
import os
import sys

import iAnnotateSV

here = os.path.abspath(os.path.dirname(__file__))


def read(*filenames, **kwargs):
    encoding = kwargs.get("encoding", "utf-8")
    sep = kwargs.get("sep", "\n")
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)


# long_description = read('README.rst', 'CHANGES.rst')

setup(
    name="iAnnotateSV",
    version=iAnnotateSV.__version__,
    description="The module helps to annotate structural variants called using NGS on human.",
    long_description="See Readme",
    include_package_data=True,
    url="https://github.com/rhshah/iAnnotateSV",
    download_url="https://github.com/rhshah/iAnnotateSV/tarball/1.2.0",
    author=iAnnotateSV.__author__,
    author_email="rons.shah@gmail.com",
    license=iAnnotateSV.__license__,
    platforms="any",
    packages=["iAnnotateSV"],
    install_requires=[
        "numpy",
        "openpyxl",
        "pandas",
        "nose",
        "codecov",
        "coverage",
        "pillow",
        "biopython",
        "reportlab",
        "coloredlogs",
        "bump2version",
    ],
    zip_safe=False,
    classifiers=(
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Development Status :: 5 - Production/Stable",
    ),
)
