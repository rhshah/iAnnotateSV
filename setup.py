from __future__ import print_function 
try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup, find_packages 
import io
import codecs
import os
import sys

import iAnnotateSV

here = os.path.abspath(os.path.dirname(__file__))

def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)

long_description = read('README.rst', 'CHANGES.rst')

setup(
      name='iAnnotateSV',
      version=iAnnotateSV.__version__,
      description='The module helps to annotate structural variants called using NGS on human.',
      long_description=long_description,
      include_package_data=True,
      url='https://github.com/rhshah/iAnnotateSV',
       download_url = 'https://github.com/rhshah/iAnnotateSV/tarball/1.0.9', 
      author=iAnnotateSV.__author__,
      author_email='rons.shah@gmail.com',
      license=iAnnotateSV.__license__,
      platforms='any',
      packages=['iAnnotateSV'],
      install_requires=[
          'numpy==1.14.0',
          'openpyxl==1.9.0',
          'pandas==0.16.2',
          'nose==1.3.7',
          'codecove==2.0.5',
          'coverage==4.3.4'
          'pillow==3.4.2',
          'biopython==1.65',
          'reportlab==3.3.0',
          'coloredlogs==5.2'
      ],
      zip_safe=False,
      classifiers=(
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'License :: OSI Approved :: Apache License',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Development Status :: 3'
        ),
      )