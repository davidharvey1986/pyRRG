#!/usr/local/bin/python3                                                        
import sys,os,string,glob,subprocess

from setuptools import setup,Extension
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install

import numpy

long_description = """\
This module uses the RRG method to measure the shapes of galaxies
in Hubble Space Telescope data
"""
#sudo python3 setup.py sdist upload -r pypi

version='0.1.1'
         
    
INCDIRS=['.']

packages = ['pyRRG', 'RRGtools','asciidata']
package_dir = {'RRGtools':'./lib/RRGtools',
                   'pyRRG':'./src',
               'asciidata':'./lib/asciidata'}
package_data = {'pyRRG': ['psf_lib/*/*','sex_files/*','*.pkl']}





setup   (       name            = "pyRRG",
                version         = version,
                author          = "David Harvey",
                author_email    = "david.harvey@epfl.ch",
                description     = "pyRRG module",
                license         = 'MIT',
                packages        = packages,
                package_dir     = package_dir,
                package_data    = package_data,
                scripts         = ['scripts/pyRRG'],
                url = 'https://github.com/davidharvey1986/pyRRG', # use the URL to the github repo
                download_url = 'https://github.com/davidharvey1986/pyRRG/archive/'+version+'.tar.gz',
                install_requires=['scikit-learn',\
                                   'pyfits', \
                                   'numpy', \
                                   'ipdb', 'pyraf',\
                                    'scipy'],                          
        )


