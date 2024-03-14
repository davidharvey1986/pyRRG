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
#twine upload dist/dist.tar.gz

version='0.4.1'
         
    
INCDIRS=['.']

packages = ['pyRRG', 'RRGtools','asciidata','stilts']
package_dir = {'RRGtools':'./lib/RRGtools',
                   'pyRRG':'./src',
               'asciidata':'./lib/asciidata',
               'stilts':'./lib/stilts/'}
package_data = \
    {'pyRRG': ['psf_lib_jwst/*/*','psf_lib/*/*','sex_files/*','*.pkl'],\
     'stilts':['*.jar']}





setup   (       name            = "pyRRG",
                version         = version,
                author          = "David Harvey",
                author_email    = "david.harvey@epfl.ch",
                description     = "pyRRG module",
                license         = 'MIT',
                packages        = packages,
                package_dir     = package_dir,
                package_data    = package_data,
                scripts         = ['scripts/pyRRG','scripts/stilts.sh'],
                url = 'https://github.com/davidharvey1986/pyRRG', # use the URL to the github repo
                download_url = 'https://github.com/davidharvey1986/pyRRG/archive/'+version+'.tar.gz',
                install_requires=['scikit-learn',\
                                   'numpy', 'tqdm', \
                                    'scipy'],                          
        )


