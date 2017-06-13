import sys,os,string,glob,subprocess

from setuptools import setup,Extension
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install

import numpy

long_description = """\
This module uses the RRG method to measure the shapes of galaxies
in Hubble Space Telescope data
"""
#python setup.py register -r pypi
#sudo python setup.py sdist upload -r pypi

version='0.0.13'
         
    
INCDIRS=['.']

packages = ['pyRRG', 'RRGtools']
package_dir = {'RRGtools':'./lib/RRGtools',
                   'pyRRG':'./src'}
package_data = {'pyRRG': ['psf_lib/*/*',
                              'sex_files/*']}






setup   (       name            = "pyRRG",
                version         = version,
                author          = "David Harvey",
                author_email    = "david.harvey@epfl.ch",
                description     = "pyRRG module",
                license         = 'MIT',
                packages        = packages,
                package_dir     = package_dir,
                package_data    = package_data,
                url = 'https://github.com/davidharvey1986/pyRRG', # use the URL to the github repo
                download_url = 'https://github.com/davidharvey1986/pyRRG/archive/'+version+'.tar.gz',
                install_requires=['idlsave','pyfits>=3.3']
                          
        )


