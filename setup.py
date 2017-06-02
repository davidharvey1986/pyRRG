import sys,os,string,glob

from distutils.sysconfig import *
from distutils.core import setup,Extension
from distutils.command.build_ext import build_ext
from distutils.command.install import install
from distutils.command.install_data import install_data

import numpy

long_description = """\
This module uses the RRG method to measure the shapes of galaxies
in Hubble Space Telescope data
"""

      
INCDIRS=['.']

packages = ['pyRRG', 'RRGtools']
package_dir = {'RRGtools':'./lib/RRGtools',
                   'pyRRG':'./src'}
package_data = {'pyRRG': ['psf_lib/*/*.moms',
                              'sex_files/*']}


setup   (       name            = "pyRRG",
                version         = "0.0.1",
                author          = "David Harvey",
                author_email    = "david.harvey@epfl.ch",
                description     = "pyRRG module",
                
                packages        = packages,
                package_dir     = package_dir,
                package_data    = package_data,
                          
        )
