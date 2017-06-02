import sys,os,string,glob,subprocess

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

class checkStilts(install):

    def run(self):
        install.run(self)
        try:
            stilts_path = subprocess.check_output(['which','stilts.sh'])
        except:
            raise ValueError('Cannot find STILTS please install and ensure it is in the shell path')
    
INCDIRS=['.']

packages = ['pyRRG', 'RRGtools']
package_dir = {'RRGtools':'./lib/RRGtools',
                   'pyRRG':'./src'}
package_data = {'pyRRG': ['psf_lib/*/*',
                              'sex_files/*']}



# in the setup function:
cmdclass={'install': checkStilts}


setup   (       name            = "pyRRG",
                version         = "0.0.1",
                author          = "David Harvey",
                author_email    = "david.harvey@epfl.ch",
                description     = "pyRRG module",
                platform        = 'MAC OS X',
                cmdclass        = cmdclass,
                packages        = packages,
                package_dir     = package_dir,
                package_data    = package_data,
                          
        )


