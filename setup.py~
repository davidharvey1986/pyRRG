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

class checkModules(install):

    def run(self):
        install.run(self)
        try:
            stilts_path = subprocess.check_output(['which','stilts.sh'])
        except:
            raise ValueError('Cannot find STILTS please install and ensure it is in the shell path')

        try:
            import pyfits as pyfits
            if pyfits.__version__ < 3.3:
                raise ImportError('Code only tested on Pyfits 3.3 or later')
        except:
            raise ImportError('Cant find PyFITS')

        try:
            stilts_path = subprocess.check_output(['which','sex'])
        except:
            raise ImportError('Cannot find SExtractir please install and ensure it can be called with "sex"')

        try:
            import pickle as pkl
        except:
            raise ImportError('Cannot find pickle, plesae run easy_install pickle')

        try:
            import idlsave as idl
        except:
            raise ImportError('Cannot find idlsave, plesae run easy_install idlsave')
            
    
INCDIRS=['.']

packages = ['pyRRG', 'RRGtools']
package_dir = {'RRGtools':'./lib/RRGtools',
                   'pyRRG':'./src'}
package_data = {'pyRRG': ['psf_lib/*/*',
                              'sex_files/*']}



# in the setup function:
cmdclass={'install': checkModules}


setup   (       name            = "pyRRG",
                version         = "0.0.2",
                author          = "David Harvey",
                author_email    = "david.harvey@epfl.ch",
                description     = "pyRRG module",
                cmdclass        = cmdclass,
                packages        = packages,
                package_dir     = package_dir,
                package_data    = package_data,
                url = 'https://github.com/davidharvey1986/pyRRG', # use the URL to the github repo
                download_url = 'https://github.com/davidharvey1986/pyRRG/archive/0.0.2.tar.gz',
                          
        )


