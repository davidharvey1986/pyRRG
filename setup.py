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
#python setup.py sdist upload -r pypi

version='0.0.6'
class checkModules(install):

    def run(self):
        install.run(self)
        try:
            stilts_path = subprocess.check_output(['which','stilts.sh'])
        except:
            raise ValueError('Cannot find STILTS please install and ensure it is in the shell path')

    
        try:
            stilts_path = subprocess.check_output(['which','sex'])
        except:
            raise ImportError('Cannot find SExtractir please install and ensure it can be called with "sex"')

        try:
            import pickle as pkl
        except:
            raise ImportError('Cannot find pickle, plesae install')

            
    
INCDIRS=['.']

packages = ['pyRRG', 'RRGtools']
package_dir = {'RRGtools':'./lib/RRGtools',
                   'pyRRG':'./src'}
package_data = {'pyRRG': ['psf_lib/*/*',
                              'sex_files/*']}



# in the setup function:
cmdclass={'install': checkModules}


setup   (       name            = "pyRRG",
                version         = version,
                author          = "David Harvey",
                author_email    = "david.harvey@epfl.ch",
                description     = "pyRRG module",
                cmdclass        = cmdclass,
                packages        = packages,
                package_dir     = package_dir,
                package_data    = package_data,
                url = 'https://github.com/davidharvey1986/pyRRG', # use the URL to the github repo
                download_url = 'https://github.com/davidharvey1986/pyRRG/archive/'+version+'.tar.gz',
                install_requires=['idlsave','pyfits>=3.3']
                          
        )


