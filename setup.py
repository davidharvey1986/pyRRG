#!/usr/bin/env python3
from setuptools import setup, find_packages

# Avoid importing numpy at the top-level of setup.py
long_description = """\
This module uses the RRG method to measure the shapes of galaxies
in Hubble Space Telescope data
"""

version = '0.4.2'

packages = ['pyRRG', 'RRGtools', 'asciidata', 'stilts']
package_dir = {
    'RRGtools': './lib/RRGtools',
    'pyRRG': './src',
    'asciidata': './lib/asciidata',
    'stilts': './lib/stilts'
}
package_data = {
    'pyRRG': [
        'psf_lib/tinytim/ACS/*/*',
        'psf_lib/webbpsf/NIRCAM/*/*',
        'sex_files/*'
    ],
    'stilts': ['*.jar']
}

setup(
    name="pyRRG",
    version=version,
    author="David Harvey",
    author_email="david.harvey@epfl.ch",
    description="pyRRG module",
    long_description=long_description,
    license='MIT',
    packages=packages,
    package_dir=package_dir,
    package_data=package_data,
    scripts=['scripts/pyRRG', 'scripts/stilts.sh'],
    url='https://github.com/davidharvey1986/pyRRG',
    download_url=f'https://github.com/davidharvey1986/pyRRG/archive/{version}.tar.gz',
    install_requires=[
        'numpy',
        'scikit-learn',
        'tqdm',
        'scipy',
        'astropy',
        'matplotlib'
    ],
    python_requires='>=3.8'
)
