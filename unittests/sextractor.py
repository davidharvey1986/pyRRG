'''
This unit test to make sure sextractor works
'''

import pysex as pysex
import pyfits as fits
import subprocess 
import numpy as np
import os as os

def test_all():
    '''
    Test to ensure that sextractor works
    '''
    
    #1. Test it is in the path
    path_test()
    print 'Path test passed'
    #2. Test that it can find something and return a fits record
    image_test()
    print 'Image test passed'


def path_test():
    '''
    Test sextractor is in the path
    '''
    try:
        sex_dir = '/'.join(subprocess.check_output(['which','sex']).split('/')[:-1] )
    except:
        raise ImportError("Path Test Failed: Cant find SExtractor in path")
    

def image_test():
    '''
    Create an image
    Sextract it and return what should be a fits record
    '''
    create_test_image()
    source_extract_image()
    os.system('rm -fr test.fits')
def source_extract_image():
    '''
    Get the position from the image
    '''
    try:
        image =  pysex.run( 'test.fits')
    except:
        raise ValueError('Image Test Failed: '+\
                             'Pysex failed for an unknown reason')
    if not type(image) == fits.fitsrec.FITS_rec:
        raise ValueError('Image test Failed: '+\
                             'Pysex has not return a fits record.'+\
                             'Make sure user is using provided pysex.py')

    if len(image.data) == 0:
        raise ValueError('Image test Failed: '+\
                             'Pysex did not return or find the object')
        

def create_test_image():
    '''
    Create a test image of a Gaussian
    '''

    sigma_sqr = 20.**2

    imagex = np.linspace(0,200)-100.
    imagey = np.linspace(0,200)-100.
    xgrid, ygrid = np.meshgrid( imagex, imagey)
    
    rgrid_square = xgrid**2 + ygrid**2

    image = np.exp( -rgrid_square/(2.*sigma_sqr))

    fits.writeto( 'test.fits', image, clobber=True)
