'''
This script will take the output from rrg and convert 
the shears in the format required to be used in lenstool
(using option 7, and the format a, b, theta)

'''

import pyfits as py
import os as os
#import mask_catalogue as mask_catalogue
import numpy as np

def rrg_to_lenstool( rrg_catalogue,
                         image_file,
                         output_catalogue=None,
                         lenstool_catalogue=None,
                         mask='mask.reg',
                         reference=None,
                         default_src_redshift=1.0):
    '''
    Take the input fits file that is the rrg_catalogue
    and calcualte the a, b and theta.

    Add these to the file and return the catalogue with added
    tagnames

    Also create a ascii file with the catalogue ready to go
    into lentool

    INPUTS :
    rrg_catalogue : the string of the rrg catalogue 
    KEYWORDS :
    reference : a two scalar vector with the ra and dec of the reference of the lenstool 
    output_catalogue : the name of the fits file which has a_lenstool, b_lenstool, theta_lenstool added to the fields
    lenstool_catalogue : a catalogue in the format x, y, a, b, theta, z, mag
    to be used in lenstool
    default_src_redshift : a scalar  of redshifts that the sources galaxies are at
    '''
    image = py.open(image_file)
    rrgCat = py.open( rrg_catalogue )
    
    if os.path.isfile( mask ):
        MaskedRRGCat = mask_catalogue.mask( rrgCat, mask )
    else:
        print 'WARNING: NO MASK FILE FOUND'
        MaskedRRGCat = rrgCat[1].data

    nGalaxies = len( MaskedRRGCat )

    theta = np.arctan2( MaskedRRGCat.gamma2, MaskedRRGCat.gamma1 )*180./np.pi/2.
    rotang = image[0].header['ORIENTAT']
    theta += rotang
    print image[0].header['ORIENTAT']
    gamma = np.sqrt(  MaskedRRGCat.gamma2**2 + MaskedRRGCat.gamma1**2)
    
    size = np.sqrt(( MaskedRRGCat.xx + MaskedRRGCat.yy)/2.)
    semi_major = size*np.sqrt(1.+gamma)
    semi_minor = size*np.sqrt(1.-gamma)

    if reference is None:
        reference = [ np.mean( MaskedRRGCat.RA ), 
                          np.mean( MaskedRRGCat.DEC ) ]


    if lenstool_catalogue is None:
        lenstool_catalogue = rrg_catalogue.split('.')[0]+'.lenstool'
        
    lenstoolCat = open( lenstool_catalogue, 'wb' )

    ngal = np.arange( nGalaxies )+1
    mag = MaskedRRGCat.MAG_AUTO
   

    redshifts = np.zeros( nGalaxies, float)+default_src_redshift
        
    write_array = np.transpose([ ngal, MaskedRRGCat.ra, MaskedRRGCat.dec, \
                                  semi_major, semi_minor, \
                                  theta, \
                                  redshifts,
                                     mag])
    #Remove nans from the file
    rows = np.arange(write_array.shape[0])[ (np.isnan(np.sum(write_array,axis=1)))]

    write_array = np.delete( write_array, rows, axis=0)
            
        
    
    np.savetxt( lenstool_catalogue, write_array,\
			fmt='%i %f %f %f %f %f %f %f', \
                    header='REFERENCE 3 %0.7f %0.7f' % (reference[0], reference[1]))
