import numpy as np
from astropy.io import fits as py
from astropy import wcs
import os as os

def deg2pix( fits, ra, dec, extension=None):
    '''
    Given a fits file, convert from ra and dec in degrees to
    x and y pix

    '''
    wcs_obj = wcs.WCS( py.open(fits)[extension].header )
    if isinstance(ra, float):
        radec = [ra, dec]
        return wcs_obj.all_world2pix(np.array(radec)[np.newaxis,:], 1, quiet=True)[0]
    else:
        return wcs_obj.all_world2pix(ra, dec, 1, quiet=True)

def pix2deg( fits_image, x_image, y_image,  extension=0):
    '''
    Given a fits file, convert from ra and dec in degrees to
    x and y pix

    '''
    wcs_obj = wcs.WCS( py.open(fits_image)[extension].header)
  
    return wcs_obj.all_pix2world( x_image, y_image, 1)
    
def deg2pix_flt( fits, ra, dec, postage_stamp=0, cut=False):
    '''
    Run the deg2pix function but for 2 separate chips
    concatenate them and remove any objects outside the chip

    ostage_stamp allows boardering of
    '''
    
    #chip1
    x_chip1, y_chip1 =  deg2pix( fits, ra, dec, ext=['sci',1])
    
    #chip2
    x_chip2, y_chip2 =  deg2pix( fits, ra, dec, ext=['sci',2])
    
    
    #Check that the objects lie within the chip and
    #which chip they are in
    inchip1 = (y_chip1 < 2048-postage_stamp) & (y_chip1 > postage_stamp) & \
        (x_chip1 > postage_stamp) & (x_chip1 < 4096-postage_stamp)
    inchip2 = (y_chip2 < 2048-postage_stamp) & (y_chip2 > postage_stamp) & \
        (x_chip2 > postage_stamp) & (x_chip2 < 4096-postage_stamp)
            
    x = np.append( x_chip1[ inchip1 ], x_chip2[ inchip2 ])
    y = np.append( y_chip1[ inchip1 ], y_chip2[ inchip2 ]+2048)
    print(x, y)
    return x, y

 



def pix2deg_flt( fits, x, y):
    '''
    Run the deg2pix function but for 2 separate chips
    concatenate them and remove any objects outside the chip
    '''

    #chip1
    ra_chip1, dec_chip1 =  pix2deg( fits+'[sci,1]', \
                                    x[ y < 2048 ], \
                                    y[ y < 2048 ])

    #chip2
    ra_chip2, dec_chip2 =  pix2deg( fits+'[sci,2]', \
                                    x[ y > 2048], \
                                    y[ y > 2048]-2048., \
                                    coordfile="xy2sky.par")


    #Check that the objects lie within the chip and
    #which chip they are in
            
    ra = np.append(ra_chip1, ra_chip2)
    dec = np.append( dec_chip1, dec_chip2)

    return ra, dec

def hmstodd( ra, dec ):
    '''
    Convert from hms to dd
    INPUTS : RA AND DEC ARE STRINGS IN THE FORMAT
    H:M:S SEPARATED BY COLONS
    '''
    ra_float = np.array(ra.split(':')).astype(float)
    dec_float = np.array(dec.split(':')).astype(float)

    
    ra_deg = (ra_float[0] + ra_float[1]/60. + ra_float[2]/3600.)/24.*360.
    if '-' in dec.split(':')[0]:
        hem = -1.
    else:
        hem = 1.
    
    dec_deg = dec_float[0]+hem*dec_float[1]/60. + hem*dec_float[2]/3600.
     
    return ra_deg, dec_deg
    
