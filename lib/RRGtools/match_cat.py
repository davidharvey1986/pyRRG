import os as os
import pyfits as py
import numpy as np

import color_color as color
from numpy.lib.recfunctions import append_fields as append_rec

def run_match( cat_A, cat_B, \
                   search_rad=2 ):
    
    command_str ='stilts.sh tmatch2 in1="'+cat_A+'" in2="'+\
        cat_B+'" matcher=sky values1="RA DEC" values2="RA DEC" params="'\
        +str(search_rad)+'" out=matched_A_B.fits'


    os.system(command_str)

    matched_cat = py.open('matched_A_B.fits')

    return matched_cat

def match_cat( cat_A, cat_B, cleanup=True, \
               search_rad=2.):
    '''
    Run the stilts command tmatch2 to match two catalogues
    since this is hardcorded for VST it will have to combined
    all 32 chips into one chip and then match.
    '''
    catalogue_A = py.open( cat_A )
    catalogue_B = py.open( cat_B )

    exptime_A = catalogue_A[0].header['EXPTIME']
    exptime_B = catalogue_B[0].header['EXPTIME']
    
    for i in range(1, len(catalogue_A)):


        #Since I am concatenating cats I need keywords
        #to calc mag, so will calc mag here and append to the
        #numpy rec array
        fluxA = catalogue_A[i].data['Aper_flux_3']
        fluxB = catalogue_B[i].data['Aper_flux_3']
        
        apcor_A = catalogue_A[i].header['APCOR3']
        apcor_B = catalogue_B[i].header['APCOR3']
        
        zpt_A = catalogue_A[i].header['NIGHTZPT']
        zpt_B = catalogue_B[i].header['NIGHTZPT']
        
        mag_A = color.magnitude( fluxA, zpt_A, exptime_A, apcor_A )
        mag_B = color.magnitude( fluxB, zpt_B, exptime_B, apcor_B )

        
        catalogue_A[i].data = append_rec(catalogue_A[i].data,\
                                            'MAG', mag_A)
        catalogue_B[i].data = append_rec(catalogue_B[i].data, \
                                            'MAG', mag_B)
    
    
        if i == 1:
            total_cat_A = catalogue_A[i].data
            total_cat_B = catalogue_B[i].data
        else:
            total_cat_A = np.append( total_cat_A, catalogue_A[i].data )
            total_cat_B = np.append( total_cat_B, catalogue_B[i].data )


        
    HDU_cat_A = py.BinTableHDU( total_cat_A )
    HDU_cat_B = py.BinTableHDU( total_cat_B )

    HDU_list_A = py.HDUList( [ catalogue_A[0], HDU_cat_A ] )
    HDU_list_B = py.HDUList( [ catalogue_B[0], HDU_cat_B ] )

    
    HDU_list_A.writeto('total_cat_A.fits', clobber=True)
    HDU_list_B.writeto('total_cat_B.fits', clobber=True)

    run_match( 'total_cat_A.fits', 'total_cat_B.fits', \
                   search_rad=search_rad  )
    if cleanup:
        cleanup_str = 'rm -fr total_cat_A.fits total_cat_A.fits matched_A_B.fits'
        os.system(cleanup_str)
    return matched_cat
