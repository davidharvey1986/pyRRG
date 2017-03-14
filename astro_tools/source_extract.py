import pysex as sex
import pyfits as fits
import ipdb as pdb
import numpy as np
import os as os
import match_cat as mc
from numpy.lib.recfunctions import append_fields as append_rec


def source_extract( image_name, weight_file, zero_point=None,
                    dataDir='PWD', outfile=None, return_sources=True,
                    stilts_dir='.', conf_path='.'):
    '''
    Given that source extration is a difficult thing for the photometry
    I will use the Mathilde mehtod, in Python to do a wrapper
    for source_extractor.


    The process is the same as in mathilde_setract.pro
    THe weight file is optional however should be used.
    zero_point is the zero point of the iamge, if not used, it
    will guess it from the header of the file and assume it is ACS

    This will run two runs of sextractor, cold and hot, where cold only finds
    big things and hot finds small things

    Then it combines the two catalogues
    
    '''
    check_sex_files( conf_path)
    
    if dataDir == 'PWD':
        dataDir = os.getcwd()
        
    if zero_point is None:
        zero_point = acs_zero_point ( fits.open( image_name )[0].header )
        
    conf_args = {'WEIGHT_IMAGE': dataDir+'/'+weight_file,
                 'MAG_ZEROPOINT': zero_point,
                 'WEIGHT_TYPE':'MAP_WEIGHT',
                 'PARAMETERS_NAME':conf_path+'/rrg.param',
                 'STARNNW_NAME':conf_path+'/default.nnw',
                 'FILTER_NAME':conf_path+'/gauss_5.0_9x9.conv'}
        



    #First run cold run
    cold_conf = conf_path+'/HFF_cold.param'
    

        
    cold_sources = sex.run( image_name, \
                                conf_file=cold_conf, \
                                conf_args=conf_args )

    cold_names = np.array(cold_sources.columns.names)
    cold_names[ cold_names == 'X_WORLD' ] = 'RA'
    cold_names[ cold_names == 'Y_WORLD' ] = 'DEC'
    cold_sources.columns.names = tuple( cold_names )
    #Second hot 
    hot_conf = conf_path+'/HFF_hot.param'
        
    hot_sources = sex.run( image_name, \
                               conf_file=hot_conf, \
                               conf_args=conf_args )

    hot_names = np.array(hot_sources.columns.names)

    hot_names[ hot_names == 'X_WORLD' ] = 'RA'
    hot_names[ hot_names == 'Y_WORLD' ] = 'DEC'
    hot_sources.columns.names = tuple( hot_names )

    #The NYMBER is a weird thing
    hot_sources['NUMBER'] = np.arange( len(hot_sources['NUMBER'])) +1
    cold_sources['NUMBER'] = np.arange( len(cold_sources['NUMBER'])) +1
    
    fits.writeto( 'cold_sources.fits', cold_sources, clobber=True )
    fits.writeto( 'hot_sources.fits', hot_sources, clobber=True )


    
    matched_sources= mc.run_match( 'cold_sources.fits',
                                'hot_sources.fits',
                                       stilts_path=stilts_dir)
    
    for iField in hot_names:
        hot_sources[iField][ matched_sources[1].data['NUMBER_2'] -1 ] = \
          cold_sources[iField][ matched_sources[1].data['NUMBER_1'] - 1]

  
    if outfile is not None:
        fits.writeto( outfile, hot_sources)
    if return_sources:
        return hot_sources
    
    
def acs_zero_point( header ):
    '''
    PURPOSE : TO GET THE ZERO POINT OF THE ACS IMAGE USING
    THE EQUATIONS FROM THE INTERNET

    THIS IS THE MAGNITUDE TO CONVERT FROM CPS TO
    AB MAGNITUDES

    SO AB_MAG = -2.5 LOG(E^-) + ZP

    INPUTS : THE HEADER FILE OF THE ACS IMAGE

    OUTPUT : A SCALAR OF THE ZERO-POINT MAG IN AB
    '''
    return -2.5*np.log10(header['PHOTFLAM']) \
        +header['PHOTZPT'] - 5.0*np.log10(header['PHOTPLAM'])+18.6921

def check_sex_files( sex_dir ):
    '''
    Check all the sex files to see if they exist
    '''
    
    if (not os.path.isfile(sex_dir+'/HFF_hot.param')):
        raise ValueError('HFF_hot.param not found')
    
    if (not os.path.isfile(sex_dir+'/HFF_cold.param')):
        raise ValueError('HFF_cold.param not found')
                
    if (not os.path.isfile(sex_dir+'/rrg.param')):
        raise ValueError('rrg.param not found')
