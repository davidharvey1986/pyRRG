import numpy as np

import pyfits as py

def calc_shear( corrected_moments, galaxies, outfile,
                min_rad=6, mult=2,
                size_cut=[0., 100.],
                mag_cut=[22.5, 30],
                signal_noise_cut=3.,
                rhodes_factor=0.86,
                dataDir='./',
                expThresh=4,
                stat_type='median'):
    '''
    PURPOSE : TO WORK OUT VARIOUS FACTORS AND ADJUST THE PSF
    CORRECTED ELLIPICITY FOR THESE AND GET AN ESTAIMTE OF THE SHEAR


    INPUT  :
        CORRECTED_MOMENTS : A fitsrec structure of the corrected
                            and uncorrected moments to measure
                            and estimate the shear.
        GALAXIES : A integer vector of indexes that correspond to those
                   CORRECTED_MOMENTS that have been classified as a galaxy
        OUTFILE : The name of the fits file in which the shears, gamma1 and gamma2
                have been estimated as


    OPTIONAL INPUTS:
         MIN_RAD: The minimum measured radios that is used to weight each galaxies in the
                  measurement.
        SIZE_CUT [LO,HI] : scalar, the cut in galaxy size that is allowed to measure shear on
        MAG_CUT [LO,HI] : two scalar array, the cut in the galaxy magnitude that is allowed to measure shear on
        SIGNAL_NOSIE_CUT : scalar, the lowest acceptable signal to nosie allowed in the catalogue
        rhodes_factor : The RRG rhodes factor (see paper)
        dataDir : the directory in which all the data exists
        expThresh : the minimum number of exposures a galaxy has to have an acceptably robust shape.
        stat_type :  the stat type used over the field of galaxies to estimate the corrections to
                     from elliptcitiy to shear.
        
    OUTPUT :
         OUTFILE : A fitsfile named outfile with the same fields as corrected_moments, except with
                   two extra fields. GAMMA1 and GAMMA2, the two componenets of the estimated shear.

    '''
    #Need to filter as i determine mean quantities that
    #shoudlnt be used from bad galaxies
    signal_noise = corrected_moments['FLUX_AUTO'] / \
      corrected_moments['FLUXERR_AUTO']
  
    uncut_ell_sqr = corrected_moments['e1']*corrected_moments['e1'] + \
          corrected_moments['e2']*corrected_moments['e2']

    uncor_size = np.sqrt( 0.5*(corrected_moments.xx_uncorrected + \
                               corrected_moments.yy_uncorrected))
          
    good = np.zeros(len(corrected_moments.x))
    good[ galaxies ] = 1
    good[ (corrected_moments.xx + corrected_moments.yy < 0)] = 0
    good[ (uncut_ell_sqr > 2 ) ]  = 0
    good[ (uncor_size < size_cut[0] )]  = 0 
    good[ (uncor_size > size_cut[1] )]  = 0
    good[( corrected_moments.MAG_AUTO < mag_cut[0] )]  = 0 
    good[( corrected_moments.MAG_AUTO > mag_cut[1] )]  = 0
    good[ (signal_noise < signal_noise_cut)]  = 0 
    
    #good[ corrected_moments.nExposures < expThresh ] = 0
    good[  (~np.isfinite(corrected_moments.xx)) ] = 0
    good[  corrected_moments.prob != 0 ] = 0

    
    momc = corrected_moments[good == 1]
    size = np.sqrt( 0.5*(momc.xx + momc.yy))

    nObjects=len(momc.xx)

    weight = momc['radius']
    weight[ momc['radius'] < min_rad] = min_rad
    
    beta = 1./(2.*momc['gal_size']**2*(momc['shear']**2+weight**2))
    
    u1 = beta*(-momc['xxxx']+momc['yyyy']) 
    u2 = -2.*beta*(momc['xxxy']+momc['xyyy'])
    
    gal_lambda=beta*(momc['xxxx']+2.*momc['xxyy']+momc['yyyy'])
  
    
    ellipticity_sqr = momc['e1']**2+momc['e2']**2
    e_dot_u = momc['e1']*u1+momc['e2']*u2
    e_cross_u = momc['e1']*u2-momc['e2']*u1

    if stat_type == 'mean':
        #These are the mean G1, G2
        G2 = 0.5*np.nanmean(e_cross_u)
        G1 = 2-np.nanmean(ellipticity_sqr)-\
            0.5*np.nanmean(gal_lambda)-\
            0.5*np.nanmean(e_dot_u)
    elif stat_type =='median':
        #The median
        G2 = 0.5*np.nanmedian(e_cross_u)
        
        G1 = 2-np.nanmedian(ellipticity_sqr)- 0.5*np.nanmedian(gal_lambda)-0.5*np.nanmedian(e_dot_u)
    else:
        raise ValueError("Stat type not recognised")
    
    gamma1=momc['e1']/G1/rhodes_factor
    gamma2=momc['e2']/G1/rhodes_factor

    fits_cols = []
    for iName in momc.columns.names:
        fits_cols.append( py.Column(name=iName, format=momc[iName].dtype, array=momc[iName] ) )
        
    newcol = [ py.Column(name='gamma1', format=gamma1.dtype, array=gamma1),
               py.Column(name='gamma2', format=gamma2.dtype, array=gamma2) ]

 
    hdu = py.BinTableHDU.from_columns(fits_cols + newcol)
    hdu.writeto(outfile, clobber=True)

