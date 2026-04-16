from astropy.io import fits
import RRGtools as at
import numpy as np

def remove_galaxy_members( 
    galaxy_filename, 
    shear_filename, 
    outfile=None,
    verbose=False) :
    '''
    Remove a fits table of galaxy members
    
    '''
    
    try:
        galaxy_cat = fits.open( galaxy_filename )[1].data
    except:
        raise ValueError("Input galaxy catalogue is not a fits file")
    
    try:
        shear_cat = fits.open( shear_filename )[1].data
    except:
        raise ValueError("Input shear catalogue is not a fits file")
    
    matched_cat = at.run_match( 
        galaxy_filename, 
        shear_filename, 
        search_rad=2.)[1].data

    nGalaxies = shear_cat.shape[0]

    if np.any('NUMBER' == np.array(matched_cat.dtype.names)):
        rm_gals = np.concatenate([ 
            np.arange(nGalaxies)[shear_cat['NUMBER'] == i]
            for i in matched_cat['NUMBER']
        ])
    else:

        rm_gals = np.concatenate([ 
            np.arange(nGalaxies)[shear_cat['NUMBER'] == i]
            for i in matched_cat['NUMBER_2a']
        ])

    if verbose:
        print(
            "Removing %i cluster members from the final catalogue "
             % (len(rm_gals) )
        )
    keep_gals = np.ones(nGalaxies)
    keep_gals[ rm_gals ] = 0
                        
    
    shear_cat = shear_cat[ keep_gals == 1]
    
    if outfile is None:
        fits.writeto( shear_filename, shear_cat, overwrite=True)
    else:
        fits.writeto( outfile, shear_cat, overwrite=True)
        
        
        