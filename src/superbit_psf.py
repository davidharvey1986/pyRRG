from scipy.interpolate import SmoothBivariateSpline
from .acs_3dpsf import moments
import numpy as np
from matplotlib import pyplot as plt


def superbit_psf( galaxy_moms, star_moms, degree=2 ):
    
    radius = np.sqrt( ( galaxy_moms.xx + galaxy_moms.yy)/2.)
    
    moms = moments( galaxy_moms['x'], galaxy_moms['y'], galaxy_moms['radius'], degree )
    star_ell = np.sqrt( star_moms['e1']**2 + star_moms['e2']**2 )
  
    for iMom in moms.keys():
        
        if iMom in ['x','y','radius','degree']:
            continue
            

        no_ninety_nines = (star_moms[iMom] != -99) & ( star_moms[iMom] < np.quantile( star_moms[iMom], 0.68)) \
            & ( star_ell < np.quantile( star_ell, 0.90))
        
                                                      
        interpolate_fct = SmoothBivariateSpline( star_moms['x'][no_ninety_nines], star_moms['y'][no_ninety_nines], 
                                                star_moms[iMom][no_ninety_nines], kx=degree, ky=degree)
        
        moms[iMom][:] = interpolate_fct.ev(galaxy_moms['x'], galaxy_moms['y'] )[:][:]
    
    
    fig, ax = plt.subplots( 3, 1)
    scale=1.
    
    
    model_e = np.sqrt( moms['e1']**2 + moms['e2']**2)
    model_pa = np.arctan2( moms['e2'], moms['e1'])/2.
    
    ax[0].quiver(moms['x'], moms['y'], 
                 model_e*np.cos(model_pa), 
                 model_e*np.sin(model_pa)  , scale=scale)   
    
    true_e = np.sqrt( star_moms['e1']**2 + star_moms['e2']**2)
    true_pa = np.arctan2( star_moms['e2'], star_moms['e1'])/2.
    
    ax[1].quiver(star_moms['x'][no_ninety_nines],star_moms['y'][no_ninety_nines], 
                true_e[no_ninety_nines]*np.cos(true_pa[no_ninety_nines]), 
                 true_e[no_ninety_nines]*np.sin(true_pa[no_ninety_nines])  , scale=scale) 
    
    ax[2].hist(star_ell[ no_ninety_nines ], bins=np.linspace(0,1,100) )
    
    '''
    ax[2].plot(star_moms['e1'], moms['e1'], '.')
    ax[3].plot(star_moms['e2'], moms['e2'], '.')

    ax[2].set_ylim(0.,0.2)
    ax[2].set_xlim(0.,0.2)
    ax[2].plot([0,0.2],[0,0.2],'-')
    
    ax[3].set_ylim(0.,0.2)
    ax[3].set_xlim(0.,0.2)
    ax[3].plot([0,0.2],[0,0.2],'-')
    '''
    ax[0].set_title("Degree = %i" % degree)
    
    plt.savefig("psf.pdf")
    
    return moms