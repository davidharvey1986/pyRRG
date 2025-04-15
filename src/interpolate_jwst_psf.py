
from scipy.interpolate import RectBivariateSpline, NearestNDInterpolator
import numpy as np


def interpolate_jwst_psf_moms(
        x, y,
        radius,
        scat,
        degree=3,
        star_moms=None,
        verbose=False,
        interpolation='spline'):
    '''
    
    Use the same interpolation function as photoutils use for the webb psf grid
    evaulate model.
    
    
    See here fore more information:
    https://photutils.readthedocs.io/en/stable/_modules/photutils/psf/models.html#GriddedPSFModel
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RectBivariateSpline.html
    
    
    '''
    
    moms = moments( x, y, radius, degree)

    x_vector = np.unique(scat['X'])
    y_vector = np.unique(scat['Y'])
    grid_size = int(len(x_vector))

    
    psf_grid_size = (grid_size, grid_size)
    if star_moms is not None:

        calibrate_to_these = ( star_moms['xx'] > 0 )
        
        all_corrs = 1./np.array([
            np.median(scat['xx'])/np.median(star_moms['xx'][ calibrate_to_these ]),
            np.median(scat['yy'])/np.median(star_moms['yy'][ calibrate_to_these ])
            ])
            
        corr = np.sqrt(np.nanmean(all_corrs))
        
        if verbose:
            print("APPLING PSF MOMENT CORRECT OF %0.4f" % corr)
        if corr == 0:
            raise ValueError("Zero correction not valid")
        if ~np.isfinite(corr):
            raise ValueError("Non finite correction")
            
    else:
        corr = 1
    
    for iMom in moms.keys():
        
        if iMom in ['x','y','radius','degree']:
            continue

        z_vector = scat[iMom].reshape(psf_grid_size)

        if interpolation.lower() == 'spline':        
            interpolate_fct = RectBivariateSpline(
                x_vector,
                y_vector,
                z_vector,
                kx=degree, ky=degree)
            this_mom = interpolate_fct.ev(x, y)
            
        else:
            interpolate_fct = NearestNDInterpolator(
                (scat['X'], scat['Y']),
                scat[iMom]
            )
            this_mom = interpolate_fct(x, y)

            
            

        if 'e' in iMom.lower():
            moms[iMom] = this_mom
        else:
            moms[iMom] = this_mom*corr**len(iMom)

        
    return moms
    
    
    
    
    
    
class moments( dict ):

    def __init__(self, x, y, radius, degree ):
        n_objects = len(x)
        self.__dict__['x'] = x
        self.__dict__['y'] = y
        self.__dict__['e1']=np.zeros(n_objects)
        self.__dict__['e2']=np.zeros(n_objects)
        self.__dict__['xx']=np.zeros(n_objects)
        self.__dict__['xy']=np.zeros(n_objects)
        self.__dict__['yy']=np.zeros(n_objects)
        self.__dict__['xxxx']=np.zeros(n_objects)
        self.__dict__['xxxy']=np.zeros(n_objects)
        self.__dict__['xxyy']=np.zeros(n_objects)
        self.__dict__['xyyy']=np.zeros(n_objects)
        self.__dict__['yyyy']=np.zeros(n_objects)
        self.__dict__['radius'] = radius
        self.__dict__['degree'] = degree

    def __setitem__(self, key, item): 
        self.__dict__[key] = item
        
    def keys(self):
        return list(self.__dict__.keys())

    def __getitem__(self, key): 
        return self.__dict__[key]
    
