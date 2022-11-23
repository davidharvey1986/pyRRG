
from scipy.interpolate import RectBivariateSpline
import numpy as np


def interpolate_jwst_psf_moms( x, y, radius, scat, degree=3):
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
    grid_size = np.int(len(x_vector))

    
    psf_grid_size = (grid_size, grid_size)
    
    for iMom in moms.keys():
        
        if iMom in ['x','y','radius','degree']:
            continue
        z_vector = scat[iMom].reshape(psf_grid_size)
        
        interpolate_fct = RectBivariateSpline( x_vector, y_vector, z_vector, kx=degree, ky=degree)
        
        moms[iMom] = interpolate_fct.ev(x, y)


        
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
    