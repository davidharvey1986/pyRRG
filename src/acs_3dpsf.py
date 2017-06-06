import numpy as np
import acs_map_xy as acs_map
def acs_3dpsf_basisfunctions( degree, x, y, focus ):




    # Generate relevant basis functions
    n_stars=np.max( np.array([len(x),len(y),len(focus)]))

    basis_function_order=np.zeros((1,3)) # All zeros
    for k in xrange(degree[2]+1):
        for j in xrange(degree[1]+1):
            for i in xrange(degree[0]+1):
                if (i+j+k > 0) & ((i+j) <= np.max(degree[0:2])):
                    
                    basis_function_order=np.vstack((basis_function_order, [i,j,k]))
    
    n_basis_functions= basis_function_order.shape[0]
    basis_function_value = np.zeros( (n_basis_functions, n_stars))
    
    for i in xrange(n_basis_functions):
        basis_function_value[i,:] = x**basis_function_order[i,0]*\
            y**basis_function_order[i,1] * \
            focus**basis_function_order[i,2]
                    
 

    return basis_function_value




# **********************************************************************
# **********************************************************************
# **********************************************************************

def acs_3dpsf_fit( scat, degree=np.array([3,2,2]), 
                   mag_cut=np.array([20.5,22]),
                   e_cut=1, size_cut=np.array([-np.inf,3])
                   ):

    # Fit the PSF from data in a SCAT catalogue
    # F814 I magnitude catalogue cut 
    

    if len(degree) < 3 :
        print "DEGREE must be 3D"
    degree[ degree > 0 ] = np.min(degree[ degree > 0 ])

    # Find the line dividing CCDs 1 and 2
    ccd_boundary = acs_map.acs_map_xy( np.array([0, 4095, 0, 4095]),
                                    np.array([2047, 2047, 2048, 2048]),
                                    pixel_scale=scat.pixscale)
    
    x1=np.mean([ccd_boundary.x[0],ccd_boundary.x[2]])
    x2=np.mean([ccd_boundary.x[1],ccd_boundary.x[3]])
    y1=np.mean([ccd_boundary.y[0],ccd_boundary.y[2]])
    y2=np.mean([ccd_boundary.y[1],ccd_boundary.y[3]])
    
    ccd_boundary_x1=np.mean([ccd_boundary.x[0],ccd_boundary.x[2]])
    ccd_boundary_x2=np.mean([ccd_boundary.x[1],ccd_boundary.x[3]])
    ccd_boundary_y1=np.mean([ccd_boundary.y[0],ccd_boundary.y[2]])
    ccd_boundary_y2=np.mean([ccd_boundary.y[1],ccd_boundary.y[3]])
    
    ccd_boundary_m=(ccd_boundary_y2-ccd_boundary_y1)/(ccd_boundary_x2-ccd_boundary_x1)
    ccd_boundary_c=ccd_boundary_y1-ccd_boundary_m*ccd_boundary_x1  


    # Find the centre of each CCD
    ccd_centre = acs_map.acs_map_xy( np.array([2048,2048]),
                                     np.array([3072,1024]), pixel_scale=scat.pixscale)


    # Select only the well-behaved stars
    good= np.isfinite(scat.field_focus[0][scat.field_id[0]]) & \
        np.isfinite(scat.e1_uncor_unrot[0]) & \
        np.isfinite(scat.e2_uncor_unrot[0]) & \
        np.isfinite(scat.xx_uncor[0]) & \
        np.isfinite(scat.xy_uncor[0]) & \
        np.isfinite(scat.yy_uncor[0]) & \
        np.isfinite(scat.xxxx_uncor[0]) & \
        np.isfinite(scat.xxxy_uncor[0]) & \
        np.isfinite(scat.xxyy_uncor[0]) & \
        np.isfinite(scat.xyyy_uncor[0]) & \
        np.isfinite(scat.yyyy_uncor[0])

    n_good = len(np.arange( len( good ))[good])
    print "Found a total of "+str(len(scat.x[0]))+" real stars, of which "+str(n_good)+" look well-behaved"
    

    # Store quantities to be fitted in local variables
    x=scat.x[0][good]
    y=scat.y[0][good]
    focus=scat.field_focus[0][scat.field_id[0]][good]
    ixx=scat.xx_uncor[0][good]
    ixy=scat.xy_uncor[0][good]
    iyy=scat.yy_uncor[0][good]
    ixxxx=scat.xxxx_uncor[0][good]
    ixxxy=scat.xxxy_uncor[0][good]
    ixxyy=scat.xxyy_uncor[0][good]
    ixyyy=scat.xyyy_uncor[0][good]
    iyyyy=scat.yyyy_uncor[0][good]
    e1=scat.e1_uncor_unrot[0][good]
    e2=scat.e2_uncor_unrot[0][good]
    
    # Work on each CCD separately
    init_coeffs_flag = True
    for ccd in xrange(2):
    
        # Report which CCD is being considered
        if ccd +1  == 1: 
            in_ccd = np.arange(len(y))[ y >= ccd_boundary_m*x+ccd_boundary_c]
            n_in_CCD = len(in_ccd)
        if ccd + 1 == 2:
            in_ccd = np.arange(len( y))[ y < ccd_boundary_m*x+ccd_boundary_c]
            n_in_CCD = len(in_ccd)
  
  
        if n_in_CCD > 0:
            #Compute matrix necessary for matrix inversion
            
            print "Fitting moments of "+str(n_in_CCD)+" real stars in CCD#"+str(ccd+1)
            basis_function_value=acs_3dpsf_basisfunctions(degree, 
                                                  x[in_ccd]-ccd_centre.x[ccd], 
                                                  y[in_ccd]-ccd_centre.y[ccd], 
                                                  focus[in_ccd])
           
            ls_matrix = np.dot( np.linalg.inv(np.dot(basis_function_value, basis_function_value.T)), basis_function_value)
        # Create global arrays to contain the answers
    
        n_basis_functions=np.shape(np.array(ls_matrix))[0]

        if init_coeffs_flag:
            acs_3dpsf_coeffs=basis_coeffs( ccd_centre,
                                    ccd_boundary_m, ccd_boundary_c,
                                    n_basis_functions, degree )
            init_coeffs_flag = False 
        print "Using "+str(n_basis_functions)+" basis functions"
        
        
        # Fit data to basis functions using least-squares inversion
        #these are all matrices
        acs_3dpsf_coeffs.ixx_fit[ccd, :]   = np.dot(ls_matrix ,ixx[in_ccd])
        acs_3dpsf_coeffs.ixy_fit[ccd, :]   = np.dot(ls_matrix , ixy[in_ccd])
        acs_3dpsf_coeffs.iyy_fit[ccd, :]   = np.dot(ls_matrix , iyy[in_ccd])
        acs_3dpsf_coeffs.ixxxx_fit[ccd, :] = np.dot(ls_matrix , ixxxx[in_ccd])
        acs_3dpsf_coeffs.ixxxy_fit[ccd, :] = np.dot(ls_matrix , ixxxy[in_ccd])
        acs_3dpsf_coeffs.ixxyy_fit[ccd, :] = np.dot(ls_matrix , ixxyy[in_ccd])
        acs_3dpsf_coeffs.ixyyy_fit[ccd, :] = np.dot(ls_matrix , ixyyy[in_ccd])
        acs_3dpsf_coeffs.iyyyy_fit[ccd, :] = np.dot(ls_matrix , iyyyy[in_ccd])
        acs_3dpsf_coeffs.e1_fit[ccd, :]    = np.dot(ls_matrix , e1[in_ccd])
        acs_3dpsf_coeffs.e2_fit[ccd, :]    = np.dot(ls_matrix , e2[in_ccd])
        
    else:
        print "No real stars found in CCD#"+str(ccd +1)
        


    return acs_3dpsf_coeffs



# **********************************************************************
# **********************************************************************
# **********************************************************************

def acs_3dpsf_reconstruct( acs_3dpsf_coeffs, x, y, focus, radius=None):


    # Create arrays to contain the final answer
    n_galaxies=np.max( np.array([len(x), len(y), len(focus)]) )
    
    if len(focus) == 1:
        focus_local = np.zeros(len(n_galaxies)) + focus
    else:
        focus_local=focus
        
    print "Found a total of "+str(n_galaxies)+" galaxies"
    if radius is None:
        radius=np.zeros(len(n_galaxies))+6
        
    moms=moments( x, y, radius[:n_galaxies],
                     acs_3dpsf_coeffs.degree )

        
    for ccd in xrange(2):
        #Report which CCD is being considered

        if ccd +1  == 1: 
            in_ccd = np.arange(len( y))[ y >= acs_3dpsf_coeffs.ccd_boundary_m*x+acs_3dpsf_coeffs.ccd_boundary_c]
            n_in_CCD = len(in_ccd)
        if ccd + 1 == 2:
            in_ccd = np.arange(len( y))[ y < acs_3dpsf_coeffs.ccd_boundary_m*x+acs_3dpsf_coeffs.ccd_boundary_c]
            n_in_CCD = len(in_ccd)
        
        if n_in_CCD > 0:
            print "Interpolating model PSF moments to the position of "+str(n_in_CCD)+" galaxies in CCD#"+str(ccd+1)

            #Fit the PSF
            basis_function_value=acs_3dpsf_basisfunctions(acs_3dpsf_coeffs.degree[0], \
                                                           x[in_ccd]-acs_3dpsf_coeffs.ccd_centre.x[ccd], \
                                                           y[in_ccd]-acs_3dpsf_coeffs.ccd_centre.y[ccd], \
                                                           focus_local[in_ccd] )
            
            moms.xx[in_ccd]   = np.dot(acs_3dpsf_coeffs.ixx_fit[ccd, :], basis_function_value)
            moms.xy[in_ccd]   =  np.dot(acs_3dpsf_coeffs.ixy_fit[ccd, :], basis_function_value)
            moms.yy[in_ccd]   =  np.dot(acs_3dpsf_coeffs.iyy_fit[ccd, :], basis_function_value)
            moms.xxxx[in_ccd] =  np.dot(acs_3dpsf_coeffs.ixxxx_fit[ccd, :], basis_function_value)
            moms.xxxy[in_ccd] =  np.dot(acs_3dpsf_coeffs.ixxxy_fit[ccd, :], basis_function_value)
            moms.xxyy[in_ccd] =  np.dot(acs_3dpsf_coeffs.ixxyy_fit[ccd, :], basis_function_value)
            moms.xyyy[in_ccd] =  np.dot(acs_3dpsf_coeffs.ixyyy_fit[ccd, :], basis_function_value)
            moms.yyyy[in_ccd] =  np.dot(acs_3dpsf_coeffs.iyyyy_fit[ccd, :], basis_function_value)
            moms.e1[in_ccd]   =  np.dot(acs_3dpsf_coeffs.e1_fit[ccd, :], basis_function_value)
            moms.e2[in_ccd]   =  np.dot(acs_3dpsf_coeffs.e2_fit[ccd, :], basis_function_value)
            
           
        else:
            print "No galaxies in CCD#"+str(ccd)

    # Work out PSF ellipticities at positions of galaxies properly. Tsk!
    moms.e1 = (moms.xx-moms.yy)/(moms.xx+moms.yy)
    moms.e2 = 2*moms.xy/(moms.xx+moms.yy)

    return moms


# **********************************************************************
# **********************************************************************
# **********************************************************************

def acs_3dpsf( x, y, focus, radius, scat, 
                acs_3dpsf_coeffs=None, 
                degree=np.array([3,2,2])):

    # Fit the PSF
    if acs_3dpsf_coeffs is None:
        acs_3dpsf_coeffs=acs_3dpsf_fit(scat, degree=degree)
    #Reconstruct the PSF
    acs_moms=acs_3dpsf_reconstruct(acs_3dpsf_coeffs, x, y, focus, radius)
    return acs_moms
  


class basis_coeffs:

    def __init__( self, ccd_centre, ccd_boundary_m, \
                  ccd_boundary_c, n_basis_functions, degree ):
        self.degree = degree,
        self.ccd_centre = ccd_centre
        self.ccd_boundary_m = ccd_boundary_m
        self.ccd_boundary_c = ccd_boundary_c
        self.ixx_fit = np.zeros((2,n_basis_functions))
        self.ixy_fit = np.zeros((2,n_basis_functions))
        self.iyy_fit = np.zeros((2,n_basis_functions))
        self.ixxxx_fit = np.zeros((2,n_basis_functions))
        self.ixxxy_fit = np.zeros((2,n_basis_functions))
        self.ixxyy_fit = np.zeros((2,n_basis_functions))
        self.ixyyy_fit = np.zeros((2,n_basis_functions))
        self.iyyyy_fit = np.zeros((2,n_basis_functions))
        self.e1_fit = np.zeros((2,n_basis_functions))
        self.e2_fit = np.zeros((2,n_basis_functions))
    
    
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


    def keys(self):
        return self.__dict__.keys()

    def __getitem__(self, key): 
        return self.__dict__[key]

