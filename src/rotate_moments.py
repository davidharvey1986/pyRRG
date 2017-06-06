import numpy as np
import combination as com
import copy as cp


def rotate_moments( moms, angle):
    
    '''

    ;PURPOSE : TO TAKE THE MOMENTS OF A PSF AND ROTATE THEM SO THEY
    ;          ARE IN THE FRAME OF REFERENCE OF THE DRIZZLED IMAGE
    
    ;         THIS SHOULD ROTATE ALL MOMENTS AND ELLIPTICITIES
    
    
    ;INPUTS :
    ;          MOMS : The structure with the moments to be rotated
    ;                 - xx, xy, yy
    ;                 - xxxx, xxxy, xxyy, xyyy, yyyy
    ;                 - e1, e2
    ;          ANGLE : Angle for the moments to be rotated through
    
    
    ;OUTPUTS :
    ;         MOMS_rotated : The same structure are moms, except
    ;                        moments now rotated
    
    ;METHOD : It will use the equation (36) found in 

    ;http://www.opticsinfobase.org/DirectPDFAccess/71A5781F-D53F-D862-FC6B61687B808297_57703/josa-70-8-920.pdf?da=1&id=57703&seq=0&mobile=no
    '''
    
    #rotate the ellipticities
    #need to make sure angles are in radians
    angle=angle*np.pi/180.
  
    n_galaxies=len(moms.xx)
    #Rotate the quadropole moments first
    #in all case
    
    #Create a matrix of the various quadropole moments
    # a=[[00,10,20],[01,11,21],[02,12,22]]
    #mu_jk = [[0,0,moms.xx],[0,moms.xy,0],[moms.yy,0,0]]
    #create a n array such that the mu_jk correspond to xx, xy, and yy
    mu_jk = np.zeros((3,3,n_galaxies))
    mu_jk[2,0,:]=moms.xx 
    mu_jk[1,1,:]=moms.xy 
    mu_jk[0,2,:]=moms.yy 
  
  
    mu_jk_rot=np.zeros((3, 3, n_galaxies))
    for j in xrange(3):
        for k in xrange(3):
            for r in xrange(j+1):
                for s in xrange(k+1):
                    #Only do those moments we care about
                    if j+k != 2:
                        continue
                    
                        
                    mu_jk_rot[j,k,:] += (-1.)**(k-s)*\
                        com.combination(j, r)*com.combination(k, s)*\
                        np.cos(angle)**(j-r+s)*\
                        np.sin(angle)**(k+r-s)*\
                        mu_jk[j+k-r-s,r+s,:]
              


    #Create a matrix of the various fourth order moments
    # a=[[00,10,20,30,40],$
    #    [01,11,21,31,41],$
    #    [02,12,22,32,42],$
    #    [03,13,23,33,43],$
    #    [04,14,24,34,44]]
    
    
    #mu_jk = [[0,0,moms.xx],[0,moms.xy,0],[moms.yy,0,0]]
    #create a n array such that the mu_jk correspond to xx, xy, and yy
    mu_jklm = np.zeros((5, 5, n_galaxies))
    mu_jklm[4,0,:] = moms.xxxx
    mu_jklm[3,1,:] = moms.xxxy
    mu_jklm[2,2,:] = moms.xxyy
    mu_jklm[1,3,:] = moms.xyyy
    mu_jklm[0,4,:] = moms.yyyy
    
    
    #then rotate them and put them in a new matrix
    mu_jklm_rot=np.zeros((5, 5, n_galaxies))
    for j in xrange(5):
        for k in xrange(5):
            for r in xrange(j+1):
                for s in xrange(k+1):
                    
 
                    if j+k != 4:
                        continue

                    mu_jklm_rot[j,k,:] += (-1)**(k-s)*\
                        com.combination(j, r)*com.combination(k, s)* \
                        np.cos(angle)**(j-r+s)*\
                        np.sin(angle)**(k+r-s)*\
                        mu_jklm[j+k-r-s,r+s, :]
              


    #e1 =(xx-yy)/(xx+yy) e2=2*xy/(xx+yy)
    e1_rotated= (mu_jk_rot[2,0,:]-mu_jk_rot[0,2,:])/\
              (mu_jk_rot[2,0,:]+mu_jk_rot[0,2,:])

    e2_rotated= (2.0*mu_jk_rot[1,1,:])/\
        (mu_jk_rot[2,0,:]+mu_jk_rot[0,2,:])

    e1_rotated[ (mu_jk_rot[2,0,:]+mu_jk_rot[0,2,:]) == 0] = 0
    e2_rotated[ (mu_jk_rot[2,0,:]+mu_jk_rot[0,2,:]) == 0] = 0
    
    MOMS_rotated= cp.deepcopy( moms )

    MOMS_rotated.xx = mu_jk_rot[2,0, :]
    MOMS_rotated.xy = mu_jk_rot[1,1, :]
    MOMS_rotated.yy = mu_jk_rot[0,2, :]
    MOMS_rotated.xxxx = mu_jklm_rot[4,0,:]
    MOMS_rotated.xxxy = mu_jklm_rot[3,1, :]
    MOMS_rotated.xxyy = mu_jklm_rot[2,2, :]
    MOMS_rotated.xyyy = mu_jklm_rot[1,3, :]
    MOMS_rotated.yyyy = mu_jklm_rot[0,4, :]
    MOMS_rotated.e1 = e1_rotated
    MOMS_rotated.e2 = e2_rotated
    MOMS_rotated.radius = moms.radius
    MOMS_rotated.degree = moms.degree
               

    return MOMS_rotated


  

