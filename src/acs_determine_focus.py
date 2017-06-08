import star_galaxy_separation as sgs
import numpy as np
import measure_moms as mm
import acs_model_e as acs_model
import drizzle_position as dp
import os as os
import pyfits as py
import directories
def acs_determine_focus_metric( true, model ):
    '''
    SO A is the flag to check if the quadrupole moments
    match
    '''

    dof = np.max( [np.float((len(model.e1)-1)), 1])
    goodness_of_fit = np.sum( (true['e1']-model.e1)**2+(true['e2']-model.e2)**2) /dof
    
    return goodness_of_fit


# **************************************************************************

def acs_determine_focus( unknown_focus_image,
                         observed_moments,
                         drizzle_file,
                         wavelength,
                         data_dir=None,
                         psf_model_dir=None,
                         pixel_scale=0.03,
                         r_match=600):
    '''
    ;+      
    ; NAME:
    ;      ACS_DETERMINE_FOCUS
    ;
    ; CATEGORY:
    ;      Reduction of ACS COSMOS data.
    ;
    ; PURPOSE:
    ;      Decides the actual focus value of HST during observations, due to 
    ;      thermal fluctuations that change the distance between the primary
    ;      and secondary mirrors. The offset from the nominal focus position
    ;      is returned, in microns.
    ;
    ; INPUTS:
    ;      moms : Meaured catalogue of moments from the image
    ;
    ; OPTIONAL INPUTS:
    ;      CATALOGUES - Structure containing model TinyTim PSF catalogues. Can
    ;                   be set for speed. If not set, ;read from disk on first run.
    ;
    ; OUTPUTS:
    ;      Focus
    ;
    ; KEYWORDS:
    ;      SINGLE : If set work out the focus for a single exposure
    ;               this should be set to the name of the single exp
    ;
    ; EXAMPLE USE:
    ;      acs_determine_focus, cluster, filter, results
    ;
    ; MODIFICATION HISTORY:
    ;      Feb 05 - All shape moments needed by RRG used, by RM.
    ;      Jan 05 - Written by Richard Massey.
    ;      Aug 13 - Changed by DRH for single exposure
    ;      Sep 19 - Dont use the positions from single image
    ;-
    '''

 
    dirs = directories.return_dirs()

    
    #I need to get the positions of just the stars
    galaxies, stars = sgs.star_galaxy_separation( observed_moments, \
                                                      restore=True,\
                                                    savefile=dirs.data_dir+'/galStar.locus')
    
    #Filter the stars out
    observed_moments_stars = observed_moments[ stars ]
    
    
    #Now I need measure the moments of the stars in the individual image
    #But the x and y here are in the frame of the drizzled frame
    #So need to change these to position in the new fram
    
    image_name = unknown_focus_image.split('/')[-1] 
    inframe_stars = observed_moments_stars[observed_moments_stars[image_name+'_INFRAME'] == 1]
    

    if not os.path.isfile( unknown_focus_image[:-5]+'_uncor.cat'):
        mm.measure_moms( unknown_focus_image,  'NOCAT_NEED',
                    unknown_focus_image[:-5]+'_uncor.cat',
                    object_catalogue=inframe_stars,
                    xGal=inframe_stars[image_name+'_X_IMAGE'],
                    yGal=inframe_stars[image_name+'_Y_IMAGE'],
                    mult=2, min_rad=6)
        
    star_moments = py.open( unknown_focus_image[:-5]+'_uncor.cat' )[1].data
    
    #Now get an array of moments for all the possible focus positions
    model_e, focus = acs_model.acs_model_e(star_moments[image_name+'_X_IMAGE'],
                                           star_moments[image_name+'_Y_IMAGE'], 
                                           wavelength=wavelength )

    #Number of focus positions
    n_focus, nobjects = model_e.xx.shape
    focus = np.arange(16)-10
    average_distance = np.sum(model_e.offset_model, axis=0)/n_focus

    
    # Select only those stars with a suitably close model (since we're not interpolating)
    close_match = np.arange(nobjects)[average_distance < r_match]
    
    n_stars=len(close_match)
    if n_stars < 2:
        print "No stars found in field "
        global_focus=0.
        global_focus_error=100.
        focus = raw_input("No stars are found so please either input focus or cancel : ")
    
    # Tabulate model ellipticities
    model_e1=model_e.e1[:,close_match]
    model_e2=model_e.e2[:,close_match]

    # Tabulate measured ellipticities
    true = star_moments[close_match]
        
     
  #Find best-fit focus
  

    best_fit_focus_indiv = np.zeros( (2, n_stars ))

    for i in xrange(len(close_match)):
        chisq=np.zeros(n_focus) #Absolute chi squared
        for f in xrange(n_focus):
        #Tabulate model ellipticities
            model = moments( 1 )
            model['e1'][0] = model_e.e1[f,close_match][i]
            model['e2'][0] = model_e.e2[f,close_match][i]
            model['x'][0] = model_e.x[close_match][i]
            model['y'][0] = model_e.y[close_match][i]
            model['xx'][0] = model_e.xx[f,close_match][i]
            model['xy'][0] = model_e.xy[f,close_match][i]
            model['yy'][0] = model_e.yy[f,close_match][i]
            model['xxxx'][0] = model_e.xxxx[f,close_match][i]
            model['xxxy'][0] = model_e.xxxy[f,close_match][i]
            model['xxyy'][0] = model_e.xxyy[f,close_match][i]
            model['xyyy'][0] = model_e.xyyy[f,close_match][i]
            model['yyyy'][0] = model_e.yyyy[f,close_match][i]
            
            chisq[f]=acs_determine_focus_metric(true[i], model)
            

        best_fit_focus_indiv[ 0, i ] = focus[ np.argmin( chisq )]
        best_fit_focus_indiv[ 1, i ] = np.min(chisq)

    if n_stars < 5:
        new_best_fit_focus =  best_fit_focus_indiv[0, np.argmin(best_fit_focus_indiv[1,:])]
        #Just for the purposes of the test! MAKE SURE I DELTE

    else:
        best_fit_focus = np.median( best_fit_focus_indiv[0,:] )
        #If the median is between two take the best fit out of those two
        if np.round( best_fit_focus ) != best_fit_focus:
            floor_chi = best_fit_focus_indiv[1, best_fit_focus_indiv[0,:] ==  \
                                                     np.floor( best_fit_focus )]
            ceil_chi = best_fit_focus_indiv[1, best_fit_focus_indiv[0,:] ==  \
                                                     np.ceil( best_fit_focus )]
                                                     
            if (len(ceil_chi) == 0) | (len(floor_chi) == 0):
                chi_weighted_mean = np.sum( best_fit_focus_indiv[0,:]/best_fit_focus_indiv[1,:])/\
                  np.sum(1./best_fit_focus_indiv[1,:])
                
                new_best_fit_focus = np.round(chi_weighted_mean)
             
            else:
                min_floor_chi = np.min( floor_chi )
                min_ceil_chi = np.min( ceil_chi )
                focii = np.append( np.floor( best_fit_focus ), np.ceil( best_fit_focus ))
                chi = np.append( min_floor_chi, min_ceil_chi )
                new_best_fit_focus = focii[ np.argmin( chi ) ]
        else:
            
            return best_fit_focus
                                                                                                
                                                                          
    #Store the best one
    
    return new_best_fit_focus


class moments( dict ):

    def __init__( self, n_objects ):

        self.__dict__['x'] = np.zeros(n_objects)
        self.__dict__['y'] = np.zeros(n_objects)
        self.__dict__['e1'] = np.zeros(n_objects)
        self.__dict__['e2'] = np.zeros(n_objects)
        self.__dict__['xx'] = np.zeros(n_objects)
        self.__dict__['xy'] = np.zeros(n_objects)
        self.__dict__['yy'] = np.zeros(n_objects)
        self.__dict__['xxxx'] = np.zeros(n_objects)
        self.__dict__['xxxy'] = np.zeros(n_objects)
        self.__dict__['xxyy'] = np.zeros(n_objects)
        self.__dict__['xyyy'] = np.zeros(n_objects)
        self.__dict__['yyyy'] = np.zeros(n_objects)



    def keys(self):
        return self.__dict__.keys()

    def __getitem__(self, key): 
        return self.__dict__[key]
