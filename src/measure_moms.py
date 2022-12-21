from . import mmm as mmm
from astropy.io import fits
import numpy as np
import sys
import RRGtools as at
import json
from tqdm import tqdm

def measure_moms(fits_image, sex_catalog, outfile, verbose=False, quiet=False, **kwargs):
                     
    '''
    ;
    ; NAME:                  rrg_measure_moms
    ;
    ; PURPOSE:
    ;    Uses the RRG method to measure the moments of objects in a fits image.  
    
    ;
    ; INPUTS:
    ;   fits_image- fits image with single extension
    ;   sex_catalog-SExtractor catalog in fits format
    ;   outfile- output file to store IDL array of variables
    ;
    ; OUTPUTS:
    ; stores moments and other info in the output file as an idl structure  
    ;
    ; KEYWORD PARAMETERS:
    ;   width- gaussian window function width, if not input is calculated from area
    ;   saturation,bad_value are saturation level and value of pixels not to use (outside image)
    ;   cut_off- cut_off for measuring objects default is 2.5x width
    ;   mult- multiplier used to find gaussian width from area
    ;   min_rad- stellar radius
    ;   startx, starty- lower left corner of image in pixel values, default to zero
    ;   verbose set to False to supress messages
    ;   stellar set to 1 to evaluate all objects with the same radius (min_rad)
    ;   weight=use a weight image to determine good/bad pixels (yes=1 0=no/default)
    ;   wt_ext= extension of weight image in fits file
    ;  nocenter= if 1 then uses Sextractor centers and doesnt centroid
    ;  name_weight_im=name of the weight image
    ;  sgm_im Sextractor segmentation image to get sky background 
    ;  souces = predetermined sources to find
    ;MODIFICATION HISTORY:
    ;    30 Dec 2002 jrhodes
    ;    11 March 2003 jrhodes
    

    IF N_PARAMS(0) LT 1 THEN BEGIN
    PRINT, 'rrg_measure_moms,fits_image,sex_catalog,outfile,width=width,saturation=saturation,'
    PRINT, '    badval=badval,cut_off=cut_off,fieldback=fieldback,mult=mult,min_rad=min_rad,'
    PRINT, '    startx=startx,starty=starty,verbose=verbose,stellar=stellar,'
    PRINT, '    weight=weight,wt_ext=wt_ext,nocenter=nocenter,name_weight_im=name_weight_im'
    return
    ENDIF
    ;
    '''
    if 'saturation' not in kwargs.keys():
        kwargs['saturation'] = 800000000
    if 'bad_val' not in kwargs.keys():
        kwargs['bad_val'] = -99
    if 'cut_off' not in kwargs.keys():
        kwargs['cut_off'] = 2.5
    if 'mult' not in kwargs.keys():
        kwargs['mult'] = 1               
    if 'min_rad' not in kwargs.keys():
        kwargs['min_rad'] = 1.5           
    if 'start' not in kwargs.keys():
        startx, starty = 0, 0
    else:
        startx, starty = kwargs['start']
        
    if 'weight_file' not in kwargs.keys():
        kwargs['weight_file'] =  None  
    if 'nocenter' not in kwargs.keys():
        kwargs['nocenter'] = 0. 
            
    if 'sgm_im' not in kwargs.keys():
        kwargs['sgm_im'] = 'null'
   
    if 'object_catalogue' not in kwargs.keys():
        object_catalogue = None
    else:
        object_catalogue = kwargs['object_catalogue']
     
    if 'regfile' not in kwargs.keys():
        kwargs['regfile'] = None

    if 'sky' not in kwargs.keys():
        skymed, skysd, skysw = None, None, None
    else:
        skymed, skysd, skysw = kwargs['sky']['skymed'], kwargs['sky']['skysd'], kwargs['sky']['skysw']
        
    if 'return_moms' not in kwargs.keys():
        kwargs['return_moms'] = True
        
    if 'error' not in kwargs.keys():
        kwargs['error'] =  0.01
        
    if 'min_it' not in kwargs.keys():
        kwargs['min_it'] =  500
        
                    
                    
                    
                    
    img_file = fits.open( fits_image )


    if 'fits_extension' not in kwargs.keys():
        kwargs['fits_extension'] = 0
    if 'expTimeName'  not in kwargs.keys():
        kwargs['expTimeName'] = 'EXPTIME'
    
    imhead = img_file[kwargs['fits_extension']].header
    img = img_file[kwargs['fits_extension']].data
    
    exp_time = imhead[kwargs['expTimeName']]
 
    #What is mmm?
    if (skymed is None )| \
        (skysd is None):
        skymed, skysd, skysw = mmm.mmm(img)
    if verbose:
        print((' % f skymed and %f skysd' % (skymed,skysd)))
    
    if object_catalogue is None:
        cat_file=fits.open(sex_catalog)
        object_catalogue = cat_file[1].data
        header=cat_file[1].header
     
    nGalaxies = len(object_catalogue['RA'])

    
    try:
        if not 'xGal' in kwargs:
            xGal = object_catalogue['X_IMAGE']
        else:
            xGal = kwargs['xGal']
        if not 'yGal' in kwargs:
            yGal = object_catalogue['Y_IMAGE']
        else:
            yGal = kwargs['yGal']
    except:
        xGal, yGal = at.deg2pix( fits_image, object_catalogue['RA'], object_catalogue['DEC'] )

  
    ysize= imhead['NAXIS2']
    xsize= imhead['NAXIS1']

    if kwargs['weight_file'] is None:
        wt_image = np.ones( img.shape)
    else:
        wt_image = fits.open( kwargs['weight_file'] )[0].data



    #The background for each galaxy. if none take from sex cata;pgue
    try:
        back=object_catalogue.BACKGROUND
    except:
        back =  np.zeros( nGalaxies) + skysd

    

    if (kwargs['sgm_im'] != 'null'):
        seg_file = fits.open(kwargs['sgm_im'])
        seg = seg_file[0].data
        shdr = seg_file[0].header
        '''
        sel=where((seg eq 0)and (wt_image gt 0))
        mmm,img(sel),seg_sky,seg_skysd
        back=replicate(seg_sky,num) 
        '''
    if not 'radius' in kwargs:
        area=np.pi*object_catalogue['A_IMAGE']*object_catalogue['B_IMAGE']
        radius = np.sqrt( area / np.pi )*kwargs['mult']
    else:
        radius = kwargs['radius']
    
    radius[ radius < kwargs['min_rad'] ] = kwargs['min_rad']

        
    cut_rad = radius*kwargs['cut_off'] # cut off radius of object
    offedge=0 #;number of objects off the endge initially
    centerprob=0 # number of objects with centroiding problems
    badpix_prob=0 #number of objects with bad pixels

    
    galaxy_moments = moms( nGalaxies, radius=radius )

    print("Measuring Object Moments")
    for i in tqdm(range( nGalaxies )):
        
        #following changed by jrhodes to account for different indexing in SExtractor and IDL
       

        deltax=1.0
        deltay=1.0
        count=1
        sum_int=0
        int_moms = moms( 1 )
        xc=xGal[i]-0.5
        yc=yGal[i]-0.5
        
        blank=0
        go_on=1
        badpix_centroid='no'
    
        if (xc-cut_rad[i]-1.< startx) | \
            (yc-cut_rad[i]-1 < starty) | \
            (xc+cut_rad[i]+1 > xsize ) | \
            (yc+cut_rad[i]+1 > ysize ):

            if verbose:
                print((' %i %f %f too close to edge at iteration 1\n' % \
                          (i,xGal[i],yGal[i])))
            offedge=offedge+1
            go_on=0
            galaxy_moments.prob[i] += 1
            

        if kwargs['nocenter']  == 0:
            while  ((np.abs(deltax) > kwargs['error']) |  \
              (np.abs(deltay) > kwargs['error'])) & \
               (count < kwargs['min_it']) & (go_on == 1) :
                #These are needed for all moments
                #Cut out a postage stamp of the image
                begin_x = np.round(xc-cut_rad[i]-1).astype(int)
                begin_y = np.round(yc-cut_rad[i]-1).astype(int)
                end_x = np.round(xc+cut_rad[i]+2).astype(int)
                end_y = np.round(yc+cut_rad[i]+2).astype(int)
                #Make sure the postage stamp doesnt go out the end
                #of the imgae.
                begin_x = np.max( [0, begin_x])
                begin_y = np.max( [0, begin_y])
                end_x = np.min( [img.shape[1], end_x])
                end_y = np.min( [img.shape[0], end_y])
                
            
            
                xvec = np.arange(begin_x, end_x)
                yvec = np.arange(begin_y, end_y)
            
                #this will be used throughotu the loop as well
                uncut_xgrid, uncut_ygrid = np.meshgrid( xvec, yvec)
                
                #Dist will also be used throughout the loop
                uncut_dist = np.sqrt( (uncut_xgrid - xc + 0.5)**2 + \
                                    (uncut_ygrid - yc + 0.5)**2)
            
                #Also cut the images at the cut_radi
                dist = uncut_dist[uncut_dist < cut_rad[i]]
                xgrid = uncut_xgrid[uncut_dist < cut_rad[i]]
                ygrid = uncut_ygrid[uncut_dist < cut_rad[i]]
    
                #These will be used throughout the loop
                uncut_postage_stamp = img[begin_y:end_y, begin_x:end_x]
                postage_stamp = uncut_postage_stamp[uncut_dist < cut_rad[i]]
            
                uncut_postage_weight = wt_image[begin_y:end_y, begin_x:end_x]
                postage_weight = uncut_postage_weight[uncut_dist < cut_rad[i]]
    
                #This might be
                g_f = np.exp(-(dist**2)/(2.*radius[i]**2))
                weight_gal = g_f*( np.abs(postage_stamp -  back[i]))
                sum_int = np.sum( weight_gal )
                

                #Then finally change the dist ]
                #find the centroid
 
   
                #changed from above line for SNAP sims jrhodes 7/09/07
                #used below for SNAP sims
                checkx=xc
                checky=yc
                
                #Get the centre
                xc = np.sum( weight_gal*(xgrid+0.5))/sum_int
                yc = np.sum( weight_gal*(ygrid+0.5))/sum_int
                
                count += 1
                deltax=xc-checkx
                deltay=yc-checky;
            
                if ( xc-cut_rad[i]-1 < startx) | \
                    ( yc-cut_rad[i]-1 < starty) | \
                    ( xc+cut_rad[i]+1 > xsize ) | \
                    ( yc+cut_rad[i]+1 > ysize) :
                    if verbose:
                        print(( ' %i %f %f too close to edge at iteration %i' % \
                            (i,xGal[i],yGal[i],(count-1)) ))
                   
                    go_on=0
                    galaxy_moments.prob[i] += 1
    
                if (badpix_centroid == 'yes'):
                    go_on=0
                    if verbose:
                        print(( '%i %f %f  bad pix in centroid at iteration %f' %
                                    (i,xGal[i],yGal[i],(count-1)) ))
                                
                    galaxy_moments.prob[i] += 4
                    badpix_prob += 1
                        #End this if
                    #Ends the xc_ci of
                #this ends the while
            if (count > 100):
                centerprob=centerprob+1
                if verbose:
                    print(( '%i %f %f Too many centering iterations %i' % \
                            (i,xGal[i],yGal[i],count)))
                galaxy_moments.prob[i] += 2
                
            if np.sqrt( (xc-xGal[i])**2 + (yc-yGal[i])**2) > radius[i]:
                go_on=0
                centerprob += 1
                if verbose:
                    print(( '%i %f %f Centroid shift too large ' % (i,xGal[i],yGal[i]) ))
                
                galaxy_moments.prob[i] += 2
      
        
            #this ends the centre befor the while
            #Now back into normal main galaxy loop
        badpix_mom='no'    
    # check for saturated pixels and bad pixels
        int_moms.x = xc
        int_moms.y = yc
       
        if go_on == 1:


            if (np.any(postage_stamp > kwargs['saturation'])) | \
                (np.any(postage_stamp == kwargs['bad_val'])) | \
                (np.any(postage_weight == 0)):
                go_on = 0
                galaxy_moments.prob[i] += 4
                badpix_mom='yes'
                
            if badpix_mom == 'yes':
                badpix_prob += 1 
                
                print(('%i %f %f Bad pixel(s) in centroiding' %\
                        (i,xGal[i],yGal[i])))

        
    #find moments and center error
    

        if go_on == 1:
            rel_xgrid = xgrid + 0.5 - int_moms.x
            rel_ygrid = ygrid + 0.5 - int_moms.y
            
            int_moms.xx = np.sum(rel_xgrid*rel_xgrid*weight_gal) / sum_int
            int_moms.yy = np.sum(rel_ygrid*rel_ygrid*weight_gal) / sum_int
            int_moms.xy = np.sum(rel_xgrid*rel_ygrid*weight_gal) / sum_int



            int_moms.xxxx = np.sum(weight_gal*(rel_xgrid**4)) / sum_int
            int_moms.xxxy = np.sum(weight_gal*(rel_xgrid**3)*rel_ygrid)  / sum_int
            int_moms.xxyy = np.sum(weight_gal*(rel_xgrid**2)*(rel_ygrid**2)) / sum_int
            int_moms.xyyy=  np.sum(weight_gal*(rel_ygrid**3)*rel_xgrid)/ sum_int
            int_moms.yyyy = np.sum(weight_gal*(rel_ygrid**4)) / sum_int

      
            #find moment errors and covariances

            #Error in internsity
            I_err_sq =  skysd**2+(postage_stamp-back[i])/( exp_time) / sum_int**2

            #Error in x and y
            sum_xxdI = np.sum(g_f*g_f*I_err_sq*rel_xgrid*rel_xgrid) / sum_int
            sum_yydI = np.sum(g_f*g_f*I_err_sq*rel_ygrid*rel_ygrid) / sum_int
            int_moms.error.x = np.sqrt(sum_xxdI)/sum_int
            int_moms.error.y = np.sqrt(sum_yydI)/sum_int

            #Error in second order
            int_moms.error.xx = np.sum(g_f*g_f*I_err_sq*(rel_xgrid**2-int_moms.xx)**2)/sum_int
            int_moms.error.yy = np.sum(g_f*g_f*I_err_sq*(rel_ygrid**2-int_moms.yy)**2)/sum_int
            int_moms.error.xy = np.sum(g_f*g_f*I_err_sq*(rel_xgrid*rel_ygrid-int_moms.xy)**2)/sum_int

            #error in 4th order
            int_moms.error.xxyy = np.sum(g_f*g_f*I_err_sq*\
                (rel_xgrid*rel_xgrid-int_moms.xx)*\
                (rel_ygrid*rel_ygrid-int_moms.yy))/sum_int**2
                          
            int_moms.error.xxxy = np.sum(g_f*g_f*I_err_sq*\
                (rel_xgrid*rel_xgrid-int_moms.xx)*\
                (rel_ygrid*rel_xgrid-int_moms.xy))/sum_int**2
                        
            int_moms.error.xyyy = np.sum(g_f*g_f*I_err_sq*\
                          (rel_xgrid*rel_xgrid-int_moms.xy)*\
                          (rel_ygrid*rel_ygrid-int_moms.yy))/sum_int**2
            
            
        else:
            
            int_moms.xx, int_moms.yy, int_moms.xy = [-99, -99, -99]
            int_moms.xxxx, int_moms.xxxy =  [0, 0]
            int_moms.xxyy, int_moms.xyyy, int_moms.yyyy = [0, 0, 0]
            int_moms.error.xx, int_moms.error.yy, int_moms.error.xy = [0, 0, 0]
            int_moms.error.xxxx, int_moms.error.xxxy =  [0, 0]
            int_moms.error.xxyy, int_moms.error.xyyy, int_moms.error.yyyy = [0, 0, 0]
            xc_err, yc_err = [0, 0]




        galaxy_moments.append( i, int_moms )
        
    #Append some needed things
    #Can just append the RA as this is not the same as the pixel position
    #Need to find the re-centred RA

    recentred_ra, recentred_dec = \
              at.pix2deg( fits_image, galaxy_moments.x, galaxy_moments.y, extension=kwargs['fits_extension'])

    galaxy_moments.ra = recentred_ra
    galaxy_moments.dec = recentred_dec
   
    galaxy_moments.calc_e1e2( mult_rad=kwargs['mult'])     
    #log this stuff as i want this for the star galaxy separation
    galaxy_moments['skymed'][:] = skymed
    galaxy_moments['skysd'][:] = skysd
    galaxy_moments['skysw'][:] = skysw
    galaxy_moments['exp_time'][:] = exp_time


    galaxy_moments.write_to_fits( object_catalogue, outfile )

    if kwargs['regfile'] is not None:
        galaxy_moments.fits_to_ellipse( kwargs['regfile'] )
    
    if kwargs['return_moms']:
        return fits.open(outfile)[1].data

class moms( dict ):

    def __init__( self, ngalaxies, radius=None, mk_error=True ):
        self.__dict__['x'] = np.zeros( ngalaxies )
        self.__dict__['y'] = np.zeros( ngalaxies )
        self.__dict__['xx'] = np.zeros( ngalaxies )
        self.__dict__['yy'] = np.zeros( ngalaxies )
        self.__dict__['xy'] = np.zeros( ngalaxies )
        self.__dict__['xxxx'] = np.zeros( ngalaxies)
        self.__dict__['yyyy'] = np.zeros( ngalaxies)
        self.__dict__['xxyy'] = np.zeros( ngalaxies)
        self.__dict__['xyyy'] = np.zeros( ngalaxies)
        self.__dict__['xxxy'] = np.zeros( ngalaxies)
        self.__dict__['prob'] = np.zeros( ngalaxies)

        
        #more information for the SVM star galaxy divider
        self.__dict__['skymed'] = np.zeros( ngalaxies)
        self.__dict__['skysd'] =  np.zeros( ngalaxies)
        self.__dict__['skysw'] =  np.zeros( ngalaxies)
        self.__dict__['exp_time'] =  np.zeros( ngalaxies)

        if radius is not None:
            self.__dict__['radius'] = radius

        if mk_error:
            self.error = moms( ngalaxies, mk_error=False )
    
    
    def append( self, index, imom ):

        mom_names = list(self.keys())
        for i in mom_names:
            if i == 'radius':
                continue
            self[i][index] = imom[i]
            if i == 'error':
                error_names = list(imom.error.keys())
                for i in error_names:
                    self.error[i][index] = imom.error[i]


    def calc_e1e2( self, mult_rad=1 ):
        self.__dict__['e1']=(self.xx-self.yy)/ (self.xx+self.yy)
        self.e2=2.*self.xy/(self.xx+self.yy)

        self.__dict__['e1_err'] = np.sqrt( ((self.xx+self.yy)**(-2))*\
                               (self.error.xx**2*(1-self.e1)**2 + \
                                self.error.yy**2*(1+self.e1)**2-\
                                2.*(1-self.e1*self.e1)*self.error.xxyy))
                                
        self.__dict__['e2_err'] = np.sqrt( ((self.xx+self.yy)**(-2))*\
                               ((self.error.xx**2+self.error.yy**2+\
                                 2.*self.error.xxyy)*self.e2**2+\
                                 4.*(self.error.xy**2-self.e2*\
                                     (self.error.xxxy*self.error.xyyy))))

        self.__dict__['ell'] = np.sqrt( self.e1**2 + self.e2**2)
        self.__dict__['pa'] = 0.5*np.arctan2(self.e2, self.e1)*180./np.pi
        
        self.__dict__['gal_size'] = np.sqrt( (self.xx + self.yy)/2.)
        self.prob[ self.xx + self.yy < 0 ] += 5
        self.prob[ self.ell > 1.0 ] += 6
        self.__dict__['a'] = self.gal_size*(1.+self.ell)
        self.__dict__['b'] = self.gal_size*(1.-self.ell)
        

    def write_to_fits( self, sex_catalog, filename):
        '''
        Append the moments to input sex catalogue
        '''
        
        mom_names = list(self.keys())
        fits_columns = []
        
        for iColumn in mom_names:
            if iColumn == 'error':
                error_names = list(self.error.keys())
                for iColumn_err in error_names:
                    fits_columns.append(fits.Column( name=iColumn_err+'_err',
                                                   format=self.error[iColumn_err].dtype,
                                                   array=self.error[iColumn_err] ))
            else:
                fits_columns.append(fits.Column( name=iColumn, \
                                               format=self[iColumn].dtype, \
                                               array=self[iColumn] ))
        
        for iCatalog in sex_catalog.columns.names:
            if (iCatalog in mom_names) | ('err' in iCatalog):
                continue
            else:
                fits_columns.append(fits.Column( name=iCatalog,
                                                   format=sex_catalog[iCatalog].dtype,
                                                   array=sex_catalog[iCatalog] ))

            
        fits_table = fits.BinTableHDU.from_columns( fits_columns )
        fits_table.writeto( filename, overwrite=True,output_verify='ignore' )



    def fits_to_ellipse( self, filename):

        regionFile = open( filename, "wb")
        regionFile.write('# Region file format: DS9 version 4.1\n')
        regionFile.write('# Filename: dummy.fits\n')
        regionFile.write("global color=green dashlist=8 3 width=1 font='helvetica 10 normal roman' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n")
        regionFile.write("physical\n")
        
        for i in range(len(self.ra)):
            
            regionFile.write('ellipse(%0.4f,%0.4f,%0.4f,%0.4f,%0.4f)\n' %
                             (self.x[i], self.y[i], self.a[i], self.b[i], self.pa[i]))


        

    def keys(self):
        return list(self.__dict__.keys())
    
    def __getitem__(self, key): 
        return self.__dict__[key]
