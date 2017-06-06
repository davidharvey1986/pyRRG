import mmm as mmm
import pyfits as py
import numpy as np

from matplotlib import pyplot as plt
import RRGtools as at
def measure_moms(fits_image, sex_catalog, outfile,
                     width=None,
                     saturation=800000000, badval=-99,
                     cut_off=2.5, fieldback=0, mult=1,
                     min_rad=1.5, startx=0, starty=0,
                     silent=0, weight_image=None,
                     wt_ext=0, nocenter=0,
                     sgm_im='null',
                     object_catalogue=None,
                     bad_val=-99, regfile=None,
                     skymed=None, skysd=None,
                     return_moms=True, **kwargs):
                     
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
    ;   fieldback -1 means use a field backgroup as calculated by IDL SKY, 0(default) uses object back from SEx
    ;   mult- multiplier used to find gaussian width from area
    ;   min_rad- stellar radius
    ;   startx, starty- lower left corner of image in pixel values, default to zero
    ;   silent set to 1 to supress messages
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
    PRINT, '    startx=startx,starty=starty,silent=silent,stellar=stellar,'
    PRINT, '    weight=weight,wt_ext=wt_ext,nocenter=nocenter,name_weight_im=name_weight_im'
    return
    ENDIF
    ;
    '''
                     
    img_file = py.open( fits_image )
    img = img_file[0].data
    imhead = img_file[0].header
    exp_time = imhead['EXPTIME']

    
    #What is mmm?
    if (skymed is None )| \
        (skysd is None):
        skymed, skysd, skysw = mmm.mmm(img)
    
    print(' % f skymed and %f skysd' % (skymed,skysd))
    
    if object_catalogue is None:
        cat_file=py.open(sex_catalog)
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

    if weight_image is None:
        wt_image = np.ones( img.shape)
    else:
        wt_image = py.open( weight_image )[0].data



    #The background for each galaxy. if none take from sex cata;pgue
    try:
        back=object_catalogue.BACKGROUND
    except:
        back =  np.zeros( nGalaxies) + skysd

    

    if (sgm_im != 'null') & (fieldback == 1):
        seg_file = py.open(sgm_im)
        seg = seg_file[0].data
        shdr = seg_file[0].header
        '''
        sel=where((seg eq 0)and (wt_image gt 0))
        mmm,img(sel),seg_sky,seg_skysd
        back=replicate(seg_sky,num) 
        '''
    if not 'radius' in kwargs:
        area=np.pi*object_catalogue['A_IMAGE']*object_catalogue['B_IMAGE']
        radius = np.sqrt( area / np.pi )*mult
    else:
        radius = kwargs['radius']
    
    radius[ radius < min_rad ] = min_rad

    
    if width is not None:
        radius=np.zeros(nGalaxies) + width
        
    cut_rad = radius*cut_off # cut off radius of object
    offedge=0 #;number of objects off the endge initially
    centerprob=0 # number of objects with centroiding problems
    badpix_prob=0 #number of objects with bad pixels

    
    galaxy_moments = moms( nGalaxies, radius=radius )
    
    for i in xrange( nGalaxies ):
        
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

            if not silent:
                print(' %i %f %f too close to edge at iteration 1\n' % (i,xGal[i],yGal[i]))
            offedge=offedge+1
            go_on=0
            galaxy_moments.prob[i] += 1
            

        if nocenter == 0:
            while  ((np.abs(deltax) > 0.01) |  \
              (np.abs(deltay) > 0.01)) & \
               (count < 500) & (go_on == 1) :
                #These are needed for all moments
                #Cut out a postage stamp of the image
                begin_x = np.round(xc-cut_rad[i]-1).astype(int)
                begin_y = np.round(yc-cut_rad[i]-1).astype(int)
                end_x = np.round(xc+cut_rad[i]+2).astype(int)
                end_y = np.round(yc+cut_rad[i]+2).astype(int)
            
            
                xvec = np.arange(begin_x, end_x)
                yvec = np.arange(begin_y, end_y)
            
                #this will be used throughotu the loop as well
                uncut_xgrid, uncut_ygrid = np.meshgrid( xvec, yvec)
                
                #Dist will also be used throughout the loop
                uncut_dist = np.sqrt( (uncut_xgrid - xc + 0.5)**2 + \
                                    (uncut_ygrid-yc+0.5)**2)
            
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
                weight_gal = g_f*( postage_stamp -  back[i])
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
                    if not silent:
                        print( ' %i %f %f too close to edge at iteration %i' % \
                            (i,xGal[i],yGal[i],(count-1)) )
                   
                    go_on=0
                    galaxy_moments.prob[i] += 1
    
                if (badpix_centroid == 'yes'):
                    go_on=0
                    if not silent:
                        print( '%i %f %f  bad pix in centroid at iteration %f' %
                                    (i,xGal[i],yGal[i],(count-1)) )
                                
                    galaxy_moments.prob[i] += 4
                    badpix_prob += 1
                        #End this if
                    #Ends the xc_ci of
                #this ends the while
            if (count > 100):
                centerprob=centerprob+1
                if not silent:
                    print( '%i %f %f Too many centering iterations %i' % \
                            (i,xGal[i],yGal[i],count))
                galaxy_moments.prob[i] += 2
                
            if np.sqrt( (xc-xGal[i])**2 + (yc-yGal[i])**2) > mult*radius[i]:
                go_on=0
                centerprob += 1
                if not silent:
                    print( '%i %f %f Centroid shift too large ' % (i,xGal[i],yGal[i]) )
                
                galaxy_moments.prob[i] += 2
      
        
            #this ends the centre befor the while
            #Now back into normal main galaxy loop
        badpix_mom='no'    
    # check for saturated pixels and bad pixels
        int_moms.x = xc
        int_moms.y = yc
       
        if go_on == 1:


            if (np.any(postage_stamp > saturation)) | \
                (np.any(postage_stamp == bad_val)) | \
                (np.any(postage_weight == 0)):
                go_on = 0
                galaxy_moments.prob[i] += 4
                badpix_mom='yes'
                
            if badpix_mom == 'yes':
                badpix_prob += 1 
                
                print('%i %f %f Bad pixel(s) in centroiding' %\
                        (i,xGal[i],yGal[i]))

        
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
                (rel_ygrid*rel_xgrid-int_moms.yy))/sum_int**2
                        
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
      at.pix2deg( fits_image, galaxy_moments.x, galaxy_moments.y)

    galaxy_moments.ra = recentred_ra
    galaxy_moments.dec = recentred_dec
   
    galaxy_moments.calc_e1e2( mult_rad=mult)     

    galaxy_moments.write_to_fits( object_catalogue, outfile )

    if regfile is not None:
        galaxy_moments.fits_to_ellipse( regfile)
    
    if return_moms:
        return py.open(outfile)[1].data

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
        if radius is not None:
            self.__dict__['radius'] = radius

        if mk_error:
            self.error = moms( ngalaxies, mk_error=False )
    
    
    def append( self, index, imom ):

        mom_names = self.keys()
        for i in mom_names:
            if i == 'radius':
                continue
            self[i][index] = imom[i]
            if i == 'error':
                error_names = imom.error.keys()
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
        
        mom_names = self.keys()
        fits_columns = []
        
        for iColumn in mom_names:
            if iColumn == 'error':
                error_names = self.error.keys()
                for iColumn_err in error_names:
                    fits_columns.append(py.Column( name=iColumn_err+'_err',
                                                   format=self.error[iColumn_err].dtype,
                                                   array=self.error[iColumn_err] ))
            else:
                fits_columns.append(py.Column( name=iColumn, \
                                               format=self[iColumn].dtype, \
                                               array=self[iColumn] ))
        
        for iCatalog in sex_catalog.columns.names:
            if (iCatalog in mom_names) | ('err' in iCatalog):
                continue
            else:
                fits_columns.append(py.Column( name=iCatalog,
                                                   format=sex_catalog[iCatalog].dtype,
                                                   array=sex_catalog[iCatalog] ))

            
        fits_table = py.BinTableHDU.from_columns( fits_columns )
        fits_table.writeto( filename, clobber=True )



    def fits_to_ellipse( self, filename):

        regionFile = open( filename, "wb")
        regionFile.write('# Region file format: DS9 version 4.1\n')
        regionFile.write('# Filename: dummy.fits\n')
        regionFile.write("global color=green dashlist=8 3 width=1 font='helvetica 10 normal roman' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n")
        regionFile.write("physical\n")
        
        for i in xrange(len(self.ra)):
            
            regionFile.write('ellipse(%0.4f,%0.4f,%0.4f,%0.4f,%0.4f)\n' %
                             (self.x[i], self.y[i], self.a[i], self.b[i], self.pa[i]))


        

    def keys(self):
        return self.__dict__.keys()
    
    def __getitem__(self, key): 
        return self.__dict__[key]
