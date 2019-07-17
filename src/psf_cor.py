import pyfits as py
import numpy as np
import os as os
import glob as glob
from . import drizzle_position as dp
from . import acs_determine_focus as adf
from . import acs_3dpsf as acs_3dpsf
from scipy.io.idl import readsav 
from . import rotate_moments as rm
import copy as cp
from . import directories
import sys
def psf_cor(    mom_file,
                outfile,
                drizzle_file,
                wavelength,
                mult=1, min_rad=1.5, chip=1,
                constantpsf=0, mscale=0, 
                order=3,
                n_chip=2):
    '''
    ;
    ; NAME:                  rrg_psf_cor
    ;
    ; PURPOSE:
    ;    Uses the RRG method to correct moments for PSF effects.  
    ;  corrects using a 3rd order polynomial fit to the moments
    ;the PSF moments should be weighted, this programs corrects for the weighting of the PSF moments
    ;min_rad is the size of the weighting function used to measure the PSF
    ;
    ; INPUTS:
    ;   mom_file- moment file IDL structure produced by rrg_measure_mom
    ;   sex_catalog-SExtractor catalog in fits format
    ;   outfile- output file to store IDL array of variables
    ;   mult multiplier used to find gaussian width from area   
    ;   min_rad- stellar radius
    ;   coeff_file text file with the coefficients of the PSF correction
    ;OUTPUTS:
    ; stores moments and other info in the output file as an idl structure  
    ;

    
    ; MODIFICATION HISTORY:
    ;   
    ;    April 2003 jrhodes
    ;  May 2005 jrhodes added mscale option, only works with constantpsf
    ;  July 2011 dharvey modified to allow for rotated images,(see section
    ;  for more details)
    
    ; TO DO :    ;CHECK THAT ANGLES ARE CORRECT  PLEASE
    
    ;mom_file has following structure
    
    ;moms={x:x,y:y,xx:xx,yy:yy,xy:xy,xxxx:xxxx,xxxy:xxxy,$
    ;xxyy:xxyy,xyyy:xyyy,yyyy:yyyy,radius:radius,sn:sn,back:back,$
    ;class:class,area:area,area_ab:area_ab,flags:flags,$
    ;a:a,b:b,theta:theta,mag:mag,prob:prob}
    '''
    dirs = directories.return_dirs( )

    moms = py.open(mom_file)[1].data

    radius = np.sqrt( ( moms.xx + moms.yy)/2.)
    
    
 
    sigma =  cp.copy(moms.radius[moms['galStarFlag']==1])
   
    sigma[ sigma < min_rad ] = min_rad
    

    w2=(1./(sigma*sigma));
    w4=(1./(2.*sigma*sigma*sigma*sigma))
    



    #; Evaluate moments of PSF
    #;I will work out the focus myself! muwagahaha
    #;Before i get the psf moments i need to get the scat catalogue, which
    #;is found by running tintim_make_scat.pro (for tiny_tim models)
    print('Getting the psf models from tinytim, cheers Tim!')
    #need to think about this
    #tinytim_make_scat, data_dir=dirs.model_dir, wavelength=filter[0], scat=scat

    scat = readsav( dirs.psf_model_dir+'/TinyTim'+wavelength+'.scat' )['scat']

    #so this function interpolates.
 
    '''
    ;Okay so the way i am going to get the correct psf for each position
    ;in the drizzled image is.
    ; 1. Create a grid of positions which covers the entire field of view
    ;of the drizzles image
    ;2. For each position in the drizzled image, work out how many images
    ;cover this position
    ;3. For each of the images covering this position, work out the psf
    ;4. Then take the average of the moments for each 
    '''

    images = glob.glob( dirs.data_dir+'/j*_drz_sci.fits')
    if len(images) == 0:
        raise ValueError('Cant find single exposures of field')

    nImages = len(images)


    #Now get the positions in the drizzle frame of ref in the individual
    #frame of ref

    print("Getting position of stars & galaxies in each exposure")

    momsWithDrizzlePosition =  \
      dp.drizzle_position( drizzle_file, images,  moms, dataDir=dirs.data_dir)
    galaxy_moms = cp.copy(momsWithDrizzlePosition[momsWithDrizzlePosition['galStarFlag'] == 1])
    star_moms = cp.copy(momsWithDrizzlePosition[momsWithDrizzlePosition['galStarFlag'] == 0])

    uncorrected_xx = galaxy_moms.xx
    uncorrected_yy = galaxy_moms.yy
    
   
   
    
    #Also get the Orientations in terms of the drizzled image, not
    #WCS


    #Now loop through each position in the psf grid, 
    #and check if the positionsreffram(i,3,*) is non zero! and interpolate
    #to that x and y and then average moms
    #check that the position in the grid is in, then loop through each
    #image, interpolating the psf to this point, 

    #PsfMoms is an vector of many classes of moments
    nGalaxies=len(galaxy_moms.x)
    psf_moms = moments( galaxy_moms.x, galaxy_moms.y, nGalaxies )

    FocusArray = np.zeros(nImages)

    sys.stdout.write("\n")
    for iImage in range(nImages):
        sys.stdout.write("Getting PSF for image: %i/%i\r" % \
                                 (iImage+1,nImages))
        sys.stdout.flush()
        #Which positions are in the cluster frame
        iImage_name = images[iImage].split('/')[-1][0:8]
        inFrame = galaxy_moms[iImage_name+'_INFRAME'] == 1

        #before i determine the psf moments i need to get the focus
        #position of the image in question

        #So get the focus position by fitting the true image stars to the
        #model
        
        focus = adf.acs_determine_focus(  images[iImage], star_moms, \
                                              drizzle_file, wavelength)

        #Just keep track of the focii i have used through out
        FocusArray[iImage] = focus

        #For all the points in the main drizzled field that are within
        #iImage, interpolate the psf from the ref fram of the single
        #image  to the X,Y of the drizzled image
        
       
        iPsfMoms=\
          acs_3dpsf.acs_3dpsf( galaxy_moms[iImage_name+'_X_IMAGE'][inFrame], 
                                   galaxy_moms[iImage_name+'_Y_IMAGE'][inFrame],
                                    np.zeros(len(galaxy_moms[iImage_name+'_INFRAME'][inFrame]))+focus, \
                                    radius, scat, degree=[3,2,2] )
        #now rotate the moments according to the angle in orient
        iPsfMomsRot = rm.rotate_moments( iPsfMoms, galaxy_moms[iImage_name+'_ORIENTAT'][inFrame])
        
        #CHECK THAT ANGLES ARE CORRECT HERE PLEASE
        
        mom_names = list(iPsfMoms.keys())
        for iMom in mom_names:
            
            #I need to now rotate each moment according to the axis orient 
            #with the drizzled image before stacking.
            if (iMom != 'degree') & \
                (iMom != 'nExposures') &\
                (iMom != 'radius') &\
                (iMom != 'x') & (iMom != 'y'):
                psf_moms[iMom][inFrame] += iPsfMomsRot[iMom]
   
        #then keep count how many images per position
        psf_moms['nExposures'][ inFrame] += 1

        
        #then give the position the value of the averaged psf_moms.

    #Save the focus array
    focuslist = open(dirs.data_dir+'/FocusArray.txt', "w")
    
    for i in range(nImages):
        ExpName = images[i].split('/')[-1].split('_')[0]
        focuslist.write( "%s %3.1f \n" % \
                             (ExpName, FocusArray[i]))


    #Now take the mean of each moment
    


    for iMom in mom_names:
        if (iMom != 'degree') & \
          (iMom != 'nExposures') & \
          (iMom != 'radius') &\
          (iMom != 'x') & (iMom != 'y'):
            psf_moms[iMom] /= psf_moms['nExposures']
            psf_moms[iMom][psf_moms['nExposures'] == 0] = 0.
        
   
        
    #Refind e1 and e2, assuming we want  <q11>-<q22>/<q11>+<q22> not <e1>
    psf_moms.e1=(psf_moms.xx-psf_moms.yy)/ \
        (psf_moms.xx+psf_moms.yy)
    psf_moms.e2=2.*psf_moms.xy/\
        (psf_moms.xx+psf_moms.yy)
    #plot the psf model
    '''
    ps, DrizzleFile+'/Psf_Model.ps', /carre, /psfont
    plot, psf_moms.x, psf_moms.y, title='PSF MODEL FOR '+cluster, $
    xtitle='X [PIXELS]', ytitle='Y [PIXELS]', /iso, /nodata
    plt_evec, psf_moms.x, psf_moms.y, $
    psf_moms.e1, psf_moms.e2, /e1e2, $
    xscale=1000, yscale=1000
    unps
    '''
    #and continue on as per usual...

    pxxw=psf_moms.xx
    pxyw=psf_moms.xy
    pyyw=psf_moms.yy
    pxxxx=psf_moms.xxxx
    pxxxy=psf_moms.xxxy
    pxxyy=psf_moms.xxyy
    pxyyy=psf_moms.xyyy
    pyyyy=psf_moms.yyyy
    w=min_rad
    
    


    #correction for weighgting

    pxx=(pxxw-(pxxw*pyyw+pxxw*pxxw-pxxxx-pxxyy)/(2.*w**2))
    pyy=(pyyw-(pyyw*pyyw+pyyw*pxxw-pyyyy-pxxyy)/(2.*w**2))
    pxy=(pxyw-(pxyw*pxxw+pxyw*pyyw-pxyyy-pxxxy)/(2.*w**2))
  
    


    trace=pxx+pyy
    a=pxx-trace/2
    b=pxy
    c=pxy
    d=pyy-trace/2
    e=pxxxx
    f=pxxxy
    z=pxxyy
    h=pxyyy
    k=pyyyy




    #So the problem is that i am timesing a which is a 
    # [1,10000] psf_vector, we moms which is a vector [a,b]
    ixxa=a*(1-w2*2*galaxy_moms.xx+w4*\
                (galaxy_moms.xxxx-galaxy_moms.xx*galaxy_moms.xx));
    ixxb=b*(w4*(galaxy_moms.xxxy-galaxy_moms.xx*galaxy_moms.xy));
    ixxc=c*(-w2*(2*galaxy_moms.xy)+\
                w4*(galaxy_moms.xxxy-galaxy_moms.xx*galaxy_moms.xy));
    ixxd=d*(w4*(galaxy_moms.xxyy-galaxy_moms.xx*galaxy_moms.yy));




    iyya=a*(w4*(galaxy_moms.xxyy-galaxy_moms.xx*galaxy_moms.yy));
    iyyb=b*(w4*(galaxy_moms.xyyy-galaxy_moms.yy*galaxy_moms.xy));
    iyyc=c*(-w2*(2*galaxy_moms.xy)+\
                w4*(galaxy_moms.xyyy-galaxy_moms.yy*galaxy_moms.xy));
    iyyd=d*(1-w2*2*galaxy_moms.yy+\
                w4*(galaxy_moms.yyyy-galaxy_moms.yy*galaxy_moms.yy));
    
    ixya=a*(-w2*galaxy_moms.xy+\
                w4*(galaxy_moms.xxxy-galaxy_moms.xy*galaxy_moms.xx));
    ixyb=b*(1-w2*galaxy_moms.xx+\
                w4*(galaxy_moms.xxyy-galaxy_moms.xy*galaxy_moms.xy));
    ixyc=c*(-w2*galaxy_moms.yy+\
                w4*(galaxy_moms.xxyy-galaxy_moms.xy*galaxy_moms.xy));
    ixyd=d*(-w2*galaxy_moms.xy+\
                w4*(galaxy_moms.xyyy-galaxy_moms.xy*galaxy_moms.yy));
    
    ixx_corr=galaxy_moms.xx-ixxa-ixxb-ixxc-ixxd;
    iyy_corr=galaxy_moms.yy-iyya-iyyb-iyyc-iyyd;
    ixy_corr=galaxy_moms.xy-ixya-ixyb-ixyc-ixyd;



    shear=np.sqrt(0.5*(a+d+trace))
    
    if constantpsf:
         shear=np.zeros(number)+np.sqrt(0.5*(a(0)+d(0)+trace(0)))


    corrected_moments = moments( galaxy_moms.x, galaxy_moms.y, \
                                     len(galaxy_moms.x))

    gw=np.sqrt( (shear*shear*sigma*sigma)/(shear*shear+sigma*sigma));
    corrected_moments.xx = ((shear/gw)**4)*(ixx_corr-gw*gw);
    corrected_moments.yy = ((shear/gw)**4)*(iyy_corr-gw*gw);
    corrected_moments.xy = ((shear/gw)**4)*(ixy_corr);
    pxxc=a+trace/2
    pyyc=d+trace/2
    pxyc=b
    


    corrected_moments.xxxx = galaxy_moms.xxxx-e-\
      6*pxxc*galaxy_moms.xx+6*pxxc*pxxc;
    corrected_moments.xxxy = galaxy_moms.xxxy-f-\
      3*(pxyc*galaxy_moms.xx+pxxc*galaxy_moms.xy)+6*pxyc*pxxc;
    corrected_moments.xxyy = galaxy_moms.xxyy-z-\
      pxxc*galaxy_moms.yy-pyyc*galaxy_moms.xx-\
      4*pxyc*galaxy_moms.xy+2*pxxc*pyyc+4*pxyc*pxyc;
    corrected_moments.xyyy = galaxy_moms.xyyy-h-\
      3*(pxyc*galaxy_moms.yy+pyyc*galaxy_moms.xy)+6*pxyc*pyyc;
    corrected_moments.yyyy = galaxy_moms.yyyy-\
      k-6*pyyc*galaxy_moms.yy+6*pyyc*pyyc;

    corrected_moments['e1'] = (corrected_moments.xx-corrected_moments.yy)/\
      (corrected_moments.xx+corrected_moments.yy)
    corrected_moments['e2'] = (2*corrected_moments.xy)/\
      (corrected_moments.xx+corrected_moments.yy)
  
    #Those moments that were originally zero and -99 make them again hgere
    for i in list(corrected_moments.keys()):
        if i in moms.columns.names:
            corrected_moments[i][galaxy_moms[i] == -99] = -99
            corrected_moments[i][galaxy_moms[i] == 0] = 0
            galaxy_moms[i] = corrected_moments[i]
            
    

    
    galaxy_moms['gal_size'] = np.sqrt( (corrected_moments.xx +corrected_moments.yy)/2.)

    #SOMETHIGN STUPID TO SEE WHAT HAPPENS
    #SOme weird shit going on, i have to write this out, and then
    #read it back in
    galaxy_moms = writeAndRemoveUnusedColums( galaxy_moms)
    

    newcol = [ py.Column(name='shear', format=shear.dtype, array=shear),
               py.Column(name='nExposures', format=psf_moms.nExposures.dtype, \
                         array=psf_moms.nExposures),
                py.Column('xx_uncorrected', format=galaxy_moms.xx.dtype, array=uncorrected_xx),
                py.Column('yy_uncorrected', format=galaxy_moms.yy.dtype, array=uncorrected_yy)]
    
    orig_cols = galaxy_moms.columns
    new_cols = py.ColDefs(newcol)
    
    hdu = py.BinTableHDU.from_columns(orig_cols+new_cols)
    hdu.writeto( outfile, clobber=True)

class moments( dict ):

    def __init__(self, x, y, n_objects):
        self.__dict__['x'] = x
        self.__dict__['y'] = y
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
        self.__dict__['nExposures'] = np.zeros(n_objects)

    def __setitem__(self, key, item): 
        self.__dict__[key] = item
        
    def keys(self):
        return list(self.__dict__.keys())

    def __getitem__(self, key): 
        return self.__dict__[key]


def   writeAndRemoveUnusedColums( moments):

    momentNames = moments.columns.names
    columns = []
    for i in momentNames:
        if (not 'INFRAME' in i) & \
            (not 'fits_X_IMAGE' in i) & \
            (not 'fits_Y_IMAGE' in i) & \
             (not 'ORIENTAT' in i ):
            iColumn = \
              py.Column(i, format=moments[i].dtype, \
                            array=moments[i])
            columns.append(iColumn)
                
    new_cols = py.ColDefs(columns)
    
    hdu = py.BinTableHDU.from_columns(new_cols)            
    hdu.writeto('galaxies.fits',clobber=True)
    return py.open('galaxies.fits')[1].data
