import numpy as np
import acs_limits as al
import RRGtools as at
import pyfits as py

def drizzle_position(      drizzle_file,
                           individual_files,
                           moments,
                           dataDir='.'):
    
    '''
    ;PURPOSE : RETURN THE POSITIONS OF EACH POINT IN THE DRIZZLED
    ;          IMAGE IN THE FRAME OF REFERENCE OF THE INDIVIDUAL EXP
    
    ;INPUTS : CLUSTER & FILTER : THE NAME AND FILTER OF CLUSTER
    ;         ASTROMETRYARR : ARRAY OF THE ASTROMETRY STRUCTURE FOR EACH
    ;         IMAGE 
    
    ;KEYWORD : NGRID : SQRT OF THE NUMBER O FPOINTS IN THE FINAL PSF GRID 
    ;          PSFGRID[X,Y] = THE GRID ON WHICH THE DRIZZLED IMAGE LAID
    ;          DRIZZLEDIMAGEORIENTATION : THE ORIENTATION OF THE DRIZZLED
    ;                                     IMAGE WRT WCS
    ;          [X,Y]_IMAGE : VECTORS OF THE INPUT X AND Y IF NOT WANT A 
    ;                        GRID IN THE FRAME OF THE DRIZZLED IMAGE
    ;          
    
    ;OPTIONAL OUTPUTS : NUM_EXPOSURES : NUMBER OF EXPOSURES EACH GALAXY HAS
    
    ;OUTPUTS : InDrizzleFrame : ARRAY OF NIMAGES X NX X NY X 3 
    ;                           where [*, 0, *] = xpositions
    ;                                 [*, 0, *] = ypositions
    ;                                 [*, 2, *] = is it in the image?
    
    ;UPDATES : 7/8/2013 : ALLOW RETURN OF DRIZZLED IMAGE ORIENTATION TO
    ;                     WCS
    ;          24/9/2013 : Returning keyword, num_exposures
    '''
    
    nGalaxies = len(moments.X_IMAGE)
    nImages = len( individual_files )

    InDrizzleFrame= np.zeros( (nImages, 3, nGalaxies**2))
    
    drizzle_obj = py.open(dataDir+'/'+drizzle_file)
    ImageData = drizzle_obj[0].data
    header =  drizzle_obj[0].header
  
    #Get orientation of the drizzled image
    drizzle_orientation=header['ORIENTAT']  
  
    Dimensions = ImageData.shape
    
    Orientation = np.zeros( nImages )
 
    InDrizzleFrame= np.zeros(( nImages,  3, nGalaxies) )

    #Now convert into wcs
    
    ra, dec = at.pix2deg( drizzle_file, moments.X_IMAGE, moments.Y_IMAGE)

    newcol = []
    for iImage in xrange(nImages):
        print( "%i/%i" %(iImage+1, nImages))
        #now see where each of our positions lie on each of the individual images
        
        SingleImageX, SingleImageY = at.deg2pix( individual_files[iImage], ra, dec)    
        #Get the limits of the field and find which one are
        #in the FOV and on the chip
        isin, image_orientation = al.acs_limits( SingleImageX, SingleImageY, \
                                                 individual_files[iImage])

        
        Orientation = np.zeros( len(SingleImageX)) +  image_orientation - drizzle_orientation
     
        isinArr= np.zeros( len( moments.X_IMAGE) )
        
        isinArr[ isin ] = 1
     
     

        InDrizzleFrame[iImage, 0, :] = SingleImageX
        InDrizzleFrame[iImage, 1, :] = SingleImageY
        InDrizzleFrame[iImage, 2, :] = isinArr

        iFilename = individual_files[iImage].split('/')[-1]
        x_column = py.Column( name=iFilename+'_X_IMAGE', \
                                format=SingleImageX.dtype, \
                                array=SingleImageX )

        y_column = py.Column( name=iFilename+'_Y_IMAGE', \
                                format=SingleImageY.dtype, \
                                array=SingleImageY )
        
        inFrame = py.Column( name=iFilename+'_INFRAME', \
                                format=isinArr.dtype, \
                                array=isinArr )
                                
        orientat = py.Column( name=iFilename+'_ORIENTAT', \
                                format=Orientation.dtype, \
                                array=Orientation )

        newcol.append( x_column )
        newcol.append( y_column )
        newcol.append( inFrame )
        newcol.append( orientat )
        
    orig_cols = moments.columns
    new_cols = py.ColDefs(newcol)
    
    hdu = py.BinTableHDU.from_columns(orig_cols + new_cols)
    
    num_exposures = np.sum(InDrizzleFrame[:, 2, :], axis=0).shape

    return hdu.data



   

   
