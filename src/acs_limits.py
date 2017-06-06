import numpy as np
import pyfits as py


def acs_limits( x, y, filename):
    '''
    ;PURPOSE : TO RETURN THE INDEX OF THOSE X AND Y POSITIONS THAT
    ;          ARE WITHIN THE FIELD OF VIEW AND NOT ZERO IN THE SINGLE
    ;          IMAGE FRAM
    
    ;INPUTS :
    ;    X : THE X POSITIONS of the sources in pixel coordinates
    ;    Y : Y POSITIONS of the sources in pixels coordinates
    ;    ASTRO : A STRCUTURE FROM EXTAST INCLUDING FILENAME
    
    ;OUTPUTS : A BOOLEAN OF THE INDEXES OF THE X AND Y POSITIONS THAT ARE IN THE FIELD
    '''
    
    image_obj=py.open(filename)
    imagedata = image_obj[0].data
    header = image_obj[0].header
    
    zeroarr = imagedata == 0

    
    isin = (x > 0) & ( x < header['NAXIS1']) & \
        (y > 0) & ( y < header['NAXIS2'])

    zeroarr = np.arange( len( y[isin] ))[imagedata[ y[isin].astype(int), x[isin].astype(int)] == 0]
    #Create an array the size of the entire field even
    #the positions outside the FOV
    SourceArr=np.ones(len(x))  
    
    #All the positions in question make 1 in the array
    #SourceArr[np.round(y-min(y)).astype(int)-1, np.round(x-min(x)).astype(int)-1] = 1
    SourceArr[ isin ] += 1
    
    isin_index = np.arange(len(y))[isin]
    
    #Deduct one off all the positions that arent in the chip
    SourceArr[ isin_index[zeroarr] ] -= 1
    
    #then all the chip positions that are in the FOV AND CHIp are value 2
   
    return SourceArr == 2, header['ORIENTAT']


     
