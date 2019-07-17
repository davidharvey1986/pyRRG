import numpy as np
import pyfits as py
import os as os

def deg2pix( fits, ra, dec, coordfile=None, postage_stamp=0,
             cut=False, ext=None):
    '''
    Given a fits file, convert from ra and dec in degrees to
    x and y pix

    '''
    os.system('rm -fr pyraf')
    from pyraf import iraf

    os.system( 'mkdir .tmp')
    rand_str = str(np.round(np.random.rand()*1000).astype(int))
    random_name = ".tmp/sky2xy"+rand_str+".fits"
    
    if ext is not None:
        data = py.open( fits )[ext[0],ext[1]]
        py.writeto( random_name, data.data, \
                    header=data.header, clobber=True)
    else:
        os.system("cp "+fits+" "+random_name)
        
    
    

    
    if coordfile is None:
        os.system("rm -fr sky2xy.par")
        coordfile = "sky2xy.par"
       
        np.savetxt( coordfile, np.array((ra,dec)).T, fmt="%0.7f")
    
    iraf.wcsctran.unlearn()
    
    iraf.wcsctran( coordfile, ".tmp/sky2xy.results", random_name, "world", "physical",verbose=True)
    x_chip, y_chip = np.loadtxt( ".tmp/sky2xy.results", unpack=True)
    os.system("rm -fr .tmp")
    if cut:
        try:
            image_size = py.open( fits.split('[')[0] )[1].data.shape
        except:
            image_size = py.open( fits.split('[')[0] )[0].data.shape
            
        
        
        inchip = (y_chip < image_size[1]-postage_stamp) & (y_chip > postage_stamp) & \
        (x_chip > postage_stamp) & (x_chip < image_size[0]-postage_stamp)
    
        return x_chip[inchip], y_chip[inchip]
    else:
        return x_chip, y_chip

def pix2deg( fits, x_image, y_image, coordfile=None, ):
    '''
    Given a fits file, convert from ra and dec in degrees to
    x and y pix

    '''
    os.system('rm -fr pyraf')
    from pyraf import iraf

    os.system("cp "+fits+" xy2sky.fits")
    os.system("rm -fr xy2sky.results")
    if coordfile is None:
        coordfile = "xy2sky.par"
        skypar = open(coordfile,"w") 
        for i in range(len(x_image)):
            skypar.write('%0.5f %0.5f\n' % \
                             (x_image[i],y_image[i]))
        skypar.close()

    iraf.wcsctran.unlearn()
    iraf.wcsctran( coordfile, "xy2sky.results", "xy2sky.fits",  "physical", "world")
    ra, dec = np.loadtxt( "xy2sky.results", unpack=True)
    
  
    return ra, dec
    
def deg2pix_flt( fits, ra, dec, postage_stamp=0, cut=False):
    '''
    Run the deg2pix function but for 2 separate chips
    concatenate them and remove any objects outside the chip

    ostage_stamp allows boardering of
    '''
    
    #chip1
    x_chip1, y_chip1 =  deg2pix( fits, ra, dec, ext=['sci',1])
    
    #chip2
    x_chip2, y_chip2 =  deg2pix( fits, ra, dec, ext=['sci',2])
    
    
    #Check that the objects lie within the chip and
    #which chip they are in
    inchip1 = (y_chip1 < 2048-postage_stamp) & (y_chip1 > postage_stamp) & \
        (x_chip1 > postage_stamp) & (x_chip1 < 4096-postage_stamp)
    inchip2 = (y_chip2 < 2048-postage_stamp) & (y_chip2 > postage_stamp) & \
        (x_chip2 > postage_stamp) & (x_chip2 < 4096-postage_stamp)
            
    x = np.append( x_chip1[ inchip1 ], x_chip2[ inchip2 ])
    y = np.append( y_chip1[ inchip1 ], y_chip2[ inchip2 ]+2048)
    print(x, y)
    return x, y

 



def pix2deg_flt( fits, x, y):
    '''
    Run the deg2pix function but for 2 separate chips
    concatenate them and remove any objects outside the chip
    '''

    #chip1
    ra_chip1, dec_chip1 =  pix2deg( fits+'[sci,1]', \
                                    x[ y < 2048 ], \
                                    y[ y < 2048 ])

    #chip2
    ra_chip2, dec_chip2 =  pix2deg( fits+'[sci,2]', \
                                    x[ y > 2048], \
                                    y[ y > 2048]-2048., \
                                    coordfile="xy2sky.par")


    #Check that the objects lie within the chip and
    #which chip they are in
            
    ra = np.append(ra_chip1, ra_chip2)
    dec = np.append( dec_chip1, dec_chip2)

    return ra, dec

def hmstodd( ra, dec ):
    '''
    Convert from hms to dd
    INPUTS : RA AND DEC ARE STRINGS IN THE FORMAT
    H:M:S SEPARATED BY COLONS
    '''
    ra_float = np.array(ra.split(':')).astype(float)
    dec_float = np.array(dec.split(':')).astype(float)

    
    ra_deg = (ra_float[0] + ra_float[1]/60. + ra_float[2]/3600.)/24.*360.
    if '-' in dec.split(':')[0]:
        hem = -1.
    else:
        hem = 1.
    
    dec_deg = dec_float[0]+hem*dec_float[1]/60. + hem*dec_float[2]/3600.
     
    return ra_deg, dec_deg
    
