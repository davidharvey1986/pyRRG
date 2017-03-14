import numpy as np

def ra_separation( ra1, dec1, ra2, dec2, abs=False):
  '''
  ;PURPOSE : TO DETERMINE THE SEPARATION OF TWO POSITIONS
  ;          ASSUMES SMALL ANGLES


  ;INPUTS : 
  ;    RA1 : VECTOR OR SCALAR OFTHE RIGHT ASCENSION OF 
  ;          THE FIRST POSITION IN DEGREES
  ;    RA2 : VECTOR OR SCALAR OFTHE RIGHT ASCENSION OF 
  ;          THE SECOND POSITION IN DEGREES
  ;    DEC : THE DECLINATION OF THE TWO POSITIONS


  ;KEYWORDS :
  ;    DEC1 : DECLINATION OF THE SECOND HALO, THE DEFAULT
  ;           IS TO HAVE THE HALOS AT THE SAME DECLINATION
  ;    ABS : RETURN THE ABSOLUTE VALUE 

  
  ;RETURNS : 
  ;    SEPARATION : THE ANGULAR SEPARATION OF THE TWO HALOS,
  ;                 POSITIVE IS POSITIVE IN THE SKY (negative east)
  '''
  #convert to radians first
  
  ra1_rad = ra1*np.pi/180.
  ra2_rad = ra2*np.pi/180.

  dec1_rad = dec1*np.pi/180
  dec2_rad = dec2*np.pi/180

  #using the small angle approximation

  separation = np.sqrt( ((ra1_rad-ra2_rad)*np.cos(dec1_rad))**2+
                     (dec1_rad-dec2_rad)**2)*206265.

  
  if abs == False:
    try:
        separation[ ra1 > ra2 ] *= -1
    except:
        if ra1 > ra2:
          separation *= -1


  
   
  return separation

                 
def bin_etang( radial, etang, nbins=20, \
            cut=None, xlim=None, log_bin=False,
            weight=None, **kwargs):
    '''
    PURPOSE : Program to bin up the tangential ellipiticity

    INPUTS :
     - DIST_RAD : The distance the gal is away from the x-axis
     - DIST_LONG : The distance the gal is away from the y-axis
     - E1 : The componet of the shape parallel and perp to the x-axis
     - E2 : The component of the shape 45 deg the x-axis
     - ANGLE : The angle the galaxy is wrt to x-axis
     
    KEYWORDS :
     - NBINS : Number of bins
     - CUT : The distance to cut at  

    '''
    #Get the radial range of the data
    if xlim is None:
      xlim = [np.min(radial), np.max(radial)]

    #Determine the tangential shear
    if weight is None:
        weight = np.ones( len(etang), float)
    
    if log_bin:
        logspace = np.linspace(np.log10(xlim[0]), np.log10(xlim[1]), nbins+1)
        bins = 10**(logspace)
    else:
        bins = np.linspace( xlim[0], xlim[1], nbins+1)

    e_binned = np.zeros((2, nbins), float)
    
    for i in xrange(nbins):

        ind = (radial > bins[i]) & (radial < bins[i+1]) 
        e_binned[0, i] = np.sum( etang[ ind ]*weight[ ind ] )/np.sum( weight[ ind ])
        
        e_binned[1, i] = np.std( etang[ind] )/np.sqrt(len(etang[ind]))

        
    return  (bins[:-1]+bins[1:])/2., e_binned
