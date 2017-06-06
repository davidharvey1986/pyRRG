import idlsave as idl
import numpy as np
import pyfits as py
import directories




def acs_model_e( x, y, wavelength='814', data_dir=None):
    '''
    ;+
    ; NAME:
    ;       ACS_MODEL_E
    ;
    ; PURPOSE:
    ;       Obtain the ellipticity of a TinyTim model PSF at one
    ;       position and any focus value(s) on ACS.
    ;  
    ; CATEGORY:
    ;       ACS data reduction.
    ;
    ; CALLING SEQUENCE:
    ;       result=acs_model_e(x, y, focus, catalogues=catalogues)
    ;
    ; INPUTS:
    ;       X          - Object x coordinate.
    ;       Y          - Object y coordinate.
    ;
    ; OPTIONAL INPUTS:
    ;       FOCUS      - Focus offset from optimal [um].
    ;                    Can be one number or an array.
    ;       CATALOGUES - Array of model ellipticities. If not set (e.g. on first
    ;                    use), will read in from disk. If already set, will recall
    ;                    from memory, for speed.
    ;       DITHER     - Dither pattern (PSF is averaged from exposures at all of
    ;                    these positions)
    ;
    ; KEYWORD PARAMETERS:
    ;       None.
    ;
    ; OUTPUTS:
    ;       Returns a structure containing model ellipticities and other parameters.
    ;
    ; OPTIONAL OUTPUTS:
    ;       None.
    ;
    ; COMMON BLOCKS:
    ;       None.
    ;
    ; SIDE EFFECTS:
    ;       None.
    ;
    ; RESTRICTIONS:
    ;       None.
    ;
    ; EXAMPLE:
    ;       See use in acs_determine_focus.pro.
    ;
    ; MODIFICATION HISTORY:
    ;       Feb 05 - All shape moments needed by RRG now returned. RM.
    ;       Feb 05 - Written by Richard Massey.
    ;-
    ; 
  '''

    dirs = directories.return_dirs()


    # Initialise variables
  

    filenames1=dirs.psf_model_dir+"/TinyTim_f"
    filenames2=".moms"
    n_focus = 16
    catalogues_focus = np.arange( n_focus ) - 10
    catalogues_string = np.array([ "f"+str(np.int(i)) for i in catalogues_focus])
    n_coords=len(x)
  
    # Load data (required on first use only)
    catalogues={ 'focus':catalogues_focus,
                 'focus_string':catalogues_string}
        
  
    # Prepare empty arrays to contain the answer
    answer = return_moms( n_focus, n_coords, x, y)
    
    # Loop over each focus value in turn
    for f in xrange(n_focus):
        moms = idl.read( filenames1+str(catalogues_focus[f])+filenames2).moms
        
        
        # Loop over each position in which we're interested
        for i in xrange(n_coords):
            # Find closest model star (i.e. don't inerpolate)
            rmin = np.sqrt( (moms.x[0]-x[i])**2 + (moms.y[0]-y[i])**2)
            match = rmin == np.min(rmin)
            
            # Store its ellipticity and shape moments

            answer.e1[f,i]=moms.e1[0][match]
            answer.e2[f,i]=moms.e2[0][match]
            answer.xx[f,i]=moms.xx[0][match]
            answer.xy[f,i]=moms.xy[0][match]
            answer.yy[f,i]=moms.yy[0][match]
            answer.xxxx[f,i]=moms.xxxx[0][match]
            answer.xxxy[f,i]=moms.xxxy[0][match]
            answer.xxyy[f,i]=moms.xxyy[0][match]
            answer.xyyy[f,i]=moms.xyyy[0][match]
            answer.yyyy[f,i]=moms.yyyy[0][match]
      
            #Store its other quantities
            answer.x_model[f,i]=moms.x[0][match]
            answer.y_model[f,i]=moms.y[0][match]
            answer.offset_model[f,i]=np.min(rmin)
  
 

    return answer, catalogues_focus

class return_moms( dict ):

    def __init__( self, n_focus, n_coords, x, y):

        self.__dict__['x'] = x
        self.__dict__['y'] = y
        self.__dict__['e1']=np.zeros((n_focus,n_coords))
        self.__dict__['e2']=np.zeros((n_focus,n_coords))
        self.__dict__['xx']=np.zeros((n_focus,n_coords))
        self.__dict__['xy']=np.zeros((n_focus,n_coords))
        self.__dict__['yy']=np.zeros((n_focus,n_coords))
        self.__dict__['xxxx']=np.zeros((n_focus,n_coords))
        self.__dict__['xxxy']=np.zeros((n_focus,n_coords))
        self.__dict__['xxyy']=np.zeros((n_focus,n_coords))
        self.__dict__['xyyy']=np.zeros((n_focus,n_coords))
        self.__dict__['yyyy']=np.zeros((n_focus,n_coords))
        self.__dict__['x_model']=np.zeros((n_focus,n_coords))
        self.__dict__['y_model']=np.zeros((n_focus,n_coords))
        self.__dict__['offset_model']=np.zeros((n_focus,n_coords))

    def keys(self):
        return self.__dict__.keys()

    def __getitem__(self, key): 
        return self.__dict__[key]

            
