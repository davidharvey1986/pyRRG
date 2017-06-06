'''
THIS FILE CONTAINS VARIOUS LITTLE ROUTINES TO HELP PLAY
WITH SEXTRACTOR CATALOGUES

1. SEX_TO_REG : CONVERTS A SEXTRACTOR CATALOGUE TO A REGION FILE
                ( sex_cat, outfile,
                coordinate_sys='IMAGE'
                SexStructure=None,
                cat_type='FITS',
                param_file=None):

        
2. CHECK_FIEKDS : ENSURE THAT THE CATALOGUE CONTAINS THE
                  CORRECT FIELDS

3. MASK_SEXCAT   : MASKS OUT A SEXTRACTOR CATALOGYE WITH INPUT MASK_FILE
                 ( sex_cat, mask_file, param_file=None,
                 cat_type='FITS', catalogue=None,
                 outfile=None):
        

4. WRITE_CATALOGYE : WRITE OUT THE CATALOGUE
                      ( catalogue_REC_ARRAY, outfile, cat_type='FITS')
                      
5. ASCII_TO_REC : TAKE AN INPUT ASCII SEXTRACTOR FILE WITH PARAM FILE
                   AND OUTPUT A REC ARRAY
                   ( sex_cat, param_file):

6. CHECK WHETHER A LIST OF POINTS ARE WITHIN VARIOUS REGIONS:
check_rotbox( x_point, y_point,
                   x_ell, y_ell, a_ell,
                   b_ell, t_ell):
check_ellipse( x_point, y_point,
                   x_ell, y_ell, a_ell,
                   b_ell, t_ell):
'''
import sys
import csv as c
import numpy as np
import astro_tools as at
import pyfits as py


def sex_to_reg( sex_cat, outfile,
                coordinate_sys='IMAGE',
                SexStructure=None,
                cat_type='FITS',
                param_file=None):

    '''
    PURPOSE : TO CONVERT A SEXTRACTOR STRUCTURE FROM
                     PYSEX TO DS9 REGION FILE


    ARGUMENT : SEX_CAT : THE NAME OF THE SEX CATALOGUE
               OUTFILE : NAME OF THE OUTPUT REGION FILES

    OPTIONAL KEWYWORSD
               COORDINATE_SYS : DEFAULT IMAGE, USE IMAGE, NOT WCS
               SEXSTRUCTURE : THE STRUCTURE RETURNED FROM PYSEX
                              MUST HAVE COORDINATES AND ELLIPTICITY
               CAT_TYPE : CAN BE FITS OR ASCII
               PARAM_FILE : IF ASCII MUST PROVIDE THE PARAM FILE USED
                 TO CREATE THE SEX CATALOGUE

    '''
    if (cat_type == 'ASCII' ) & (param_file== None):
        print 'IF FILE TYPE IS ASCII I MUST PROVIDE PARAMFILE'
        return 0

    if SexStructure is None:
        if cat_type == 'FITS':
            try:
                cat_list = py.open( sex_cat )
            except:
                print 'Error raised, could not be a fits file, try ascii'
                return 0
            SexStructure = cat_list[-1].data
        else:
            SexStructure = at.ascii_to_rec( sex_cat, param_file)
  
    
    
    x, y, a, b, theta, coordinate_sys = check_fields( SexStructure, \
                                                      coordinate_sys=coordinate_sys)

    regionFile = c.writer( open( outfile, 'wb'))

    #Write the necessary header
    regionFile.writerow(['# Region file format: DS9 version 4.1'])
    regionFile.writerow(['# Filename: dummy.fits'])
    regionFile.writerow(["global color=green dashlist=8 3 width=1 font='helvetica 10 normal roman' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1"])
    regionFile.writerow([coordinate_sys])
    
    for iGal in xrange(len(x)):
        galInfo = ['ellipse('+str(x[iGal]), \
                   str(y[iGal]), \
                   str(a[iGal]), \
                   str(b[iGal]), \
                   str(theta[iGal])+')']
                   
        regionFile.writerow(galInfo)

    

    


def check_fields( structure, coordinate_sys='IMAGE' ):

    '''
    PURPOSE : TO CHECK THE FIELDS HAVE CONSISTENT X, Y, A, B,
    THETA, IN IMAGE OR WORLD COORDINATES

    INPUT : STRUCTURE TO CHECK

    OUTPUT : X: THE X COORDINATES
             Y: THE Y COORDIANTES
             A: SEMI MAJOR AXI
             B: SEMI MINOR AXIS
             THETA : THE ANGLE IN RADIANS
             COORDINATE_SYS: THE COODINATED SYSTEM BEING USED
                            EITHER 'WORLD' OR 'IMAGE'


    '''
    
    try:
        if coordinate_sys == 'IMAGE':
            x = structure['X_IMAGE']
            y = structure['Y_IMAGE']
            a = structure['A_IMAGE']
            b = structure['B_IMAGE']
            theta = structure['THETA_IMAGE']
            coordinate='IMAGE'
        else:
            x = structure['X_WORLD']
            y = structure['Y_WORLD']
            a = structure['A_WORLD']
            b = structure['B_WORLD']
            theta = structure['THETA_WORLD']
            coordinate='WCS'
    except:
        try: 
            x = structure['X_WORLD']
            y = structure['Y_WORLD']
            a = structure['A_WORLD']
            b = structure['B_WORLD']
            theta = structure['THETA_WORLD']
            coordinate='WCS'
        except:
            print 'No coordinates of the right sort found'
            sys.exit(0)


    return x, y, a, b, theta, coordinate
        
def wcs_separation( ra1, dec1, ra2, dec2, mag=True):
    '''
    WORK OUT THE SEPARATION OF TWO POINTS IN THE SKY
    USING SPHEICAL TRIG

    INPUTS : THE RA AND DEC OF EACH POSITION IN RADIANS

    ASSUMES SMALL SEPARATION

    RETURNS : SEPARATION IN ARCSECS
    '''

    
    separation =  np.sqrt( ((ra1-ra2)*np.cos((dec1+dec2)/2.))**2+
                            (dec1-dec2)**2)*206265.


    if not mag:
        separation[ ra1 - ra2 > 0 ] *= -1 

    return separation


    params = np.loadtxt( param_file, dtype=[('params', object)], unpack=True)
    dtypes = [ ( iParam[0], float) for iParam in params ]
   
    return np.loadtxt( sex_cat, dtype=dtypes )




def mask_sexcat( sex_cat, mask_file, param_file=None,
                 cat_type='FITS', catalogue=None,
                 outfile=None):
    '''
    A script that will mask out the SExtractor catalogue
    Inputs:
       sex_cat : Name of the sex_cat; default the fits version
                 but can take ascii
       mask_file : name of the  mask file in physical coords
                   can take circle, ellipse, or square
       param_file : name of the param file if the catalogye is
                     ascii
       cat_type : either FITS or ASCII
       catalogue : catalogyue in the form of a rec_array

    outputs : catalogue that has removed all points lying within
              the regions

    '''

    if catalogue is None:
        if cat_type == 'ASCII' and param_file is None:
            print 'IF CATALOGUE IS ASCII PLEASE GIVE PARAM FILE'
            return
    

    
        if cat_type == 'ASCII':
            catalogue = ascii_to_rec( sex_cat, param_file )
        else:
            cat_list = py.open( sex_cat )
            catalogue = cat_list[-1].data

            
    mask_dtype = [('shape',object), ('pars', object)]
    masked_regions = np.loadtxt( mask_file, delimiter='(',
                                 dtype=mask_dtype )
    
    masked = np.ones( len( catalogue['X_IMAGE']))
    for iMask in xrange(len(masked_regions['shape'])):
        
        
        
        if masked_regions['shape'][iMask] == 'circle':
            x_mask, y_mask, r_mask = masked_regions['pars'][iMask].split(')')[0].split(',')
                
            flag = check_ellipse( catalogue['X_IMAGE'], catalogue['Y_IMAGE'],
                                  np.float(x_mask), np.float(y_mask),
                                  np.float(r_mask), np.float(r_mask), \
                                  0. )

        if masked_regions['shape'][iMask] == 'ellipse':
            x_mask, y_mask, a_mask, b_mask, t_mask = \
              masked_regions['pars'][iMask].split(')')[0].split(',')
                
            flag = check_ellipse( catalogue['X_IMAGE'], catalogue['Y_IMAGE'],
                                  np.float(x_mask), np.float(y_mask),
                                  np.float(a_mask), np.float(b_mask), \
                                  np.float(t_mask)*np.pi/180.)

        if masked_regions['shape'][iMask] == 'rotbox':
            x_mask, y_mask, a_mask, b_mask, t_mask = \
              masked_regions['pars'][iMask].split(')')[0].split(',')
            flag = check_rotbox( catalogue['X_IMAGE'], catalogue['Y_IMAGE'],
                                  np.float(x_mask), np.float(y_mask),
                                  np.float(a_mask), np.float(b_mask), \
                                  np.float(t_mask)*np.pi/180.)
                                  
        masked[ flag ] = 0

    flag = masked == 1

    if outfile is not None:
        if cat_type=='FITS':
            catlist[-1] = catalogue[flag]
            write_catalogue( cat_list, outfile, cat_type=cat_type )
        else:
            write_catalogue( catalogue[flag], outfile, cat_type=cat_type )
    
    return catalogue[flag]


def write_catalogue( catalogue, outfile, cat_type='FITS'):
    '''
    Write out the catalogue into ASCII or FITS format
    '''

    if cat_type == 'FITS':
        catalogue.writeto( outfile, clobber=True )

    else:
        if cat_type == 'ASCII':
            np.savetxt( outfile, catalogue, fmt='%.5f' )
        else:
            print 'UKNOWN FILE FORMAT'
            

        
    

def check_ellipse( x_point, y_point,
                   x_ell, y_ell, a_ell,
                   b_ell, t_ell):
    '''
    return a flag if the point is inside the ellispe region
    '''
    
    rotate_x = (x_point-x_ell)*np.cos(t_ell) + (y_point-y_ell)*np.sin(t_ell)
    rotate_y = (x_point-x_ell)*np.sin(t_ell) - (y_point-y_ell)*np.cos(t_ell)
    
    
    
    return rotate_x**2/a_ell**2 + rotate_y**2/b_ell**2 < 1

        
def check_rotbox( x_point, y_point,
                   x_ell, y_ell, a_ell,
                   b_ell, t_ell):
    '''
    return a flag if the point is inside the rotated box region
    '''
    
    rotate_x = (x_point-x_ell)*np.cos(t_ell) + (y_point-y_ell)*np.sin(t_ell)
    rotate_y = (x_point-x_ell)*np.sin(t_ell) - (y_point-y_ell)*np.cos(t_ell)

    return (np.abs(rotate_x) < a_ell/2. ) & (np.abs(rotate_y) < b_ell/2. )
    
    
def ascii_to_rec( sex_cat, param_file):
    '''
    Convert a ascii catalogue to a recarray
    '''

    params = np.loadtxt( param_file, dtype=[('params', object)], unpack=True)
    dtypes = [ ( iParam[0], float) for iParam in params ]
   
    return np.loadtxt( sex_cat, dtype=dtypes )


def fits_to_reg( fitsfile, outfile, circle_rad=3.):
    '''
    Take a fits file and turn it into a region file
    assuming that is has RA and DEC as pars

    also assuming all sources are circles and have a radius of 3 arcsecs

    '''


    regionFile = open( outfile, 'wb')

    #Write the necessary header
    regionFile.write('# Region file format: DS9 version 4.1\n')
    regionFile.write('# Filename: dummy.fits\n')
    regionFile.write("global color=green dashlist=8 3 width=1 font='helvetica 10 normal roman' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n")
    regionFile.write("fk5\n")

    data = py.open(fitsfile)[-1].data
    ra = data['RA']
    dec = data['DEC']
    for iGal in xrange(len(ra)):
        regionFile.write('circle(%0.4f,%0.4f,%0.1f")\n' % (ra[iGal],dec[iGal],circle_rad))    
