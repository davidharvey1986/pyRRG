import pysex as sex
import RRGtools as at
import os as os

def hot_cold_sextract( dirs, infile,
                       weight_file=None,
                       cluster_name, filter, 
                       noisy=False,
                       zeropoint=None, 
                       fitsext=1):

    '''
    ;PURPOSE : TO RUN THE SEXRACTOR ROUTINE OF 
    ;          MATHILDE
    
    ;INPUTS : DIRS : THE DIRECTORY STUCYURE FROM RUNME
    ;         FIELD : NAME OF THE FITS IMAGE
    ;         CLUSTER_NAME : THE NAME OF THE OBJECT
    ;         FILTER : THE HST FILTER BEING USED
    
    ;KEYWORDS : NOISY : IF THE IMAGE IS NOISY THIS METHOD
    ;                MAY ONLY PICK UP NOISE, SO USE MY 
    ;                OLD CONSERVATIVE METHOD
    ;           ZEROPOINT : THE ZPOINT OF THE IMAGE
    ;           FITSEXT : THE FITSET TO USE WHEN LOADING THE
    ;                     CAT
    '''
    #Global sex config args
    sex_config_args = \
          {"PARAMETERS_NAME":dirs.sex_files+"/rrg.param ",
           "FILTER_NAME":dirs.sex_files+"/gauss_5.0_9x9.conv"
           "STARNNW_NAME":dirs.sex_files+"/default.nnw"
           "CATALOG_NAME":catalogue_name,
           "WEIGHT_IMAGE":weight_file,
           "MAG_ZEROPOINT":str(zeropoint) }
    
    
    if weight_file is None:
        weight_file = infile[:-7]+'wht.fits'
        
    if zeropoint is None:
        zeropoint = at.acs_zero_point( infile )
  

    if os.path.isfile( infile ) | os.path.isfile( weight_file ):
        raise ValueError("EITHER SCIENCE OR WEIGHT FILE NOT FOUND")
        
    
    if noisy:
        catalogue_name = infile[:-4]+'_sex.cat'
        sex_config_file = sex_files+'/rrgpy.sex'
        
           
        sex.run( infile,
                 conf_file=sex_config_file,
                 conf_args=sex_config_args )
                   

    else:
        cold_sex_cat = infile[:-4]+'_cold_sex.cat'
        hot_sex_cat = infile[:-4]+'_hot_sex.cat'
        #FIRST RUN COLD
        if not os.path.isfile(cold_sex_cat):
            sex_config_file = dirs.sex_files+"HFF_cold.param"
            sex_config_args["CATALOG_NAME"] = cold_sex_cat
            sex.run( infile,
                    conf_file=sex_config_file,
                    conf_args=sex_config_args )
            
     
     
     
        #SECOND HOT COLD
        if os.path.isfile( hot_sex_cat ):
            sex_config_args["CATALOG_NAME"] = cold_sex_cat
            sex_config_file = dirs.sex_files+"HFF_hot.param"
            sex.run( infile,
                 conf_file=sex_config_file,
                 conf_args=sex_config_args )
     
     
        print 'Matching sex catalogues'
    
        matched_cat = at.run_match( cold_sex_cat, hot_sex_cat)

        #Inject the cold run objects in to the hot run
        hot_sex_file = py.open( hot_sex_cat )
        hot_sex_data = hot_sex_file[
        cold_sex_file = py.open( cold_sex_cat )
        
        hot_sex_file[2].data[ hot_sex_cat[matched_cat['NUMBER_2']-1]] = \
            cold_sex_cat[matched_cat['NUMBER_1']-1]


        

