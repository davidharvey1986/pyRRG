import os as os
import numpy as np
import glob as glob
import RRGtools as at
import pickle as pkl
from . import measure_moms as measure_moms
from . import star_galaxy_separation as sgs
from astropy.io import fits
from . import calc_shear as cs
from . import psf_cor as psf
from . import plot_shears as plot
from . import ellipse_to_reg as etr
from . import directories as directories
from . import rrg_to_lenstool as rtl
import subprocess
from . import masking_star as mask
from . import double_detection_removal as remove_doubles
from .setDefaultParams import setDefaultParams
from .remove_galaxy import remove_galaxy_members
import sys
import json

def main(  ):
    '''
    ;PURPOSE : RUN RRG OVER THE GIVEN CLUSTER AND FILTER, CAN TAKE
    ;          IN MULTIPLE EXPOSURES
    ;          
    ;          THEN FILTERS THE CATALOGUE AND USES ONLY THOSE
    ;          GALAXIES WHICH HAVE EXPTHRESH NUMBER OF EXPOSURES
    ;          OVER THEM AND OUTPUTS .RRG FILE WHICH IS THE CAT
    ;          IN X-Y COORDINATESD
    
    ;          AND THEN RUNS ELLCONVERTER WHICH INVOKES THE MASKS
    ;          IN THE DS9 FILE, MASK.REG AND OUTPUTS A .LENSTOOL
    ;          FILE WHICH IS READY TO BE PUT ITO LENSTOOL
    
    
    ;INPUTS : CLUSTER_NAME : THE NAME OF THE CLUSTER
    
    ;KEYWORDS : EXPTHRESH : THE MINIMUM NUMBER OF EXPOSURES ONE
    ;                       GALAXY MUST HAVE BEFORE ALLOWED TO
    ;                       BE IN SAMPLE; DEFAULT = 2 TO PREVENT
    ;                       EDGE GALAXIES IN THE DITHER TO BE INCL.
    ;           MAG_CUT : MAGNITUDE CUTS
    

    ;TO DO : 
    ;        2. CHANGE SUCH THAT FUNCTIONS DONT CONSISTENTLY READ
    ;           FITS IMAGES AND SLOW THIGNS DOWN
    
    '''    
    
    default_params = json.load(open("pyRRG.params","r"))
    
    
    #SET GLOBAL PARAMETERS TO BE USED FOR ALL
    params = setDefaultParams( default_params )
    
   
    # Define survey parameters
    #------------------------------------------
    #Now as keywords


    sex_catalogue = params['field'][:-5]+"_sex.cat"

        
    #Find objects and measure their raw shapes
    if not os.path.isfile( sex_catalogue):

        sources = at.source_extract( params['FILENAME'], params['weight_file'],
                                         outfile=sex_catalogue,
                                         conf_path=params['dirs'].sex_files,
                                         dataDir=params['dirs'].data_dir,
                                         zero_point=params['zero_point'],
                                         extension=params['fits_extension'])
    else:
        sources = fits.open( sex_catalogue )[1].data

  
    uncorrected_moments_cat = params['field'][:-5]+"_uncor.cat"
    
    if not os.path.isfile(uncorrected_moments_cat):
        measure_moms( params['FILENAME'], sex_catalogue,
                                 uncorrected_moments_cat, **params)

    uncorrected_moments = fits.open( uncorrected_moments_cat )[1].data
 

    
    sgs.star_galaxy_separation( uncorrected_moments, outfile=uncorrected_moments_cat, batch_run=params['batch_run'])
  
    corrected_moments_cat = params['field'][:-5]+"_cor.cat"

    #Correct for the PSF
    if not os.path.isfile(corrected_moments_cat):
         psf.psf_cor( uncorrected_moments_cat,
                    corrected_moments_cat,
                    params['FILENAME'], **params)
    

    corrected_moments = fits.open( corrected_moments_cat )[1].data

    #Correct zerpoint for the stacked num exposures
  
    sheared_cat = params['field'][:-5]+".shears"
    
    cs.calc_shear( corrected_moments, sheared_cat, **params)
    


    beforeDoubles_cat = params['field'][:-5]+"_clean_withDoubles.shears"
    
    mask.main( sheared_cat, uncorrected_moments_cat,
                   outFile=beforeDoubles_cat, **params)


    clean_cat = params['field'][:-5]+"_clean.shears"

    remove_doubles.remove_object(beforeDoubles_cat, \
                    clean_cat, FWHM_to_radius=params['FWHM_to_radius'])
    
    if params['cluster_members_cat'] is not None:
        print("Removing Cluster Members from %s" % params['cluster_members_cat'])
        remove_galaxy_members( params['cluster_members_cat'], clean_cat)
    
    if not params['batch_run']:
        plot.plot_shears( clean_cat )

    etr.ellipse_to_reg( clean_cat )
    
    lenstool_file =  params['field'][:-5]+".lenstool"
    rtl.rrg_to_lenstool( clean_cat, params['field'], params)
  

