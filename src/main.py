import os as os
import numpy as np
import glob as glob
import RRGtools as at
import pickle as pkl
from . import measure_moms as measure_moms
from . import star_galaxy_separation as sgs
import pyfits as py
from . import calc_shear as cs
from . import psf_cor as psf
from . import plot_shears as plot
from . import ellipse_to_reg as etr
from . import directories as directories
from . import rrg_to_lenstool as rtl
import subprocess
from . import check_external_packages as cep
from . import masking_star as mask
from . import double_detection_removal as remove_doubles
import sys
from .getHSTfilter import getHSTfilter 

def main(  infile,
            data_dir=None,
            code_dir=None,
            sex_files=None,
            psf_model_dir=None,
            expThresh = 2, 
            mag_cut=[0.,40.], 
            signal_noise_cut=4.4,
            size_cut=[3., 30.],
            min_rad=6.,
            mult=2.,\
            weight_file=None):
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

    hst_filter=getHSTfilter(infile)
    wavelength=''.join([  s for s in hst_filter if s.isdigit()])
                   
    #SET GLOBAL PARAMETERS TO BE USED FOR ALL
    #get the cwd first, and make sure these are aboslute paths!
    if code_dir is None:
        code_dir = '/'.join(os.path.abspath(__file__).split('/')[:-1])
        
    if sex_files is None:
        sex_files=code_dir+'/sex_files/'
    
    if  psf_model_dir is None:
        psf_model_dir=code_dir+'/psf_lib/'
   
    if data_dir is None:
        data_dir = os.getcwd()+'/'    
   
    
    stilts_dir = '/'+'/'.join(str(subprocess.check_output(['which','stilts.sh'])).split('/')[1:-1])
    dirs = directories.directories(data_dir,  sex_files,
                           psf_model_dir+'/'+str(wavelength)+'/',
                                       code_dir, stilts_dir)
    dirs.check_dirs()
    dirs.write_dirs()
    cep.check_external_packages()

    
    #Check files exist
    field=dirs.data_dir+infile
    if not os.path.isfile( field ):
        raise ValueError('Cant find input image (%s)' % field)

    if not os.path.isfile( field):
        raise ValueError("%s not found" % field)
  

    # Define survey parameters
    #------------------------------------------
    #Now as keywords


    sex_catalogue = field[:-5]+"_sex.cat"
    
    #Find objects and measure their raw shapes
    if not os.path.isfile( sex_catalogue):
        if weight_file is None:
            weight_file = infile[:-8]+'wht.fits'
        sources = at.source_extract( infile, weight_file,
                                         outfile=sex_catalogue,
                                         conf_path=dirs.sex_files,
                                         dataDir=dirs.data_dir)
    else:
        sources = py.open( sex_catalogue )[1].data

  
    uncorrected_moments_cat = field[:-5]+"_uncor.cat"
    
    if not os.path.isfile(uncorrected_moments_cat):
        measure_moms( infile,
                                   sex_catalogue,
                                   uncorrected_moments_cat,
                                    min_rad=min_rad, mult=mult,
                                       silent=True)

    uncorrected_moments = py.open( uncorrected_moments_cat )[1].data
 

    
    sgs.star_galaxy_separation( uncorrected_moments, outfile=uncorrected_moments_cat)
  
    corrected_moments_cat = field[:-5]+"_cor.cat"

    #Correct for the PSF
    if not os.path.isfile(corrected_moments_cat):
         psf.psf_cor( uncorrected_moments_cat,
                    corrected_moments_cat,
                    infile, wavelength,
                    mult=1, min_rad=min_rad, chip=1,
                    constantpsf=0, mscale=0, 
                    order=3,
                    n_chip=2)
    

    corrected_moments = py.open( corrected_moments_cat )[1].data

    #Correct zerpoint for the stacked num exposures
  
    sheared_cat = field[:-5]+".shears"
    
    cs.calc_shear( corrected_moments,
                   sheared_cat, 
                    min_rad=min_rad, mult=mult,
                    signal_noise_cut=signal_noise_cut,
                    size_cut=size_cut,
                    mag_cut=mag_cut,
                    dataDir=data_dir,\
                       expThresh=expThresh)
    


    beforeDoubles_cat = field[:-5]+"_clean_withDoubles.shears"
    mask.main( sheared_cat, corrected_moments_cat,
                   outFile=beforeDoubles_cat)


    clean_cat = field[:-5]+"_clean.shears"

    remove_doubles.remove_object(sheared_cat, clean_cat, FWHM_to_radius=1)
    
    plot.plot_shears( clean_cat )

    etr.ellipse_to_reg( clean_cat )
    
    lenstool_file =  field[:-5]+".lenstool"
    rtl.rrg_to_lenstool( clean_cat, field)
  

if __name__ == "__main__":
    print(len(sys.argv))
    if len(sys.argv) == 1:
        print("Please to add image name")
        print("python main.py <image_name>")
        raise ValueError()
    
    main(sys.argv[1])
