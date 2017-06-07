import os as os
import numpy as np
import glob as glob
import RRGtools as at
import pickle as pkl
import measure_moms as measure_moms
import star_galaxy_separation as sgs
import pyfits as py
import calc_shear as cs
import psf_cor as psf
import plot_shears as plot
import ellipse_to_reg as etr
import directories as directories
import rrg_to_lenstool as rtl
import subprocess
import check_external_packages as cep

def main(  infile, hst_filter=None,
            data_dir=None,
            code_dir=None,
            sex_files=None,
            psf_model_dir=None,
            expThresh = 3, 
            noisy=False, 
            nonstop=True, 
            fits_cat=None, 
            mag_cut=[0.,40.], 
            signal_noise_cut=4.4,
            size_cut=[3., 30.],
            min_rad=6.,
            mult=2.):
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
    ;         FILTER : THE HST FILTER USED
    
    ;KEYWORDS : EXPTHRESH : THE MINIMUM NUMBER OF EXPOSURES ONE
    ;                       GALAXY MUST HAVE BEFORE ALLOWED TO
    ;                       BE IN SAMPLE; DEFAULT = 2 TO PREVENT
    ;                       EDGE GALAXIES IN THE DITHER TO BE INCL.
    ;           FILTER2 : THE SECOND FILTER USED TO GET THE RED_SEQUENCE
    ;                     IF NOT SET THEN RED_SEQUENCE IS NOT FOUND
    ;           NOISY : RUN THE SOURCE EXTRACTION USING THE OLD AND
    ;                   NOT THE MATHILDE NEW ONE THAT IS MROE SENSITIVE
    ;           NONSTOP : DO NOT PAUSE FOR CONFIRATION I WANT TO CONTINUE
    ;           FITS_CAT :  MAKE A CATALOGUE THAT IS IN FITS FORMAT
    ;           MAG_CUT : MAGNITUDE CUTS
    

    ;TO DO : 
    ;        2. CHANGE SUCH THAT FUNCTIONS DONT CONSISTENTLY READ
    ;           FITS IMAGES AND SLOW THIGNS DOWN
    
    '''
    if hst_filter is None:
        hst_filter='F814W'
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
   
    stilts_dir = '/'.join(subprocess.check_output(['which','stilts.sh']).split('/')[:-1])
        
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
  
    Exposures = glob.glob( dirs.data_dir+'/j*.fits ')
    nExposures = len(Exposures)
  
    if nExposures  < expThresh:
        expThresh = nExposures
        print 'WARNING: Low number of exposures'
    
  

    # Define survey parameters
    #------------------------------------------
    #Now as keywords


    sex_catalogue = field[:-5]+"_sex.cat"
    
    #Find objects and measure their raw shapes
    if not os.path.isfile( sex_catalogue):
        weight_file = infile[:-8]+'wht.fits'
        sources = at.source_extract( infile, weight_file,
                                         outfile=sex_catalogue,
                                         conf_path=dirs.sex_files,
                                         dataDir=dirs.data_dir)
    else:
        sources = py.open( sex_catalogue )[1].data

     
        
    print sex_catalogue
  
  
    uncorrected_moments_cat = field[:-5]+"_uncor.cat"
    
    if not os.path.isfile(uncorrected_moments_cat):
        measure_moms.measure_moms( infile,
                                   sex_catalogue,
                                   uncorrected_moments_cat,
                                    min_rad=min_rad, mult=mult)

    uncorrected_moments = py.open( uncorrected_moments_cat )[1].data
 

    
    galaxies, stars = sgs.star_galaxy_separation( uncorrected_moments,
                                                  savefile='galStar.locus' )

    
    n_stars=len(stars)
  
    corrected_moments_cat = field[:-5]+"_cor.cat"

    #Correct for the PSF
    if not os.path.isfile(corrected_moments_cat):
         psf.psf_cor( uncorrected_moments_cat,
                    corrected_moments_cat,
                    infile, wavelength,
                    mult=1, min_rad=min_rad, chip=1,
                    constantpsf=0, mscale=0, 
                    num_exposures=1, order=3,
                    n_chip=2)
    

    corrected_moments = py.open( corrected_moments_cat )[1].data

    #Correct zerpoint for the stacked num exposures
  
    sheared_cat = field[:-5]+".shears"
    
    cs.calc_shear( corrected_moments, galaxies,
                   sheared_cat, 
                    min_rad=min_rad, mult=mult,
                    signal_noise_cut=signal_noise_cut,
                    size_cut=size_cut,
                    mag_cut=mag_cut,
                    dataDir=data_dir)
    

    
    plot.plot_shears( sheared_cat )
    etr.ellipse_to_reg( sheared_cat )

    lenstool_file =  field[:-5]+".lenstool"
    rtl.rrg_to_lenstool( sheared_cat, field)
  

