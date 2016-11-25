import os as os
import numpy as np
import glob as glob
import astro_tools as at
import pickle as pkl
import measure_moms as measure_moms
import star_galaxy_separation as sgs
import pyfits as py
import calc_shear as cs
import psf_cor as psf
import ipdb as pdb
import plot_shears as plot
import ellipse_to_reg as etr
def main(  infile, hst_filter=None,
            data_dir='./',
            sex_files=None,
            psf_model_dir=None,
            expThresh = 3, 
            noisy=False, 
            nonstop=True, 
            fits_cat=None, 
            mag_cut=25., 
            signal_noise_cut=4.4,
            size_cut_lo=3.,
            size_cut_hi=30,
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
    if sex_files is None:
        sex_files='~/Library/Code/python/lensing/lenstool/sex_files/'
        
    if  psf_model_dir is None:
        psf_model_dir='~/Library/Code/python/lensing/lenstool/psf_lib/'

        
    dirs = directories(data_dir,  sex_files, psf_model_dir+'/'+str(wavelength)+'/' )

    field=dirs.data_dir+infile

    if not os.path.isfile( field):
        raise ValueError("%s not found" % field)
  
    Exposures = glob.glob( dirs.data_dir+'/j*.fits ')
    nExposures = len(Exposures)
  
    if nExposures  < expThresh:
        expThresh = nExposures
        print 'WARNING: Insufficient exposures, continuing anyway'
  
    #Make this min of three unless number of exposures is only 2.
    
  

    # Define survey parameters
    #------------------------------------------
    #Now as keywords

 
    sex_catalogue = field[:-5]+"_sex.cat"
    print sex_catalogue
    #Find objects and measure their raw shapes
    if not os.path.isfile( sex_catalogue):
        weight_file = infile[:-8]+'wht.fits'
        sources = at.source_extract( infile, weight_file,
                                         outfile=sex_catalogue )
    else:
        sources = py.open( sex_catalogue )[1].data

     
        
    print sex_catalogue
  
  
    uncorrected_moments_cat = field[:-5]+"_uncor.cat"
    
    if not os.path.isfile(uncorrected_moments_cat):
        measure_moms.measure_moms( infile,
                                   sex_catalogue,
                                   uncorrected_moments_cat,
                                    min_rad=min_rad, mult=mult,
                                    skymed=-0.000178,
                                    skysd=0.002472 )

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
                    n_chip=2, dataDir=data_dir)
    

    corrected_moments = py.open( corrected_moments_cat )[1].data

  #Correct zerpoint for the stacked num exposures
  
    sheared_cat = field[:-5]+".shears"
    
    cs.calc_shear( corrected_moments, galaxies, sheared_cat, 
                                        min_rad=min_rad, mult=mult,
                                        signal_noise_cut=signal_noise_cut,
                                        size_cut_lo=size_cut_lo,
                                        size_cut_hi=size_cut_hi,
                                        dataDir=data_dir)
    

    
    plot.plot_shears( sheared_cat )
    etr.ellipse_to_reg( sheared_cat )
    '''
  if keyword_set(filter2) then begin
     red_sequence = source_gal(cluster_name, filter, filter2)
  endif else begin
     red_sequence = findgen(n_elements(momc.x))
  endelse
  redSeqFlag = fltarr( n_elements(momc.x))
  redSeqFlag[red_sequence] += 1
  
  

; Read in catalogue
  restore,field+"_cor.cat"     ;,/verbose
  sn=momc.sn/0.56              ; Correction for correlated noise (see Leauthaud et al.)
  flux=10^(-0.4*momc.mag)





; Form shear estimator
  gamma1=e1/g1/0.86             ; See Rhodes et al.
  gamma2=e2/g1/0.86
;Is this the size of the galaxy?
;If so should it not be
  d=sqrt(((momc.xx+momc.yy)/2)>0) ; size in pixels


; Weight good galaxies
  weight=fltarr(n_elements(momc.x))
  weight[galaxies]=1.
  
  ;cut out the member ( biggest galaxies ) and any galaxy near them.
  ;sizes 
  near_cut = 5.0                  ;Factor bigge to incorporate near galaxies

  indexGalMembers = get_cluster_members( momc, size_cut_upper,$
                                         near_cut=near_cut )
  
 if not keyword_set(mag_cut) then mag_cut=[0,40]
 if n_elements(mag_cut) ne 2 then begin
    print, 'Mag cut is not properly defined'
    return
 endif
  good=where(d gt size_cut $
             and momc.sn gt sn_cut and $
             momc.num_exposures gt expThresh and $
             RedSeqFlag gt 0 and momc.mag gt mag_cut[0] $
             and momc.mag lt mag_cut[1], n_good)
  
  ;if no galaxy members then indexGalMemebrs = -1
  if indexGalMembers[0] gt -1 then begin
     goodTmp = 0
     for iGal=0, n_good-1 do $
        if where( good[iGal] eq indexGalMembers) eq -1 then $
           goodTmp = [ goodTmp, good[iGal] ]
     good = goodTmp[1:*]
  endif
  
  
  


  weight[good]+=1.
  weight[stars] = 0

;Make a stick plot
  plot_stick, momc, gamma1, gamma2, weight=weight

; Write out catalogue in useable form
  openw,lun,/get_lun,field+".rrg"
  printf,lun,"           x           y         mag        size      gamma1      gamma2"
  for i=0,n_elements(momc.x)-1 do $
     if weight[i] gt 1 then $
        printf,lun,format="(6F15.6)",$
               float(momc.x[i]),float(momc.y[i]),momc.mag[i],d[i],gamma1[i],gamma2[i]
                                ;float(momc.ra[i]),float(momc.dec[i]),momc.mag[i],d[i],gamma1[i],gamma2[i]
  close,lun
  free_lun,lun

  ellconverter, cluster_name, filter

  if keyword_set(fits_cat) then begin
     header = headfits( field+'_sci.fits')
     extast, header, astro
     xy2ad, momc.x, momc.y, astro, $
            ra, dec
     good_ones = where( good gt 1, nGood)
     
     names = [ 'RA', 'DEC', 'MAG', 'SIZE', 'GAMMA1', 'GAMMA2']
     types = ['1.d0','1.d0','1.d0','1.d0','1.d0','1.d0']
     cat_str = mrd_struct( names, types, nGood)
     cat_str.RA =ra[good_ones]
     cat_str.DEC = dec[good_ones]
     cat_str.MAG = momc.MAG[good_ones]
     cat_str.SIZE = d[good_ones]
     cat_str.GAMMA1 = GAMMA1[good_ones]
     cat_str.GAMMA2 = GAMMA2[good_ones]
     
     MWRFITS, cat_str, fits_cat, Header
  endif
  
end

'''

class directories( dict ):

    def __init__( self, data_dir, sex_files, psf_mode_dir):
        self.__dict__['sex_files'] = sex_files
        self.__dict__['data_dir'] = data_dir
        self.__dict__['psf_mode_dir'] = psf_mode_dir
        self.write_dirs()

    def write_dirs( self ):
        file_obj = open('directories.cat',"wb")
        file_obj.write("DATA_DIR: %s \n" %self.data_dir)
        file_obj.write("SEX_FILES: %s \n" %self.sex_files)
        file_obj.write("PSF_MODEL_DIR: %s \n" %self.psf_mode_dir)
      
    def keys(self):
        return self.__dict__.keys()

    def __getitem__(self, key): 
        return self.__dict__[key]
