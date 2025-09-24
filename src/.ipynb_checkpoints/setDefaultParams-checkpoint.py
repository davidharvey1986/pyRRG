from .getImageData import getImageData 
import os
import json
import subprocess 
from . import directories as directories
from . import check_external_packages as cep
from astropy.io import fits

def setDefaultParams( params ):

    if params["data_dir"] is None:
        params["data_dir"] = os.getcwd()+'/'    

    #Make sure the the filename is not the entire path
    if '/' in params["FILENAME"]:
        raise ValueError("Filename seems to contain path,"+
                         " please use only the name of the"+
                         "image and the data_dir argument for the path"
                        )
        
    #Check files exist
    params["field"] = os.path.join(params["data_dir"],params["FILENAME"])
    
    if params["weight_file"] is None:
        
        params["weight_file"] = None
        print("WARNING : NO WEIGHT FILE USED - THIS IS UNADVISED")
    else:
        params["weight_file"] = os.path.join(params["data_dir"],params["weight_file"])

        if not os.path.isfile(params["weight_file"]):
            raise IOError("Can't find weight file: %s" % params["weight_file"])
            
    if params["root_name"] is None:
        #params["root_name"] = params['FILENAME'].split('.')[0]
        params["root_name"], _ = os.path.splitext(params['FILENAME'])
        
    params["root_name"] = os.path.join(params["output_dir"],params["root_name"])
    
    if not os.path.isfile( params["field"] ):
        raise ValueError('Cant find input image (%s)' % params["field"])
        
    telescope, instrument, image_filter = getImageData(params)
    
    params["telescope"] = telescope
    params["filter"] = image_filter
    params["instrument"] = instrument
    params["wavelength"] = ''.join([  s for s in params["filter"] if s.isdigit()])
    

    

    
    #get the cwd first, and make sure these are aboslute paths!
    if params["code_dir"] is None:
        params["code_dir"] = '/'.join(os.path.abspath(__file__).split('/')[:-1])
        
    if params["sex_files"] is None:
        params['sex_files']=params['code_dir']+'/sex_files/'
        
    default_psf_model = {
        'JWST':'webbPSF',
        'HST':'tinytim'
    }
    
    if params["psf_model"] is None:
        params["psf_model"] = default_psf_model[telescope].lower()
        
        print("WARNING: No PSF MODEL put so defaulting to %s (for %s)" %
              ( params["psf_model"], telescope )
             )

        
    if params["telescope"] == "JWST" and (params["psf_model"] == 'tinytim'):
        raise ValueError("You have selected TinyTim PSF model with JWST image -> this is not allowed")
        
    if  params["telescope"] == "HST" and  params["psf_model"] != 'tinytim':
        raise ValueError("You have selected HST but not a HST compatible PSF model")
    

    if params["psf_model"] != "empirical":
        if params["psf_model_dir"] is None:
            

            params["psf_model_dir"]='%s/psf_lib/%s/%s/' % \
                ( params['code_dir'], 
                  params["psf_model"].lower(), 
                  params["instrument"]
                )
        
    print("Using model from %s " % params["psf_model_dir"])

   
    
    params['stilts_dir'] = '/'+'/'.join(str(subprocess.check_output(['which','stilts.sh'])).split('/')[1:-1])
    
    params["dirs"] = directories.directories(
        params['data_dir'],
        params['output_dir'],
        params['sex_files'],
        params['psf_model_dir']+'/'+str(params['wavelength'])+'/',
        params['code_dir'], params['stilts_dir'])
    
    params["dirs"].check_dirs()
    params["dirs"].write_dirs()
    #cep.check_external_packages()

    
    
    params['zero_point'] = params["telescope"].lower()
   
    
    if params['fits_extension'] is None:
        if params["telescope"]  == 'JWST':
            params['fits_extension'] = 1
        else:
            params['fits_extension'] = 0
        try:
            test = fits.open( params["field"] )[ params['fits_extension']]
            assert len(test.data.shape) == 2, "Data does not seem to be an image, consider fits extension"
        except:
            raise ValueError("Failed to open fits file extension consider --fits_extension keyword")
            
            
                
    if params['expTimeName'] is None:
        if params["telescope"]  == 'JWST':
            params['expTimeName'] = 'XPOSURE'
        else:
            params['expTimeName'] = 'EXPTIME'
        
    if params['orientation_header'] is None:
        if params["telescope"]  == 'JWST':
            params['orientation_header'] = 'PA_V3'
        else:
            params['orientation_header'] = 'ORIENTAT'      
     
    if params['no_exposures'] and (params['expThresh'] > 0):
        print("WARNING: With --no-exposures and exposure threshold > 0"+\
                         "there will be no galaxies, so making exposure threshold=0")
        params['expThresh'] = 0

    compatibility_check = params["telescope"] == 'JWST' and (params['focus_array'] is not None)
    assert ~compatibility_check, \
        "JWST set and focus-array given - this is not compatible since JWST does not require a focus position"
    
    if params['focus_array'] is not None:
        assert os.path.isfile(params['focus_array']), f"Can't find user-defined focus array: {params['focus_array']}"
        
    json.dump(params, open("pyRRG.params","w"))
    return params
