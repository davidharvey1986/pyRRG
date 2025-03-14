from .getHSTfilter import getHSTfilter 
import os
import json
import subprocess 
from . import directories as directories
from . import check_external_packages as cep
from astropy.io import fits

def setDefaultParams( params ):

    if params["data_dir"] is None:
        params["data_dir"] = os.getcwd()+'/'    

    #Check files exist
    params["field"] = os.path.join(params["data_dir"],params["FILENAME"])
    params["weight_file"] = os.path.join(params["data_dir"],params["weight_file"])
    if params["root_name"] is None:
        params["root_name"] = params['FILENAME'].split('.')[0]
        
    params["root_name"] = os.path.join(params["output_dir"],params["root_name"])
    
    if not os.path.isfile( params["field"] ):
        raise ValueError('Cant find input image (%s)' % params["field"])
        
    
    params["hst_filter"] = getHSTfilter(params)
    params["wavelength"] = ''.join([  s for s in params["hst_filter"] if s.isdigit()])
    

    

    
    #get the cwd first, and make sure these are aboslute paths!
    if params["code_dir"] is None:
        params["code_dir"] = '/'.join(os.path.abspath(__file__).split('/')[:-1])
        
    if params["sex_files"] is None:
        params['sex_files']=params['code_dir']+'/sex_files/'

        
    if params['jwst'] and (params["psf_model"] == 'tinytim'):
        raise ValueError("You have selected TinyTim PSF model with JWST -> this is not allowed")
    
    if (not params['jwst']) and (params["psf_model"] != 'tinytim') or (params["psf_model"] == 'empirical'):
        raise ValueError("You have selected HST but not a HST compatible PSF model")

    params["psf_model"] = params["psf_model"].lower( )

    if params["psf_model"] != "empirical":
        if params["psf_model_dir"] is None:
            params["psf_model_dir"]=params['code_dir']+'/psf_lib/%s' % params["psf_model"].lower()
        
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
    cep.check_external_packages()

    
    
    if params['jwst']:
        params['zero_point'] = 'jwst'
    else:
        params['zero_point'] = 'hst'
                   
    
    if params['fits_extension'] is None:
        if params['jwst']:
            params['fits_extension'] = 1
        else:
            params['fits_extension'] = 0
        try:
            test = fits.open( params["field"] )[ params['fits_extension']]
            assert len(test.data.shape) == 2, "Data does not seem to be an image, consider fits extension"
        except:
            raise ValueError("Failed to open fits file (%s) with"
                             "extension (%i) consider --fits_extension keyword" %
                             ( params['FILENAME'],  params['fits_extension']))
            
            
                
    if params['expTimeName'] is None:
        if params['jwst']:
            params['expTimeName'] = 'XPOSURE'
        else:
            params['expTimeName'] = 'EXPTIME'
        
    if params['orientation_header'] is None:
        if params['jwst']:
            params['orientation_header'] = 'PA_V3'
        else:
            params['orientation_header'] = 'ORIENTAT'      
            
    json.dump(params, open("pyRRG.params","w"))
    return params

