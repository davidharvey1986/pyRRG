from .getHSTfilter import getHSTfilter 
import os
import json
import subprocess 
from . import directories as directories
from . import check_external_packages as cep
from astropy.io import fits

def setDefaultParams( params ):
    
    params["hst_filter"] = getHSTfilter(params)
    
    params["wavelength"] = ''.join([  s for s in params["hst_filter"] if s.isdigit()])
    
    #get the cwd first, and make sure these are aboslute paths!
    if params["code_dir"] is None:
        params["code_dir"] = '/'.join(os.path.abspath(__file__).split('/')[:-1])
        
    if params["sex_files"] is None:
        params['sex_files']=params['code_dir']+'/sex_files/'
    
    if params["psf_model_dir"] is None:
        if params['jwst']:
            params["psf_model_dir"]=params['code_dir']+'/psf_lib_jwst/'
        else:
            params["psf_model_dir"]=params['code_dir']+'/psf_lib/'

    if params["data_dir"] is None:
        params["data_dir"] = os.getcwd()+'/'    
   
    
    params['stilts_dir'] = '/'+'/'.join(str(subprocess.check_output(['which','stilts.sh'])).split('/')[1:-1])
    
    params["dirs"] = directories.directories(params['data_dir'],  params['sex_files'],
                           params['psf_model_dir']+'/'+str(params['wavelength'])+'/',
                                       params['code_dir'], params['stilts_dir'])
    params["dirs"].check_dirs()
    params["dirs"].write_dirs()
    cep.check_external_packages()

    
    #Check files exist
    params["field"] = params["dirs"].data_dir+params["FILENAME"]
    
    if not os.path.isfile( params["field"] ):
        raise ValueError('Cant find input image (%s)' % params["field"])
        
    

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
            test = fits.open( params['FILENAME'] )[ params['fits_extension']]
            assert len(test.data.shape) != 2, "Data does not seem to be an image, consider fits extension"
        except:
            raise ValueError("Failed to open fits file extension consider --fits_extension keyword")
            
            
                
    if params['weight_file'] is None:
        if params['jwst']:
            params['weight_file'] = params['FILENAME']+'[4]'
        else:
            params['weight_file'] = params['FILENAME'][:-8]+'wht.fits'
     
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