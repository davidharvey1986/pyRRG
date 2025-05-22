'''
return the HST filter of the input image
'''

from astropy.io import fits
import json

def getHSTfilter( params ):
    '''
    Take the input infile and get the filter name
    '''
    
    header = fits.open( params['field'] )[0].header

    instrument =  header['INSTRUME'].upper()

    
    if params['jwst']:
        image_filter = header['FILTER']

    else:
        if 'CLEAR' in header['FILTER1']:
            image_filter = header['FILTER2']
        else:
            image_filter = header['FILTER1']
        
    print(("Using instrument: %s with filetr: %s for image %s" %
           (instrument, image_filter,  params['field'])))
    return instrument, image_filter
