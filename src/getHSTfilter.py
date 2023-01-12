'''
return the HST filter of the input image
'''

from astropy.io import fits
import json

def getHSTfilter( params ):
    '''
    Take the input infile and get the filter name
    '''
    
    header = fits.open( params['FILENAME'] )[0].header

    if params['jwst']:
        print(("Using filter %s for image %s" % (header['FILTER'],  params['FILENAME'])))

        return header['FILTER']

    if 'CLEAR' in header['FILTER1']:
        hst_filter = header['FILTER2']
    else:
        hst_filter = header['FILTER1']
        
    print(("Using filter %s for image %s" % (hst_filter,  params['FILENAME'])))
    return hst_filter
