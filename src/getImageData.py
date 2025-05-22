'''
return the HST filter of the input image
'''

from astropy.io import fits
import json

def getImageData( params ):
    '''
    Take the input infile and get the filter name
    '''
    
    header = fits.open( params['field'] )[0].header
    telescope =  header['TELESCOP'].upper()
    instrument =  header['INSTRUME'].upper()

    
    if telescope == 'JWST':
        image_filter = header['FILTER'].upper()
    else:
        if 'CLEAR' in header['FILTER1']:
            image_filter = header['FILTER2'].upper()
        else:
            image_filter = header['FILTER1'].upper()
        
    print("Found for image %s: " %  params['field'])
    print("TELESCOPE: %s" % telescope )
    print("INSTRUMENT: %s" % instrument)
    print("FILTER: %s" % image_filter )
    
    return telescope, instrument, image_filter
