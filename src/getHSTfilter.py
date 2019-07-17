'''
return the HST filter of the input image
'''
import pyfits as fits

def getHSTfilter( infile ):
    '''
    Take the input infile and get the filter name
    '''

    header = fits.open( infile )[0].header

    if 'CLEAR' in header['FILTER1']:
        hst_filter = header['FILTER2']
    else:
        hst_filter = header['FILTER1']
        
    print(("Using filter %s for image %s" % (hst_filter, infile)))
    return hst_filter
