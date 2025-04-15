from astropy.io import fits
import warnings
warnings.filterwarnings('ignore', category=UserWarning, append=True)

def write_rrgparams_to_header( fitsname, params, outputname=None):
    '''
    Write all the rrg params to the header of the catalogue for
    future knowledge
    '''

    hdu = fits.open(fitsname)


    for iparam in params.keys():
        if len(str(params[iparam])) > 8:
            continue
        hdu[1].header['RRG-%s' % iparam] = \
            "%s" % str(params[iparam])

    if outputname is None:
        hdu.writeto( fitsname, overwrite=True)
    else:
        hdu.writeto( outputname )
    
