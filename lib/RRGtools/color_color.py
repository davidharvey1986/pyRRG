import numpy as np
from matplotlib import pyplot as plt
import pyfits as py
import nearest_neighbour as nearest

from scipy.stats import gaussian_kde
import match_cat as mc

def magnitude( flux, zpt, exptime, apcor ):
    
    return zpt - 2.5*np.log10(flux/exptime) - apcor

def correct_zpt( zpt, airmass, extinct, default_extinct):

    return zpt - (np.sec(airmass)-1.)*extinct + default_extinct - extinct 

def color( filename_A, filename_B ):
    '''
    Get the color from two cats

    color = mag(A) - mag(B) 

    '''
    data_dir = '/Users/DavidHarvey/Documents/Work/ATLAS/data/masters_student/data'

    matched_cat = mc.match_cat(  data_dir+'/'+filename_A, \
                              data_dir+'/'+filename_B )

    color_cat = np.rec.fromarrays( (matched_cat[1].data['RA_1'], \
                           matched_cat[1].data['DEC_1'], \
                           matched_cat[1].data['MAG_1']-\
                           matched_cat[1].data['MAG_2']), \
                          dtype=[('RA', float), ('DEC', float), ('COLOR', float)] )
    return  color_cat

def match_colors( cat_A, cat_B, cleanup=True, \
               stilts_path='/Users/DavidHarvey/Library/Stilts', \
               search_rad=0.5 ):
    '''
    match two recarrays together
    '''
    HDU_cat_A = py.BinTableHDU( cat_A )
    HDU_cat_B = py.BinTableHDU( cat_B )

    HDU_cat_A.writeto('color_A.fits', clobber=True)
    HDU_cat_B.writeto('color_B.fits', clobber=True)

    command_str =stilts_path+'/stilts.sh tmatch2 in1="color_A.fits" in2="color_B.fits" matcher=sky values1="RA DEC" values2="RA DEC" params="'+str(search_rad)+'" out=matched_color.fits'
    
    mc.os.system(command_str)
    matched_colors = py.open('matched_color.fits')
    if cleanup:
        cleanup_str = 'rm -fr color_A.fits color_B.fits matched_color.fits'
        mc.os.system(cleanup_str)
        
    return matched_colors
    

def color_color( cat_A, cat_B, cat_C, **kwargs ):

    '''

    plot a color color plot of B-A against A-C

    As long as i match both catalogues and they have
    neither one has an -1 then it is okay
    
    '''

    # B-A
    cat_BA = color( cat_B, cat_A )

    # A-C
    cat_AC = color( cat_A, cat_C )

    color_color_cat = match_colors( cat_BA, cat_AC )
    return color_color_cat

def main( **kwargs ):

    ex = ['2', '3']

    GR = np.array([])
    RI = np.array([])
    
    for i in ex:
        iGR_iRI =  color_color('ex'+i+'_r.fits', 'ex'+i+'_g.fits', 'ex'+i+'_i.fits')
        GR = np.append( GR, iGR_iRI[1].data['COLOR_1'])
        RI = np.append( RI, iGR_iRI[1].data['COLOR_2'])

        if i == '2':
            plt.plot( iGR_iRI[1].data['COLOR_1'], iGR_iRI[1].data['COLOR_2'], 'r,')
        if i =='3':
            plt.plot( iGR_iRI[1].data['COLOR_1'], iGR_iRI[1].data['COLOR_2'], 'g,')
    if 'xlabel' in kwargs:
        plt.xlabel( kwargs['xlabel'])
    if 'ylabel' in kwargs:
        plt.ylabel( kwargs['ylabel'])
    if 'filename' in kwargs:
        plt.savefig( kwargs['filename'], format='pdf')
        plt.draw()
    else:
        plt.show()
    
