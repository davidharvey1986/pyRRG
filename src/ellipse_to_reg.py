import numpy as np
import pyfits as py



def ellipse_to_reg( catalogue_name ):
    '''
    Creat a region file to check thta i have measured ellipses correctedly
    '''


    catalogue = py.open( catalogue_name )[1].data


    regionFile = open( 'ellipse.reg', 'wb')

    ellipticity = np.sqrt( catalogue.e1**2 +
                           catalogue.e2**2 )
    
    semi_major = (1 + ellipticity)*catalogue.gal_size
    semi_minor = (1 - ellipticity)*catalogue.gal_size
    angle = np.arctan2( catalogue.e2, catalogue.e1 )*180./np.pi/2.


    regionFile.write('# Region file format: DS9 version 4.1\n')
    regionFile.write('# Filename: dummy.fits\n')
    regionFile.write("global color=green dashlist=8 3 width=1 font='helvetica 10 normal roman' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n")
    regionFile.write("physical\n")

    for iGal in xrange(len(catalogue.x)):
        
        regionFile.write('Ellipse(%0.4f,%0.4f,%0.4f,%0.4f,%0.4f)\n'
                    % (catalogue.x[iGal], catalogue.y[iGal],
                        semi_major[iGal], semi_minor[iGal], angle[iGal]))
    
