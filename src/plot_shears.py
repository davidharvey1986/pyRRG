from matplotlib import pyplot as plt
import pyfits as py
import numpy as np
from matplotlib import gridspec as gridspec

def plot_shears( moments_catalogue, nbins=None,
                 min_gals_per_bins=100., catalogue=None ):
    '''

    Plot a quiver plot of the 'gamma1' and 'gamma2'


    '''


    fig = plt.figure()
    gs = gridspec.GridSpec( 1, 1)
    ax1 = plt.subplot( gs[0,0])

    if catalogue is None:
        catalogue = py.open( moments_catalogue )[1].data
    nGalaxies = len( catalogue.x)

    
 
    
    if nbins is None:
        nbins = np.int(nGalaxies / min_gals_per_bins)
    xbins = np.linspace( np.min(catalogue.x), np.max(catalogue.x), nbins + 1)
    ybins = np.linspace( np.min(catalogue.y), np.max(catalogue.y), nbins + 1)
    
    e1_map = np.zeros( (nbins, nbins) )
    e2_map = np.zeros( (nbins, nbins) )
    theta_map = np.zeros( (nbins, nbins) )


    xgrid, ygrid = np.meshgrid( xbins[:-1] - np.mean( xbins), ybins[:-1] - np.mean( ybins) )

        
    for i in xrange( nbins ):
        for j in xrange( nbins ):
            in_bin = (catalogue.x > xbins[i] ) & \
              (catalogue.x < xbins[i+1] ) &\
              (catalogue.y > ybins[j] ) & \
              (catalogue.y < ybins[j+1] )

            theta = np.arctan2( np.mean(catalogue.gamma2[ in_bin ]),
                                np.mean(catalogue.gamma1[ in_bin ]) )/2.
            
            
            gamma = np.sqrt( np.mean(catalogue.gamma1[ in_bin ])**2 +
                             np.mean(catalogue.gamma2[ in_bin ])**2 )
            
            e1_map[ j, i] = gamma*np.cos(theta) #np.mean(catalogue.gamma1[ in_bin ])
            e2_map[ j, i] = gamma*np.sin(theta) #np.mean(catalogue.gamma2[ in_bin ])
            theta_map[ j, i] = theta
    
    
    

    quiveropts = dict(color='white', headlength=0, pivot='middle', 
                    linewidth=.5, units='xy', width=.5, headwidth=1) 
    ax1.quiver( xgrid, ygrid, e1_map, e2_map, headaxislength=0, **quiveropts )
    plt.show()

    #bin the x and y
    
