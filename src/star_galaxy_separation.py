from matplotlib import pyplot as plt
import os as os
import pickle as pkl
from matplotlib import pyplot as plt
from matplotlib import gridspec as gridspec
import sys as sys
import numpy as np
import ipdb as pdb

def star_galaxy_separation( sources, restore=False,
                            savefile='galStar.sav',
                            NonStop=False,
                            axes=None):

    '''
    ;PURPORSE : To separate out the stars and galaxies in the input
    ;           catalogue
    
    ;INPUTS : sources : a sextractor catalogue of sources
    
    ;KEYWORDS : 
    ; MMM : If selected use the moments from
    ;                 rrg_measure_moms.pro structure
    ; stars : return the indices of the stars
    ; restore : restore a saved file of classes
    ; savefile : save the file to the given string
    
    ;RETURNS : A VECTOR OF THE GALAXIES
    
    ; UPDATES : I Shouldnt be working on the positions of
    ;           or the separation using the single images
    ;           I should be using the stacked and therefore dont
    ;           need the single keyword anymore
    
    ;           Use the moments from Measyrubg moments
    ;             and the radius size for stars use mu max for gals
    '''

           
    radius = np.sqrt(sources['ISOAREA_IMAGE']/np.pi)
    
    #Create a box that separateds the galaxies fom teh stars
    #Each box is define by
    #The line separating Stars and Galaxies
    
  
    nObjects = len(sources['X_IMAGE'])
    object_indexes = np.arange( nObjects )
  


    print savefile
    if not os.path.isfile( savefile ):
        gal_star = galStar( savefile, sources, set_defaults=True)
        return object_indexes[ gal_star.galaxies], object_indexes[ gal_star.stars ]
    else:
        if restore:
            gal_star = galStar( savefile, sources, set_defaults=False)   
            gal_star.get_galaxies( sources )
            gal_star.get_stars( sources )
            return object_indexes[ gal_star.galaxies], object_indexes[ gal_star.stars ]
        
        print 'Loading '+savefile
        gal_star = galStar( savefile, sources, set_defaults=False)   
        gal_star.get_galaxies( sources )
        gal_star.get_stars( sources )
        gal_star.generate_axes( sources )
        gal_star.plot_boundaries(  sources )
        gal_star.plot_stars_galaxies( sources, gal_star.axes, overwrite=True )
        plt.show()

        if NonStop:
            overwrite='c'
        else:
            overwrite = \
                raw_input('Save file exists, overwrite Yes (y) Continue (c) No (n)?')
                  
        
        if overwrite == 'y':
            print 'Writing over gal file and removing j*uncor.cat'
            os.system('rm -fr j*uncor.cat *_cor.cat')
            gal_star.write( savefile )
            return object_indexes[ gal_star.galaxies], object_indexes[ gal_star.stars ]
        else:
            if overwrite == 'c':
                return object_indexes[ gal_star.galaxies], object_indexes[ gal_star.stars ]
            else:
                print 'Reseparating'
                os.system('rm -fr '+savefile)
                galaxies, stars = star_galaxy_separation( sources )
                return galaxies, stars

    


class galStar():

        def __init__( self, filename, sources, set_defaults=True ):
            self.star_points = []
            self.gal_points = []
            
            if set_defaults:
                self.generate_axes( sources)
                self.defaults( filename, sources )
            else:
                self.load( filename )


        def defaults( self, filename, sources):
            #Use default parameters (first guess)
            #galaxy locuss parameters
            self.StarsLowCut = 0.
            self.StarsUpCut = 100
            self.stars = []
            self.GalLowCut = 0

            self.get_params_interactively( sources ) 

            self.write( filename )
            
        def write(self, filename):
            galstar_file = open( filename, 'wb')
            galstar_file.write('#Galaxy Parameters\n')
            galstar_file.write('GalLowCut %5.2f\n' % self.GalLowCut)
            galstar_file.write('GradGal %5.2f\n' % self.GradGal)
            galstar_file.write('IntGal %5.2f\n' % self.IntGal)
            galstar_file.write('#Star Parameters\n')
            galstar_file.write('GradStarsLowCut %5.2f\n' % self.GradStarsLowCut)
            galstar_file.write('IntStarsLowCut %5.2f\n' % self.IntStarsLowCut)
            galstar_file.write('StarsUpCut %5.2f\n' % self.StarsUpCut)
            galstar_file.write('StarsLowCut %5.2f\n' % self.StarsLowCut)

        def load( self, filename ):
            name, values = np.loadtxt( filename, unpack=True, dtype=[('name', object), ('value', float) ])
            
            self.GalLowCut = values[ name  == 'GalLowCut']
            self.GradGal = values[name == 'GradGal']
            self.IntGal = values[ name == 'IntGal']
    
            #Star locuss Parameters
            self.GradStarsLowCut = values[name == 'GradStarsLowCut']
            self.IntStarsLowCut = values[name == 'IntStarsLowCut']
            
    
            self.StarsUpCut = values[ name == 'StarsUpCut']
            self.StarsLowCut = values[ name == 'StarsLowCut']

        def generate_axes( self, sources ):
            self.fig = plt.figure( figsize=(10,10))
            gs = gridspec.GridSpec( 2, 2)
            self.ax3 = plt.subplot( gs[0,0])
            self.ax2 = plt.subplot( gs[1,0])
            self.ax1 = plt.subplot( gs[0,1])
            self.axes = [ self.ax3, self.ax1, self.ax2 ]
            self.ax1.set_xlim([10,35])
            self.ax1.set_ylim([0,25])
            self.ax1.set_xlabel("Magnitude")
            self.ax1.set_ylabel("Size")
            self.ax2.set_xlabel("Magnitude")
            self.ax2.set_ylabel("Radius")
            self.ax2.set_xlim([10,30])
            self.ax2.set_ylim([0,25])
            self.ax3.set_xlabel("Magnitude")
            self.ax3.set_ylabel("Mu Max")
            self.ax3.set_xlim([10,30])
            self.ax3.set_ylim([10,30])


            self.ax1.plot( sources['MAG_AUTO'], sources['gal_size'], 'k,')
            self.ax2.plot( sources['MAG_AUTO'], sources['RADIUS'],   'k,')
            self.ax3.plot( sources['MAG_AUTO'], sources['MU_MAX'],   'k,')
            
        def plot_boundaries( self, sources ):
            self.ax3.plot( sources['MAG_AUTO'], self.GradGal*sources['MAG_AUTO']  + self.IntGal, 'k-')
            self.ax3.plot( sources['MAG_AUTO'], np.zeros(len(sources['MU_MAX']))+ self.GalLowCut, 'k-')
            self.ax3.plot( sources['MAG_AUTO'], sources['MAG_AUTO']*self.GradStarsLowCut + self.IntStarsLowCut, 'y-') 
            self.ax3.plot( sources['MAG_AUTO'], np.zeros(len(sources['MU_MAX']))+  self.StarsLowCut, 'y-')
            self.ax3.plot( sources['MAG_AUTO'], np.zeros(len(sources['MU_MAX']))+ self.StarsUpCut, 'y-')                
            self.ax3.annotate( 'When happy, close plot window', xy=( 0.01, 0.95), \
                                                xycoords='axes fraction', fontsize=10)
                
        def get_params_interactively( self, sources ):

            self.ax3.annotate( 'Use double left click to select point', xy=( 0.01, 0.9), \
                                    xycoords='axes fraction', fontsize=10)
            self.ax3.annotate( 'Use double right click to undo point', xy=( 0.01, 0.85), \
                                    xycoords='axes fraction', fontsize=10)
            self.ann = self.ax3.annotate( 'First draw boundary between galaxies and stars', \
                                              xy=( 0.01, 0.95), weight='bold', \
                                    xycoords='axes fraction', fontsize=10)

                                    
            self.xcoords = []
            self.ycoords = []
            self.borders = []

            
            def onclick(event):
                axes = [event.inaxes, self.ax1, self.ax2]
                if not event.dblclick:
                    return
                
                if event.button == 1:
                    
                    ix, iy = event.xdata, event.ydata
                  
                
                    self.xcoords.append(ix)
                    self.ycoords.append(iy)
                    event.inaxes.plot( ix, iy, 'co')
                
                    event.inaxes.figure.canvas.draw()

                    
                    if len(self.xcoords) == 2:
                        self.borders.append(event.inaxes.plot( self.xcoords, self.ycoords, 'k-'))
                        intercept, grad = coords_to_line( self.xcoords, self.ycoords )
                        self.GradGal = grad
                        self.IntGal = intercept
                        self.ann.remove()
                        self.ann = self.ax3.annotate( 'Second draw lower boundary for galaxies', \
                                                          xy=( 0.01, 0.95), weight='bold', \
                                            xycoords='axes fraction', fontsize=10)

                        self.get_galaxies( sources )
                        self.plot_stars_galaxies(  sources, axes, overwrite=False )
                        event.inaxes.figure.canvas.draw()
                        
                    if len(self.xcoords) == 4:
                        self.borders.append(event.inaxes.plot( self.xcoords[2:], self.ycoords[2:], 'k-'))
                        intercept = (self.ycoords[2] + self.ycoords[3] )/2.
                        self.GalLowCut = intercept
                        self.get_galaxies( sources )
                        self.ann.remove()
                        self.ann = self.ax3.annotate( 'Third draw lower boundary for stars', \
                                                          xy=( 0.01, 0.95), weight='bold', \
                                            xycoords='axes fraction', fontsize=10)

                        self.get_galaxies( sources )
                        self.plot_stars_galaxies(  sources, axes )
                        event.inaxes.figure.canvas.draw()
                    
                    if len(self.xcoords) == 6:
                        
                        self.borders.append(event.inaxes.plot( self.xcoords[4:], self.ycoords[4:], 'y-'))
                        intercept, grad = coords_to_line( self.xcoords[4:], self.ycoords[4:] )
                        self.GradStarsLowCut = grad
                        self.IntStarsLowCut = intercept
                        self.get_stars(sources)
     
                        self.plot_stars_galaxies( sources, axes, overwrite=False )
                        
                        self.ann.remove()
                        self.ann = self.ax3.annotate( 'Fourth draw lower threshold for saturated stars', \
                                                          xy=( 0.01, 0.95), weight='bold', \
                                                xycoords='axes fraction', fontsize=10)
           
                        event.inaxes.figure.canvas.draw()

                    if len(self.xcoords) == 8:
                        self.borders.append(event.inaxes.plot( self.xcoords[6:], self.ycoords[6:], 'y-'))
                        self.StarsLowCut = (self.ycoords[6] + self.ycoords[7])/2.
                        self.get_stars(sources)
                        
                        self.ann.remove()
                        self.plot_stars_galaxies( sources, axes )
                        self.ann = self.ax3.annotate( 'Finally draw upper threshold for noisy stars', \
                                                          xy=( 0.01, 0.95), weight='bold', \
                                                xycoords='axes fraction', fontsize=10)
                                                
                        
                        self.plot_stars_galaxies( sources, axes )
                        event.inaxes.figure.canvas.draw()
                        
                    if len(self.xcoords) == 10:
                        self.borders.append(event.inaxes.plot( self.xcoords[8:], self.ycoords[8:], 'y-'))
                        self.StarsUpCut = (self.ycoords[8] + self.ycoords[9])/2.
                        self.get_stars(sources)
                        
                        self.plot_stars_galaxies( sources, axes)
                        self.ann.remove()
                        self.ann = self.ax3.annotate( 'When happy, close plot window', \
                                                          xy=( 0.01, 0.95), weight='bold',\
                                                xycoords='axes fraction', fontsize=10)
                                                
                   
                        event.inaxes.figure.canvas.draw()
                        
                        
                                        
                    
                else:
                    remove_line = self.borders[-1]
                    remove_line.pop(0).remove()
                    del self.borders[-1]
                    
                    if len(self.xcoords) == 2:
                        self.galaxies = []
                        self.plot_stars_galaxies( sources, axes)
                    elif len(self.xcoords) == 4:
                        self.GalLowCut = 0.
                        self.get_galaxies( sources )
                        self.plot_stars_galaxies( sources, axes)
                    elif len(self.xcoords) == 6:
                        self.stars = []
                        self.plot_stars_galaxies( sources, axes)
                    elif len(self.xcoords) == 8:
                        self.StarsLowCut = 0
                        self.get_stars(sources)
                        self.plot_stars_galaxies( sources, axes )
                    elif len(self.xcoords) == 10:
                        self.StarsUpCut = 100
                        self.get_stars(sources)
                        self.plot_stars_galaxies( sources, axes )

                    self.xcoords = self.xcoords[:-2]
                    self.ycoords = self.ycoords[:-2]
                    event.inaxes.figure.canvas.draw()
                plt.show()


            cid = self.fig.canvas.mpl_connect('button_press_event', onclick)
            plt.show()



        def plot_stars_galaxies( self, sources, axes, overwrite=True ):
            
            if len(self.star_points) > 0:
                iStar_points = self.star_points[ -1 ]
                iGal_points = self.gal_points[ -1]
   
                for i in xrange(len(iStar_points)):
                    iStar_points[i].pop(0).remove()
                    iGal_points[i].pop(0).remove()
                del self.star_points[ -1 ]
                del self.gal_points[ -1 ]
                
                    
            self.star_points.append( [ axes[0].plot( sources['MAG_AUTO'][self.stars], \
                                                         sources['MU_MAX'][self.stars], 'y*'),
                                           axes[1].plot( sources['MAG_AUTO'][self.stars], \
                                                             sources['gal_size'][self.stars], 'y*'),
                                           axes[2].plot( sources['MAG_AUTO'][self.stars], \
                                                             sources['RADIUS'][self.stars], 'y*') ])
            self.gal_points.append( [ axes[0].plot( sources['MAG_AUTO'][self.galaxies], \
                                                        sources['MU_MAX'][self.galaxies], 'r.'  ),
                                          axes[1].plot( sources['MAG_AUTO'][self.galaxies], \
                                                            sources['gal_size'][self.galaxies], 'r.' ),
                                          axes[2].plot( sources['MAG_AUTO'][self.galaxies], \
                                                            sources['RADIUS'][self.galaxies], 'r.' )])




        def get_galaxies( self, sources  ):
            self.galaxies =  (sources['MU_MAX'] > \
                 self.GradGal*sources['MAG_AUTO']  + self.IntGal) & \
                 ( sources['MU_MAX'] > self.GalLowCut )
                                
        def get_stars( self, sources):
            self.stars = (sources['MU_MAX'] < \
                    sources['MAG_AUTO']*self.GradGal + self.IntGal) & \
                            (sources['MU_MAX'] > \
                    sources['MAG_AUTO']*self.GradStarsLowCut + self.IntStarsLowCut) & \
                            (sources['MU_MAX'] > self.StarsLowCut) & \
                            (sources['MU_MAX'] < self.StarsUpCut)

                
def coords_to_line( x, y):

    grad = (y[1] - y[0]) / (x[1] - x[0])

    intercept = y[1] - grad*x[1]

    return intercept, grad
