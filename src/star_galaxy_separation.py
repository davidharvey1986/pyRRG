from matplotlib import pyplot as plt
import os as os
import pickle as pkl
from matplotlib import gridspec as gridspec
import sys as sys
import numpy as np

def star_galaxy_separation( sources, restore=False,
                            savefile='galStar.sav'):

    '''
    PURPORSE : To separate out the stars and galaxies in the input
               catalogue
    
    INPUTS : sources : a sextractor catalogue of sources (fits record)
    
    KEYWORDS : 
     restore : restore a saved file of classes without prompting user 
               to re-separate stars anda galaixes
     savefile : save the file to the given string
    
    RETURNS : A vector of N galaxies and stars 
              with the indexes associated with SOURCES
              

    METHOD  : If the savefile does not exist it will load the sources into 3 plot
               1 ) MAG vs MU MAX the main plot where user will separate stars and galaxies
               2 ) MAG vs SIZE as dedcued from mometns
               3 ) MAG vs RADUIS as deduced from source extractor

    AUTHOR : David Harvey
    
    INFO : For more information on how this works and how it is used by the user please see 
           the documentation
    '''
  
    object_indexes = np.arange( len(sources['X_IMAGE']) )
  
    if not os.path.isfile( savefile ):
        #If the file does not exist get the user inputs for stars / galaxies
        print('Save file not found, separating, please follow instructions')
        print('found in the top left plot. Please only select here')
        gal_star = galStar( savefile, sources, set_defaults=True)
        return object_indexes[ gal_star.galaxies], object_indexes[ gal_star.stars ]
    else:
        if restore:
            #Dont ask, just get the stars and galaxies and return
            print 'Restoring '+savefile
            gal_star = galStar( savefile, sources, set_defaults=False)   
            gal_star.get_galaxies( sources )
            gal_star.get_stars( sources )
            return object_indexes[ gal_star.galaxies], object_indexes[ gal_star.stars ]
        
        print 'Loading savefile: '+savefile
        gal_star = galStar( savefile, sources, set_defaults=False)   
        gal_star.get_galaxies( sources )
        gal_star.get_stars( sources )
        gal_star.generate_axes( sources )
        gal_star.plot_boundaries(  sources )
        gal_star.plot_stars_galaxies( sources, gal_star.axes, overwrite=True )
        plt.show()

        
        overwrite = \
          raw_input('Save file exists, accept or reject current selection?\n'+\
                        'Yes (y)      : And remove all files and remeasure stars and galaxies\n'+\
                        'Continue (c) : Use current selection, do not remeasure galaxies and stars\n'+\
                        'No (n)       : Reject start-gal separation and re-do\n'+\
                        '>>> ')
                  
        
        if overwrite == 'y':
            print 'Writing over gal file and removing j*uncor.cat'
            os.system('rm -fr j*uncor.cat *_cor.cat')
            gal_star.write( savefile )
            return object_indexes[ gal_star.galaxies], object_indexes[ gal_star.stars ]
        else:
            if overwrite == 'c':
                return object_indexes[ gal_star.galaxies], object_indexes[ gal_star.stars ]
            else:
                #If i reseparate i will want to remeasure all stars and galaxies
                print 'Reseparating'
                os.system('rm -fr '+savefile)
                os.system('rm -fr j*uncor.cat *_cor.cat')
                galaxies, stars = star_galaxy_separation( sources, savefile=savefile )
                return galaxies, stars

    
#This class is where all the meat of the code is run including the interactive plotting
class galStar():

        def __init__( self, filename, sources, set_defaults=True ):

            #These store the plotted points and the order they are plotted
            #so they can be removed if the user decideds to undo a line
            self.star_points = []
            self.gal_points = []
            
            if set_defaults:
                self.generate_axes( sources)
                self.defaults( filename, sources )
            else:
                self.load( filename )


        def defaults( self, filename, sources):
            '''
            This function is called if we have no savefile
            It uses interactive plotting to get the default settings
            
            INPUT :
                filename : the name of string that the settings will be saved to
                sources : a fits_record of the objects in the catalogue
            '''
            self.StarsLowCut = 0.
            self.StarsUpCut = 100
            self.stars = []
            self.GalLowCut = 0

            self.get_params_interactively( sources ) 

            self.write( filename )
            
        def write(self, filename):
            '''
            write the paramters that define the boundaries
            for the galaxies and stars 
            INPUTS :
               filename : a string for the name of the file to save
                          the parameters to
            '''
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
            '''
            Load the paramters that define the boundaries for the 
            stars and galaxies from the save file
            
            INPUTS :
                filename : string of the savefile filename
            '''
            
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
            '''
            Generate the axes and label the three plots 
            that will be used to separaet the stars and galaxies
            
            INPUTS :
                sources : a fits_record of sources from source extractor
                
            '''
            
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
            '''
            Function to draw the boundaries for the stars and agalaxies
            All on the mu max plot 
            1. ) The separation between the galaxies and stars
            2. ) The lower cut (constant in MU MAX) in MU MAX for the galaxies
            3. ) The lower cut (grad and intercept) for the stars
            4. ) THe lower cut (constant MU MAX) to remove saturated stars
            5. ) The upper cut (constant MU MAX) to remove noisy stars
            '''
            
            self.ax3.plot( sources['MAG_AUTO'], self.GradGal*sources['MAG_AUTO']  + self.IntGal, 'k-')
            self.ax3.plot( sources['MAG_AUTO'], np.zeros(len(sources['MU_MAX']))+ self.GalLowCut, 'k-')
            self.ax3.plot( sources['MAG_AUTO'], sources['MAG_AUTO']*self.GradStarsLowCut + self.IntStarsLowCut, 'y-') 
            self.ax3.plot( sources['MAG_AUTO'], np.zeros(len(sources['MU_MAX']))+  self.StarsLowCut, 'y-')
            self.ax3.plot( sources['MAG_AUTO'], np.zeros(len(sources['MU_MAX']))+ self.StarsUpCut, 'y-')                
            self.ax3.annotate( 'When happy, close plot window', xy=( 0.01, 0.95), \
                                                xycoords='axes fraction', fontsize=10)
                
        def get_params_interactively( self, sources ):
            '''
            Runs the plotting algorithm that lets the user draw lines
            For more information see the documentation

            INPUT : 
                sources : fits record of the sources from source exractor 
            '''
            
            #Write on the plot some useful directions for the user
            self.ax3.annotate( 'Use double left click to select point', xy=( 0.01, 0.9), \
                                    xycoords='axes fraction', fontsize=10)
            self.ax3.annotate( 'Use double right click to undo point', xy=( 0.01, 0.85), \
                                    xycoords='axes fraction', fontsize=10)
            self.ann = self.ax3.annotate( 'First draw boundary between galaxies and stars', \
                                              xy=( 0.01, 0.95), weight='bold', \
                                    xycoords='axes fraction', fontsize=10)


            #Set paramters that will record where  the user clicks
            self.xcoords = []
            self.ycoords = []
            #Set the list of lines that will be used so they can be undone
            self.borders = []

            
            def onclick(event):
                '''
                Function that decides what to do when the user clicks on the plot
                '''
                #Determine which axis is being used (at the moment only the first
                axes = [event.inaxes, self.ax1, self.ax2]

                #The user must double click allowng the user to zoom and scroll in the image
                if not event.dblclick:
                    return
                #It also must be a left button click
                if event.button == 1:

                    #Where the user clicks
                    ix, iy = event.xdata, event.ydata
                  
                
                    self.xcoords.append(ix)
                    self.ycoords.append(iy)
                    event.inaxes.plot( ix, iy, 'co')

                    #Plot the point clicked by the user
                    event.inaxes.figure.canvas.draw()

                    
                    if len(self.xcoords) == 2:
                        #Once there are two points plot this
                        self.borders.append(event.inaxes.plot( self.xcoords, self.ycoords, 'k-'))
                        #Determine what this is in terms of a gradient and intercept
                        intercept, grad = coords_to_line( self.xcoords, self.ycoords )
                        #Store as paramters 
                        self.GradGal = grad
                        self.IntGal = intercept
                        #Remove the first help message and replace with the next
                        self.ann.remove()
                        self.ann = self.ax3.annotate( 'Second draw lower boundary for galaxies', \
                                                          xy=( 0.01, 0.95), weight='bold', \
                                            xycoords='axes fraction', fontsize=10)
                        #Get which galaxies the user has chosen
                        self.get_galaxies( sources )
                        #Plot the galaxies on top
                        self.plot_stars_galaxies(  sources, axes, overwrite=False )
                        event.inaxes.figure.canvas.draw()
                        
                    if len(self.xcoords) == 4:
                        #Now the user has selected the lower threshold for the galaxies plot the border
                        self.borders.append(event.inaxes.plot( self.xcoords[2:], self.ycoords[2:], 'k-'))
                        #The border is constant MU MAX so take the mean of the two ycoords
                        intercept = (self.ycoords[2] + self.ycoords[3] )/2.
                        #Set  in the galaxy param file
                        self.GalLowCut = intercept
                        self.get_galaxies( sources )
                        #Remove help message, show next
                        self.ann.remove()
                        self.ann = self.ax3.annotate( 'Third draw lower boundary for stars', \
                                                          xy=( 0.01, 0.95), weight='bold', \
                                            xycoords='axes fraction', fontsize=10)
                        #update galaxy selection and show
                        self.get_galaxies( sources )
                        self.plot_stars_galaxies(  sources, axes )
                        event.inaxes.figure.canvas.draw()
                    
                    if len(self.xcoords) == 6:
                        #Now get the lower threshold of the stars
                        self.borders.append(event.inaxes.plot( self.xcoords[4:], self.ycoords[4:], 'y-'))
                        intercept, grad = coords_to_line( self.xcoords[4:], self.ycoords[4:] )
                        self.GradStarsLowCut = grad
                        self.IntStarsLowCut = intercept
                        
                        #Now get the current star selection and show
                        self.get_stars(sources)
                        self.plot_stars_galaxies( sources, axes, overwrite=False )
                        
                        #Remove help message, show next
                        self.ann.remove()
                        self.ann = self.ax3.annotate( 'Fourth draw lower threshold for saturated stars', \
                                                          xy=( 0.01, 0.95), weight='bold', \
                                                xycoords='axes fraction', fontsize=10)
           
                        event.inaxes.figure.canvas.draw()

                    if len(self.xcoords) == 8:
                        #Get the lower MU AX for saturated stars
                        self.borders.append(event.inaxes.plot( self.xcoords[6:], self.ycoords[6:], 'y-'))
                        self.StarsLowCut = (self.ycoords[6] + self.ycoords[7])/2.

                        #Update star selection and plot
                        self.get_stars(sources)
                        self.plot_stars_galaxies( sources, axes )
                        
                        #Remove help message, show next
                        self.ann.remove()
                        self.ann = self.ax3.annotate( 'Finally draw upper threshold for noisy stars', \
                                                          xy=( 0.01, 0.95), weight='bold', \
                                                xycoords='axes fraction', fontsize=10)
                                                
                        
                        self.plot_stars_galaxies( sources, axes )
                        event.inaxes.figure.canvas.draw()
                        
                    if len(self.xcoords) == 10:
                        #Finally get the noisy star threshold
                        self.borders.append(event.inaxes.plot( self.xcoords[8:], self.ycoords[8:], 'y-'))
                        self.StarsUpCut = (self.ycoords[8] + self.ycoords[9])/2.

                        #Update star selection and show
                        self.get_stars(sources)
                        self.plot_stars_galaxies( sources, axes)

                        #Remove help message, show next
                        self.ann.remove()
                        self.ann = self.ax3.annotate( 'When happy, close plot window', \
                                                          xy=( 0.01, 0.95), weight='bold',\
                                                xycoords='axes fraction', fontsize=10)
                                                
                   
                        event.inaxes.figure.canvas.draw()
                        
                        
                                        
                    
                else:
                    #If the double click was the right hand button
                    #Remove the most recent border (line)
                    remove_line = self.borders[-1]
                    remove_line.pop(0).remove()
                    del self.borders[-1]

                    #Then reset the relevant paramters
                    #Reset galaxies or stars
                    #Replot the galaxies or stars
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
                    #Remove the coords from x and ycoords
                    self.xcoords = self.xcoords[:-2]
                    self.ycoords = self.ycoords[:-2]
                    event.inaxes.figure.canvas.draw()
                plt.show()

            #tell python to have interactive plotting using the onclick function
            cid = self.fig.canvas.mpl_connect('button_press_event', onclick)
            plt.show()



        def plot_stars_galaxies( self, sources, axes, overwrite=True ):
            '''
            Plot all the objects
            Plot the stars in yellow stars
            Plot the galaxies in red points
            '''
            
            if len(self.star_points) > 0:
                iStar_points = self.star_points[ -1 ]
                iGal_points = self.gal_points[ -1]
   
                for i in xrange(len(iStar_points)):
                    #First remove any overlaid galaxy or star points alrady
                    #on the plot
                    iStar_points[i].pop(0).remove()
                    iGal_points[i].pop(0).remove()
                del self.star_points[ -1 ]
                del self.gal_points[ -1 ]
                
            #Now plot the stars and galaxies and store the plotted points in
            #star_points and gal_points
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
            '''
            Using the boundaries set determien which of the sourecs
            are galaxies. Store as a boolean array of N objects
            
            INPUT :
               sources : a fits record of the sources from source extractor
            '''
            
            self.galaxies =  (sources['MU_MAX'] > \
                 self.GradGal*sources['MAG_AUTO']  + self.IntGal) & \
                 ( sources['MU_MAX'] > self.GalLowCut )
                                
        def get_stars( self, sources):
            '''
            Using the boundaries set determien which of the sourecs
            are stars. Store as a boolean array of N objects
            
            INPUT :
               sources : a fits record of the sources from source extractor
            '''
            self.stars = (sources['MU_MAX'] < \
                        sources['MAG_AUTO']*self.GradGal + self.IntGal) & \
                            (sources['MU_MAX'] > \
                    sources['MAG_AUTO']*self.GradStarsLowCut + self.IntStarsLowCut) & \
                            (sources['MU_MAX'] > self.StarsLowCut) & \
                            (sources['MU_MAX'] < self.StarsUpCut)

                
def coords_to_line( x, y):
    '''
    Given two coordaintes, determine the gradient and intercept of the 
    lie
    INPUTS :
       x : a two point vector 
       y : a two point vector
    '''
    
    grad = (y[1] - y[0]) / (x[1] - x[0])

    intercept = y[1] - grad*x[1]

    return intercept, grad
