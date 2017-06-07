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
  
     
    if not os.path.isfile( savefile ):
        gal_star = galStar( savefile, sources, set_defaults=True)

        print 'Loading defaults '+savefile
    else:
        print 'Loading '+savefile
        gal_star = galStar( savefile, sources, set_defaults=False)   
   
  

     
    #Original number ; 1, -3.3, -9, 1, -3.9, -3
   
    galaxies = object_indexes[ (sources['MU_MAX'] > gal_star.GradGal*sources['MAG_AUTO']  + \
                                gal_star.IntGal) &
                                ( sources['MU_MAX'] > gal_star.GalLowCut) & \
                                ( sources['gal_size'] > gal_star.GradStars*sources['MAG_AUTO'] +\
                                gal_star.IntStars )]
                            
    n_galaxies = len(galaxies)
 
    
    #this is just a box aroun dhe stars in the mag vs radius for the 
    #the drizzled image
    stars = object_indexes[ (sources['MU_MAX'] < sources['MAG_AUTO']*gal_star.GradGal + gal_star.IntGal) & \
                             (sources['MU_MAX'] > sources['MAG_AUTO']*gal_star.GradStarsLowCut + gal_star.IntStarsLowCut) & \
                             (sources['gal_size'] < sources['MAG_AUTO']*gal_star.GradStars +  gal_star.IntStars) &\
                            (sources['gal_size'] > gal_star.StarsLowCut )]
    
    n_stars = len(stars)

    if restore:
        return  galaxies, stars
    
    if savefile is not None:

        print 'There are '+str(n_stars)+ ' in the field now plotting'
        #Plot the stars and galaxies to make sure i am taking the right ones
        
        ax1.plot( sources['MAG_AUTO'], sources['gal_size'], 'k,')
        ax1.set_xlabel("Magnitude")
        ax1.set_ylabel("Size")
       

        
        ax1.plot( sources['MAG_AUTO'], gal_star.GradStars*\
                    sources['MAG_AUTO'] +gal_star.IntStars, 'g-' )

        ax1.plot( sources['MAG_AUTO'], np.zeros(nObjects)+gal_star.StarsUpCut,'g-')
         
        ax1.plot( sources['MAG_AUTO'], np.zeros(nObjects)+gal_star.StarsLowCut, 'y-')
     
        ax1.plot( sources['MAG_AUTO'][galaxies], sources['gal_size'][galaxies], 'rD', label='Galaxies')
        ax1.plot( sources['MAG_AUTO'][stars], sources['gal_size'][stars],'y*', label='Stars')
        
        ax1.set_xlim([10,35])
        ax1.set_ylim([0,25])
      
        ax1.legend()
        
    
     
        ax2.plot( sources['MAG_AUTO'], sources['RADIUS'],'k,')

        ax2.plot( sources['MAG_AUTO'][galaxies], sources['RADIUS'][galaxies], 'rD')
        ax2.plot( sources['MAG_AUTO'] [stars], sources['RADIUS'][stars], 'y*')

        ax2.set_xlabel("Magnitude")
        ax2.set_ylabel("Radius")
        ax2.set_xlim([10,30])
        ax2.set_ylim([0,25])
     
     
        ax3.plot( sources['MAG_AUTO'] ,sources['MU_MAX'], 'k,')




        ax3.plot( sources['MAG_AUTO'], np.zeros(nObjects)+gal_star.GalLowCut, 'r-')
        ax3.plot( sources['MAG_AUTO'], gal_star.GradGal*\
                    sources['MAG_AUTO'] + gal_star.IntGal, 'g-')
                    
        ax3.plot( sources['MAG_AUTO'], sources['MAG_AUTO']*\
                gal_star.GradStarsLowCut + gal_star.IntStarsLowCut, 'r-')
        
        
        ax3.plot( sources['MAG_AUTO'][galaxies],sources['MU_MAX'][galaxies],'rD')
        ax3.plot( sources['MAG_AUTO'][stars],sources['MU_MAX'][stars], 'y*')
        ax3.set_xlabel("Magnitude")
        ax3.set_ylabel("Mu Max")
        ax3.set_xlim([10,30])
        ax3.set_ylim([10,30])

        plt.savefig( 'StarGalaxy.pdf' )
        if axes is None:
            plt.show( block=False)
        else:
            plt.draw()
        
       
     

        if NonStop:
            overwrite='c'
        else:
            overwrite = \
                raw_input('Save file exists, overwrite Yes (y) Continue (c) No (n)?')
                  
        
        if overwrite == 'y':
            print 'Writing over gal file and removing j*uncor.cat'
            os.system('rm -fr j*uncor.cat *_cor.cat')
            gal_star.write( savefile )

        else:
            if overwrite == 'c':
                cont=1
                plt.close()
                return galaxies, stars
            else:
                print 'Reseparating'
                galaxies, stars = star_galaxy_separation( sources, savefile=savefile, axes=[ax1,ax2,ax3] )

    plt.close()
    return galaxies, stars

    


class galStar():

        def __init__( self, filename, sources, set_defaults=True ):
            if set_defaults:
                self.defaults( filename, sources )
            else:
                self.load( filename )


        def defaults( self, filename, sources):
            #Use default parameters (first guess)
            #galaxy locuss parameters
            self.StarsLowCut = 0.
            self.StarsUpCut = 100

            self.get_params_interactively( sources ) 

            self.GalLowCut = 16
            self.GradGal = 1
            self.IntGal = -4
    
            #Star locuss Parameters
            self.GradStars = -1.1
            self.IntStars = 28.5
        
            self.GradStarsLowCut = 1.
            self.IntStarsLowCut = -5.
    
            self.write( filename )
            
        def write(self, filename):
            galstar_file = open( filename, 'wb')
            pdb.set_trace()
            galstar_file.write('#Galaxy Parameters\n')
            galstar_file.write('GalLowCut %5.2f\n' % self.GalLowCut)
            galstar_file.write('GradGal %5.2f\n' % self.GradGal)
            galstar_file.write('IntGal %5.2f\n' % self.IntGal)
            galstar_file.write('#Star Parameters\n')
            galstar_file.write('GradStars %5.2f\n' % self.GradStars)
            galstar_file.write('IntStars %5.2f\n' % self.IntStars)
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
            self.GradStars = values[name == 'GradStars']
            self.IntStars = values[name == 'IntStars']
            
            self.GradStarsLowCut = values[name == 'GradStarsLowCut']
            self.IntStarsLowCut = values[name == 'IntStarsLowCut']
    
            self.StarsUpCut = values[ name == 'StarsUpCut']
            self.StarsLowCut = values[ name == 'StarsLowCut']


        def get_params_interactively( self, sources ):
            fig = plt.figure( figsize=(10,10))
            gs = gridspec.GridSpec( 2, 2)
            ax3 = plt.subplot( gs[0,0])
            ax2 = plt.subplot( gs[1,0])
            ax1 = plt.subplot( gs[0,1])

            ax1.plot( sources['MAG_AUTO'], sources['gal_size'], 'k,')
            ax2.plot( sources['MAG_AUTO'], sources['RADIUS'],   'k,')
            ax3.plot( sources['MAG_AUTO'], sources['MU_MAX'],   'k,')

            ax1.set_xlim([10,35])
            ax1.set_ylim([0,25])
            ax1.set_xlabel("Magnitude")
            ax1.set_ylabel("Size")
            ax2.set_xlabel("Magnitude")
            ax2.set_ylabel("Radius")
            ax2.set_xlim([10,30])
            ax2.set_ylim([0,25])
            ax3.set_xlabel("Magnitude")
            ax3.set_ylabel("Mu Max")
            ax3.set_xlim([10,30])
            ax3.set_ylim([10,30])
            
            ann = ax3.annotate( 'First draw upper boundary for stars', xy=( 0.01, 0.95), \
                                    xycoords='axes fraction', fontsize=10)
            global xcoords
            global ycoords
            xcoords = []
            ycoords = []
    
            def onclick(event):
                
                if not event.dblclick:
                    return
                    
                global ix, iy
                ix, iy = event.xdata, event.ydata
                print 'x = %d, y = %d'%(
                        ix, iy)

                
                xcoords.append(ix)
                ycoords.append(iy)
                event.inaxes.plot( ix, iy, '*')
                
                event.inaxes.figure.canvas.draw()

                
                if len(xcoords) == 2:
                    event.inaxes.plot( xcoords, ycoords, '-')
                    intercept, grad = coords_to_line( xcoords, ycoords )
                    self.GradGal = grad
                    self.IntGal = intercept
                    ann.remove()
                    global ann1
                    ann1 = ax3.annotate( 'Second draw lower boundary for stars', xy=( 0.01, 0.95), \
                                            xycoords='axes fraction', fontsize=10)
                    event.inaxes.figure.canvas.draw()
                    
                    
                if len(xcoords) == 4:
                    ann1.remove()
                    event.inaxes.plot( xcoords[2:], ycoords[2:], '-')
                    intercept, grad = coords_to_line( xcoords[2:], ycoords[2:] )
                    self.GradStarsLowCut = grad
                    self.IntStarsLowCut = intercept
                    stars = get_stars(sources, self)
                    global points_ax1
                    global points_ax2
                    global points_ax3
                    points_ax3 = event.inaxes.plot( sources['MAG_AUTO'][stars], sources['MU_MAX'][stars], 'y*')
                    points_ax1 = ax1.plot( sources['MAG_AUTO'][stars], sources['gal_size'][stars], 'y*')
                    points_ax2 = ax2.plot( sources['MAG_AUTO'][stars], sources['RADIUS'][stars], 'y*')
                    
                    global ann2
                    ann2 = ax3.annotate( 'Third draw lower threshold for saturated stars', xy=( 0.01, 0.95), \
                                            xycoords='axes fraction', fontsize=10)
           
                    event.inaxes.figure.canvas.draw()

                if len(xcoords) == 6:
                    event.inaxes.plot( xcoords[4:], ycoords[4:], '-')
                    self.StarsLowCut = (ycoords[4] + ycoords[5])/2.
                    stars = get_stars(sources, self)
                    
                    points_ax3.pop(0).remove()
                    points_ax2.pop(0).remove()
                    points_ax1.pop(0).remove()

                    points_ax3 = event.inaxes.plot( sources['MAG_AUTO'][stars], sources['MU_MAX'][stars], 'y*')
                    points_ax1 = ax1.plot( sources['MAG_AUTO'][stars], sources['gal_size'][stars], 'y*')
                    points_ax2 = ax2.plot( sources['MAG_AUTO'][stars], sources['RADIUS'][stars], 'y*')
                    event.inaxes.figure.canvas.draw()
                    ann2.remove()
                    ann3 = ax3.annotate( 'Third draw upper threshold for noisy stars', xy=( 0.01, 0.95), \
                                            xycoords='axes fraction', fontsize=10)



                if len(xcoords) == 8:
                    event.inaxes.plot( xcoords[6:], ycoords[6:], '-')
                    self.StarsUpCut = (ycoords[6] + ycoords[7])/2.
                    stars = get_stars(sources, self)
                    points_ax3.pop(0).remove()
                    points_ax2.pop(0).remove()
                    points_ax1.pop(0).remove()

                    points_ax3 = event.inaxes.plot( sources['MAG_AUTO'][stars], sources['MU_MAX'][stars], 'y*')
                    points_ax1 = ax1.plot( sources['MAG_AUTO'][stars], sources['gal_size'][stars], 'y*')
                    points_ax2 = ax2.plot( sources['MAG_AUTO'][stars], sources['RADIUS'][stars], 'y*')
                    event.inaxes.figure.canvas.draw()



                    
                    fig.canvas.mpl_disconnect(cid)
                    

                    
                plt.show()
                return xcoords, ycoords

            
 #fig.canvas.mpl_connect('axes_enter_event', enter_axes)
#fig.canvas.mpl_connect('axes_leave_event', leave_axes)

            cid = fig.canvas.mpl_connect('button_press_event', onclick)
            plt.show()
            self.x = xcoords
            self.y = ycoords



def coords_to_line( x, y):

    grad = (y[1] - y[0]) / (x[1] - x[0])

    intercept = y[1] - grad*x[1]

    return intercept, grad


def get_galaxies( sources, galStar ):
    
    return  (sources['MU_MAX'] > \
                 gal_star.GradGal*sources['MAG_AUTO']  + gal_star.IntGal)

                                #&
                                #( sources['MU_MAX'] > gal_star.GalLowCut) & \
                                #( sources['gal_size'] > gal_star.GradStars*sources['MAG_AUTO'] +\
                                #gal_star.IntStars )]
def get_stars( sources, galStar):
    
    return (sources['MU_MAX'] < \
                sources['MAG_AUTO']*galStar.GradGal + galStar.IntGal) & \
            (sources['MU_MAX'] > \
                 sources['MAG_AUTO']*galStar.GradStarsLowCut + galStar.IntStarsLowCut) & \
            (sources['MU_MAX'] > galStar.StarsLowCut) & (sources['MU_MAX'] < galStar.StarsUpCut)
#& \
  #                           (sources['gal_size'] < sources['MAG_AUTO']*gal_star.GradStars +  gal_star.IntStars) &\
   #                         (sources['gal_size'] > gal_star.StarsLowCut )]
