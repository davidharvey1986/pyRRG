from matplotlib import pyplot as plt
import os as os
import pickle as pkl
from matplotlib import pyplot as plt
from matplotlib import gridspec as gridspec
import sys as sys
import numpy as np
import ipdb as pdb
def star_galaxy_separation( sources, restore=False,
                            savefile=None,
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

    #There will be three plots
    if not restore:
        if axes is None :
            fig = plt.figure( figsize=(10,10))
            gs = gridspec.GridSpec( 2, 2)
            ax1 = plt.subplot( gs[0,0])
            ax2 = plt.subplot( gs[1,0])
            ax3 = plt.subplot( gs[0,1])
        else:
            ax1, ax2, ax3 = axes
            ax1.clear()
            ax2.clear()
            ax3.clear()
        
    radius = np.sqrt(sources['ISOAREA_IMAGE']/np.pi)
    
    #Create a box that separateds the galaxies fom teh stars
    #Each box is define by
    #The line separating Stars and Galaxies
    
  
    nObjects = len(sources['X_IMAGE'])
    object_indexes = np.arange( nObjects )
  
     
    if not os.path.isfile( savefile ):
        gal_star = galStar( savefile, set_defaults=True)
        print 'Loading defaults '+savefile
    else:
        print 'Loading '+savefile
        gal_star = galStar( savefile, set_defaults=False)   
   
  

     
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

        def __init__( self, filename, set_defaults=True ):
            if set_defaults:
                self.defaults( filename )
            else:
                self.load( filename )


        def defaults( self, filename):
            #Use default parameters (first guess)
            #galaxy locuss parameters
            self.GalLowCut = 16
            self.GradGal = 1
            self.IntGal = -4
    
            #Star locuss Parameters
            self.GradStars = -1.1
            self.IntStars = 28.5
        
            self.GradStarsLowCut = 1.
            self.IntStarsLowCut = -5.
    
            self.StarsUpCut = 2.5
            self.StarsLowCut = 2.

            self.write( filename )
            
        def write(self, filename):
            galstar_file = open( filename, 'wb')
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


