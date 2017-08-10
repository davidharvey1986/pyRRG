#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Rewite the idl rountine: inpoly.pro, instar.pro to python version 
@author: sutieng
"""
import pyfits as fits
import numpy as np
import matplotlib.pyplot as plt;
import os as os


def instar(xl,yl,xs,ys,m):
    if m>mcut1: 
        l=spike1
        R=r1
    if m>mcut2 and m<=mcut2:    
        l=spike2
        R=r2
    if m>mcut3 and m<=mcut2:    
        l=spike3
        R=r2
    if m<=mcut3:    
        l=spike4
        R=r4
    
    #star x corrdinates
    xstar=np.array([xs-e,xs-e,xs+e,xs+e,xs+R,xs+l,xs+l,xs+R,xs+e,xs+e,xs-e,xs-e,xs-R,xs-l,xs-l,xs-R])
    
    #star y corrdinates
    ystar=np.array([ys+R,ys+l,ys+l,ys+R,ys+e,ys+e,ys-e,ys-e,ys-R,ys-l,ys-l,ys-R,ys-e,ys-e,ys+e,ys+e])
    inside=inpoly(xl,yl,xstar,ystar)
    return xstar, ystar

def inpoly(Px,Py,xl,yl):
#    Determines if a point P(px, py) in inside or outside a polygon.
#   The method used is the ray intersection method based on the
#   Jordan curve theorem. Basically, we trace a "ray" (a semi line)
#   from the point P and we calculate the number of edges that it 
#   intersects. If this number is odd, P is inside.
#    
#   (Px,Py)    Coordinates of point
#   xv         List of x coordinates of vertices of polygon
#   yv         List of y coordinates of vertices of polygon
#
#   The x and y lists may be given clockwise or anticlockwise.
#   The algorithm works for both convex and concave polygons.
    N=len(xl)
    xv=np.zeros(N)
    yv=np.zeros(N)
    
    for i in np.arange(N):
        xv[i]=xl[i]
        yv[i]=yl[i]
        
    nc=0        #Number of edge crossings
    N=len(xv)   #Number of vertices
    
    #test input
    if N<3:
        print "A polygon must have at least three vertices"
    if len(xv)!=len(yv):
        print 'Must have same number of X and Y coordinates'
        
    #---------------------- Change coordinate system -----------------
    #Place P at the center of the coordinate system.
    
    for i in np.arange(N):
        xv[i]=xv[i]-Px
        yv[i]=yv[i]-Py
          
    #---------------------- Calculate crossings ----------------------          
    #The positive half of the x axis is chosen as the ray
    #We must determine how many edges cross the x axis with x>0
    for i in np.arange(N):
        Ax=xv[i]    #first vertice of edge
        Ay=yv[i]
        
        if i==(N-1):
            Bx=xv[0]
            By=yv[0]
        else: 
            Bx=xv[i+1]  #Second vertice of edge
            By=yv[i+1]

    #We define two regions in the plan: R1/ y<0 and R2/ y>=0. Note that
    #the case y=0 (vertice on the ray) is included in R2.
        if Ay<0:
            signA=-1
        else:   signA=+1
        if By<0:
            signB=-1
        else:   signB=+1
    
    #The edge crosses the ray only if A and B are in different regions.
    #If a vertice is only the ray, the number of crossings will still be
    #correct.
    
        if (signA*signB<0):
            if (Ax>0 and Bx>0): nc+=1
            else:   
                x=Ax-(Ay*(Bx-Ax))/(By-Ay)
                if x>0: nc+=1
 


#if inside then uneven
#if outside then even               
    nc=nc%2
    return nc


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

##define some parameters------------------------------------------------------------------------
# Spikes of stars in pixels
e=18                          #width of spike

spike1=200      	      #length of first spike
spike2=500                    # length of second spike
spike3=900                    # length of third spike
spike4=1500           	      # length of fourth spike

r1=50               # radius of mask
r2=120
r3=200
r4=260

mcut1=16.8            		# magnitude cut
mcut2=15.9              	# magnitude cut
mcut3=13.8	      		# magnitude cut

r_close=1.5                     #radius for looking for bright objects in arcmin
#---------------------------------------------------------------------------------------------
#-------------------------select the stars from the size vs mag diagram----------------------------------------------------



savefile='galStar.locus'
if not os.path.isfile( savefile):
    gal_star = galStar( savefile, set_defaults=True)
    print 'Loading defaults '+savefile
else:
    print 'Loading '+savefile
    gal_star = galStar( savefile, set_defaults=False) 

all_sources=fits.open('final_F814W_drz_sci_uncor.cat')[1].data  ##output from SExtractor
stars = [(all_sources['gal_size'] < all_sources['MAG_AUTO']*gal_star.GradStars +  gal_star.IntStars) &
                            (all_sources['gal_size'] > gal_star.StarsLowCut )]

Star_catalogue=all_sources[stars]
fits.writeto('Star.fits',Star_catalogue)

##---------------------add a new column 'clean' to shear catalogue---------------------------------
data=fits.open('final_F814W_drz_sci_10.shears')[1].data   ##remember to change it to the name of your shear catalogue
clean=np.zeros(len(data['ra']))
cols = [] 
cols.append(
    fits.Column(name='clean', format='D', array= clean)
    )
orig_cols = data.columns
new_cols = fits.ColDefs(cols)
hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)
hdu.writeto('final_F814W_drz_sci_clean.shears')


##-------------------------------start masking-------------------------------
Shears=fits.open('final_F814W_drz_sci_clean.shears')[1].data
##go through the sources list:
for i in np.arange(len(Shears['ra'])):
    xl=Shears['X_IMAGE'][i]
    yl=Shears["Y_IMAGE"][i]
    m=Shears["MAG_AUTO"][i]
    for j in np.arange(len(Star_catalogue["ra"])):   #go through the star list
        star_x=Star_catalogue['X_IMAGE'][j]
        star_y=Star_catalogue['Y_IMAGE'][j]
        inside=instar(xl,yl,star_x,star_y,m)
        if inside==1:
            Shears['clean'][i]=1
            break
        
Shears_remove=Shears[Shears['clean']==0]
fits.writeto('Shear_remove.fits',Shears_remove )


##------------------------------polygen defined by hand------------------------------
#For large objects, please define the polygen by hand, enter the coordinate in px1,py1
'''
px1=np.array([])  ##please enter the X-image corrdinate of the polygen
py1=np.array([])  ##please enter the Y-image corrdinate of the polygen

print px1,py1
Shears_2=fits.open('Shear_remove.fits')[1].data
for i in np.arange(len(Shears_2['ra'])):
    xl=Shears_2['X_IMAGE'][i]
    yl=Shears_2["Y_IMAGE"][i]
    m=Shears_2["MAG_AUTO"][i]
    inside=inpoly(xl,yl,px1,py1)
    if inside==1:
            Shears_2['clean'][i]=1


Shears_remove=Shears_2[Shears_2['clean']==0]
fits.writeto('Shear_remove2.fits',Shears_remove )
'''


