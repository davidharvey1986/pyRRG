import pyfits as fits
import numpy as np
import matplotlib.pyplot as plt
import os as os
import star_galaxy_separation as sgs
import ipdb as pdb
import RRGtools as tools


def plot_region_maskstar( filename, star):
    regionFile = open( filename, "wb")
    regionFile.write('# Region file format: DS9 version 4.1\n')
    regionFile.write("global color=green dashlist=8 3 width=1 font='helvetica 10 normal roman' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n")
    regionFile.write("image\n")
    for j in xrange(len(star)):
        polygonStr = 'polygon('
        for i in star[j]:
            polygonStr += '%0.4f,'  % i
        polygonStr=polygonStr[:-2]  ##delete , in the end
        polygonStr += ')\n'
        regionFile.write(polygonStr)

                             


def instar(xl,yl,xs,ys,m):
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
    star=[]
    for i in np.arange(len(xstar)):
        star.append(xstar[i])
        star.append(ystar[i])
    return inside, star


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


def main(  shear_catalog, object_catalog_fits, \
         mask_file='mask.reg', outFile='Shear_remove.fits' ,plot_reg=False):
    '''
        This algoritm will do two things.
        a) Draw masks around detected stars automatically
        b) Take an additional optional mask file (defaulted to mask.reg) and remove
        any object that lies within those regions. This file MUST be a ds9 region file with a
        formatting with degrees and not sexidecimal.
        '''
    
    #-------------------------select the stars from the size vs mag diagram----------------------------------------------------
    
    object_catalog = fits.open(object_catalog_fits)[1].data
    
    galaxies, stars = sgs.star_galaxy_separation( object_catalog,
                                                 savefile='galStar.locus',
                                                 restore=True,
                                                 include_sat=True)
        
    Star_catalogue = object_catalog[ stars]
     
     
     ##---------------------add a new column 'clean' to shear catalogue---------------------------------
    data=fits.open(shear_catalog)[1].data   ##remember to change it to the name of your shear catalogue
    clean=np.zeros(len(data['ra']))
    cols = []
    cols.append(
                 fits.Column(name='clean', format='D', array= clean)
                 )
    orig_cols = data.columns
    new_cols = fits.ColDefs(cols)
    hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)
    clean_catalog = shear_catalog.split('.')[0]+'_clean.'+\
        shear_catalog.split('.')[1]
    hdu.writeto(clean_catalog, clobber=True,output_verify='ignore')
    

    star_corr=[[] for i in np.arange(len(Star_catalogue["ra"]))]
##-------------------------------start masking-------------------------------
    Shears=fits.open(clean_catalog)[1].data
    ##go through the sources list:
    for i in np.arange(len(Shears['ra'])):
        xl=Shears['X_IMAGE'][i]
        yl=Shears["Y_IMAGE"][i]
        m=Shears["MAG_AUTO"][i]
        for j in np.arange(len(Star_catalogue["ra"])):   #go through the star list
            star_x=Star_catalogue['X_IMAGE'][j]
            star_y=Star_catalogue['Y_IMAGE'][j]
            inside,star_corr_one=instar(xl,yl,star_x,star_y,m)
            star_corr[j]=star_corr_one
            if inside==1:
                Shears['clean'][i]=1
            break
    if plot_reg==True:
        plot_region_maskstar("remove_star.reg",star_corr)  ##if true, creat a ds9 region file

    Shears_remove=Shears[Shears['clean']==0]
        
    if os.path.isfile(mask_file):
        mask_obj = open(mask_file, 'rb')
        for mask in mask_obj:
            if mask[0:3] != 'box':
                continue
            else:
                
                mask_x = np.float(mask.split('(')[1].split(',')[0])
                mask_y = np.float(mask.split('(')[1].split(',')[1])
                mask_sizex = np.float(mask.split('(')[1].split(',')[2][:-1])
                mask_sizey = np.float(mask.split('(')[1].split(',')[3][:-1])
                
                mask_angle = np.float(mask.split('(')[1].split(',')[4][:-2])
                #rotate the shears into fram of reference of the mask
                shears_x_mask_ref = tools.ra_separation(Shears_remove['X_WORLD'], mask_y, mask_x , mask_y)
                shears_y_mask_ref = (Shears_remove['Y_WORLD'] - mask_y)*3600.
                
                shears_x_mask_rot = np.cos(mask_angle*np.pi/180.)*shears_x_mask_ref +  np.sin(mask_angle*np.pi/180.)*shears_y_mask_ref
                shears_y_mask_rot = -np.sin(mask_angle*np.pi/180.)*shears_x_mask_ref + np.cos(mask_angle*np.pi/180.)*shears_y_mask_ref
                
                inBox =  (shears_x_mask_rot < mask_sizex/2.) & \
                    (shears_x_mask_rot > -mask_sizex/2.) &\
                        (shears_y_mask_rot < mask_sizey/2.) &\
                            (shears_y_mask_rot > -mask_sizey/2.)
    
                Shears_remove = Shears_remove[ inBox == False ]


        fits.writeto(outFile,Shears_remove, clobber=True,output_verify='ignore' )

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

