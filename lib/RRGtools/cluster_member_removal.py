'''
Remove cluster members from catlogue

This algorithm will take 2 filters from HST and match the catalogues
It will then look and the color-magnitude plot and find the red sequence.
From this it wil remove galaxies that appear to be in the cluster from the catalogue


It requires that the pyRRG code has been run on two different bands of the same cluster
'''

import pyfits as fits
import RRGtools as tools
from matplotlib import pyplot as plt
import numpy as np
from numpy.lib.recfunctions import append_fields as append_rec
import os as os
def clusterMemberRemove( red_band, blue_band, outname='cluster_mem_rem.fits'):
    '''
    INPUTS:
        - red_band is the name of the observation in the redder band
        - blue_band is the name of the observation in the bluer band

    OUTPUTS
        - a joint catalogue with the two magnitudes, but with the cluster members 
        removed.

    NOTE: It changes the name of the magnitudes within the structure, so the 
        name of the file is assumed to be ${OBJECTNAME}_${FILTER}_drz_sci.fits
    '''
    if not ( os.path.isfile( red_band ) & \
                 os.path.isfile( blue_band ) ):
        raise ValueError("File not found")
    
    
    joint_catalogue = tools.run_match( red_band, blue_band)

    #This assumes the name of the file is in the order  ${OBJECTNAME}_${FILTER}_drz_sci.fits
    red_band_filter = red_band.split('_')[1]
    blue_band_filter = blue_band.split('_')[1]
    joint_catalogue[1].data = append_rec( joint_catalogue[1].data, \
                                              red_band_filter, \
                                              joint_catalogue[1].data['MAG_AUTO_1'], \
                                              usemask=False, asrecarray=True)
    joint_catalogue[1].data = append_rec( joint_catalogue[1].data,\
                                              blue_band_filter, \
                                              joint_catalogue[1].data['MAG_AUTO_2'], \
                                              usemask=False, asrecarray=True)

    color =  joint_catalogue[1].data[red_band_filter] - \
      joint_catalogue[1].data[blue_band_filter]
    plt.plot( joint_catalogue[1].data[red_band_filter], color, 'b*')
    plt.xlim(10,30)
    plt.ylim(-10,10)
    plt.show(block=False)


    while True:
        UpperThreshold = np.float(input('Please input the upper threshold of the red-sequence: '))
        LowerThreshold = np.float(input('Please input the lower threshold of the red-sequence: '))
        MagThreshold = np.float(input('Please input the magnitude threshold of the red-sequence: '))

        clusterMembers =  (color < UpperThreshold) & \
          (color > LowerThreshold) & \
          (joint_catalogue[1].data[red_band_filter] < MagThreshold)
        plt.plot( joint_catalogue[1].data[red_band_filter], color, 'b*')

        plt.plot(  joint_catalogue[1].data[red_band_filter][clusterMembers], \
                       color[clusterMembers], 'r*')
        plt.ylim(np.min([-2,LowerThreshold*2.]), np.max([2,UpperThreshold*2.]))
        plt.xlim(10, np.max([MagThreshold*1.5,28]))
        plt.draw()
        done = input('Are you happy with this? (Yes or No): ')
        if done == 'Yes':
            break
    plt.close()

    final_data = joint_catalogue[1].data[ clusterMembers == False ]
    final_cat_names = np.array(joint_catalogue[1].data.columns.names)
    joint_columns = []
    
    for i in range(len(final_cat_names)):
        if '_1' in final_cat_names[i]:
            final_cat_name = '_'.join(final_cat_names[i].split('_')[:-1])
        else:
            final_cat_name = final_cat_names[i]
        joint_columns.append( \
                fits.Column(name=final_cat_name, \
                            array=final_data[final_cat_names[i]], \
                            format=final_data[final_cat_names[i]].dtype))
                            
    final_cat = fits.BinTableHDU.from_columns(joint_columns)
   
    final_cat.writeto( outname, clobber=True)
    return fits.open(outname)
                       
                                  
    
    
