import pyfits as fits
import numpy as np
from scipy import spatial
import matplotlib.pyplot as plt

def remove_object(rrg_catalogue, output_catalogue, FWHM_to_radius=1):
    hdulist=fits.open(rrg_catalogue)
    data_org=hdulist[1].data
    print("num of objects in the rrg catalogue:",len(data_org))

    x=np.array(data_org['X_IMAGE'])
    y=np.array(data_org['Y_IMAGE'])
    FWHM=np.array(data_org['FWHM_IMAGE'])*FWHM_to_radius
    mask=np.ones(len(x),dtype=bool)
    ori_order=np.arange(len(x))                     ##original index in rrg_catalogue

    temp=np.vstack((x,y,FWHM,ori_order))            ##temp: x, y, FWHM, ori_order
    temp=temp.T
    temp_sort=temp[temp[:,2].argsort()[::-1]]       ##sorting by FWHM, from big objects to small one
    pos=list(zip(temp_sort[:,0],temp_sort[:,1]))          ##x, y
    tree=spatial.KDTree(pos)

    for i in range(len(x)):
        if mask[i]==False:
            continue
        rm_obj=tree.query_ball_point(pos[i],2.0*temp_sort[i][2]) #find neighbour objects
        for index in rm_obj:
            if index==i:
                continue
            delta_x=temp_sort[index][0]-temp_sort[i][0]
            delta_y=temp_sort[index][1]-temp_sort[i][1]
            distance=np.sqrt((delta_x)**2+(delta_y)**2)
            if distance<(temp_sort[i][2]+temp_sort[index][2]): #remove overlapped FWHM ojects
                mask[index]=False

    temp_sort=temp_sort[mask]
    sort_data=temp_sort[temp_sort[:, 3].argsort()]  ##sorting by the original index
    data_org=data_org[sort_data[:, 3].astype(int)]

    print(("Num of objects after removing double-detection: %i" % len(data_org)))
        
    fits.writeto(output_catalogue, data=data_org, \
                     header=hdulist[1].header, \
                     clobber=True,output_verify='ignore')

