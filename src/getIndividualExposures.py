'''
A script that will return the individual expsoures 
that made up the image

If some exist but not all, it will complain

If none of them exists it will state so and then
ask if we should use the input image as a way to get the
PSF

'''
from astropy.io import fits
import numpy as np
import os

def getIndividualExposures( inputFileName ):
    '''
    From an input file (named inputFIleName) return
    a list of names that are the input files that
    made up that image
    '''


    inputHeader = fits.open(inputFileName)[0].header

    exposureNameList = \
      np.unique([ inputHeader[i].split('_')[0]+'_drz_sci.fits'\
                      for i in inputHeader.keys() \
                      if 'DATA' in i ])

    fileCheck = []
    for iFile in exposureNameList:
        if not os.path.isfile(iFile):
            print("%s file not found" % iFile )
        else :
            print("%s found!" % iFile )
        fileCheck.append( os.path.isfile(iFile) )
    
    if np.all( np.array(fileCheck) == False):

        raise ValueError("No individual files found")

    if np.all( np.array(fileCheck) ) == False:
        contninueWithout = input("Haven't found all files, continue anyway? (y or n)")
        if contninueWithout == 'y':
            return exposureNameList[ fileCheck ]
        else:
            raise ValueError("Not all individual files found")

    return exposureNameList
                
            
                
        
