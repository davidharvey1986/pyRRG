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
import json

def getIndividualExposures(  **kwargs ):
    '''
    From an input file (named inputFIleName) return
    a list of names that are the input files that
    made up that image
    '''
    if 'FILENAME' not in kwargs.keys():
        raise ValueError("Please state FILENAME in kwargs")

    if not kwargs['no_exposures']:
        if kwargs['exposureNameList'] is None:
            if kwargs['telescope'] == 'JWST':
                exposureNameList = fits.open(kwargs['field'])[8].data['FILENAME']
            else:
                inputHeader = fits.open(kwargs['field'])[0].header
                exposureNameList = \
                    np.unique([ inputHeader[i].split('_')[0]+'_drz_sci.fits'\
                      for i in inputHeader.keys() \
                      if 'DATA' in i ])
        else:
            exposureNameList = np.loadtxt(kwargs['exposureNameList'], dtype=object)
    else:
        print("WARNING - NO EXPOSURES - USING DRIZZLE FILE TO ESTIMATE")
        exposureNameList = kwargs['field' ]  
                    
    fileCheck = []
    for iFile in  np.atleast_1d(exposureNameList):
        if not os.path.isfile(iFile):
            print("%s file not found" % iFile )
        fileCheck.append( os.path.isfile(iFile) )
    
    if np.all( np.array(fileCheck) == False):

        raise ValueError("No individual files found - if no files exist use --no-exposures flag")

    if np.all( np.array(fileCheck) ) == False:
        contninueWithout = input("Haven't found all files, continue anyway? (y or n)")
        if contninueWithout == 'y':
            return exposureNameList[ fileCheck ]
        else:
            raise ValueError("Not all individual files found")

    return np.atleast_1d(exposureNameList)
                
            
                
        
