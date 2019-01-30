'''
This program is going to generate the trainig data.

It needs a load of images at varying depth, of which the moments
have been measured.

Form this I will generate a table of features and then the final column, 
which will be galaxy (1) or star (0)

For now I will just use the SLACS images that i have downloaded.
But in the future i will require varying depth

'''
import ipdb as pdb

import glob as glob
import pyfits as fits
import numpy as np

def generateTrainingData():
    '''
    Generate a table of data using the data in the file trainingData
    '''

    trainingGalaxyFiles = glob.glob('trainingData/*galaxies*')[1:]
    trainingStarsFiles = glob.glob('trainingData/*stars*')[1:]

    trainingGalaxy = filesToRecArray( trainingGalaxyFiles )
    trainingStars = filesToRecArray( trainingStarsFiles )
    
    allTrainingData = np.vstack((trainingGalaxy,trainingStars))
    allTrainingAnswers = np.append(np.ones(len(trainingGalaxy)), \
                                       np.zeros(len(trainingStars)))
    
    return allTrainingData, allTrainingAnswers
    

def filesToRecArray( files ):
    '''
    Take a list of fits files and append them together in 
    to one rec-array
    
    They MUST have the same dtypes
    '''
    for i, iFile in enumerate(files):
        data = fits.open(iFile)[1].data
        print iFile
        if i==0:
            allData = rec2array( data )
        else:
            iFileData = rec2array( data )
            allData = np.vstack( (allData, iFileData))
    
    return allData
    
def generateTestData():
    '''
    This will pick one of the datasets
    in trainingData to test on
    
    '''


    trainingGalaxyFiles = [glob.glob('trainingData/*galaxies*')[0]]
    trainingStarsFiles = [glob.glob('trainingData/*stars*')[0]]
    print trainingStarsFiles
    #trainingGalaxyFiles = [glob.glob('trainingData/clusters/*galaxies*')[-1]]
    #trainingStarsFiles = [glob.glob('trainingData/clusters/*stars*')[-1]]
    
    
    trainingGalaxy = filesToRecArray( trainingGalaxyFiles )
    trainingStars = filesToRecArray( trainingStarsFiles )

    allTrainingData = np.vstack((trainingGalaxy,trainingStars))
    allTrainingAnswers = np.append(np.ones(len(trainingGalaxy)), \
                                       np.zeros(len(trainingStars)))
    featureLabels = getFeatureLabels(trainingGalaxyFiles[0])
    return featureLabels, allTrainingData, allTrainingAnswers

def rec2array( recArray):
    '''
    the sklearn classifier requires a normal array
    of homogenous dtype, so i need to convert

    '''

    #dont include the errors in this fit
    includeNames = [ i for i in list(recArray.columns.names) if not 'err' in i ]
    #includeNames.remove('skymed')
    #includeNames.remove('exp_time')
    #includeNames.remove('skysw')
    #includeNames.remove('skysd')

    includeNames=['MAG_AUTO','gal_size','MU_MAX','MAG_ISO','BACKGROUND']
    
    newArray = np.zeros((len(recArray),len(includeNames)), float)

    for i, iField in enumerate(includeNames):
        newArray[:,i] = recArray[iField]

    #remove nan
    #newArray[ np.isfinite(newArray) == False ] = -99
    nanCheck = np.isfinite(np.sum(newArray, axis=1))
    newArrayNansRemoved = newArray[nanCheck, :]
    Nremoved = newArray.shape[0] - newArrayNansRemoved.shape[0]
    
    print("%i/%i removed due to nans" % (Nremoved, newArray.shape[0]))

    nanCheckField = np.isfinite(np.sum(newArray, axis=0))
    
    return newArrayNansRemoved
        
def getFeatureLabels( fitsFile ):
    includeNames = fits.open(fitsFile)[1].data.columns.names
    #remove all those with err in it
    namesNoErr = [ i for i in includeNames if not 'err' in i ]
    #includeNames.remove('skymed')
    #includeNames.remove('exp_time')
    #includeNames.remove('skysw')
    #includeNames.remove('skysd')
    namesNoErr=['MAG_AUTO','gal_size','MU_MAX','MAG_ISO','RADIUS']
    print includeNames
    return np.array(namesNoErr)
