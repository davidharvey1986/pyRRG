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
import RRGtools as at
import glob as glob
import pyfits as fits
import numpy as np
import numpy.lib.recfunctions as rec
import pickle as pkl
import os as os
def generateTrainingData(allGalaxyFiles=None, \
                             trainingDataPklFile='DataTrained.pkl'):
    '''
    Generate a table of data using the data in the file trainingData
    '''
    if allGalaxyFiles is None:

        allGalaxyFiles = glob.glob('trainingData/*uncor*')





         
    if os.path.isfile( trainingDataPklFile ):
        allTrainingData, allTrainingAnswers = \
          pkl.load(open(trainingDataPklFile,'rb'))
    else:
        allTrainingData, allTrainingAnswers = \
          filesToRecArray( allGalaxyFiles )
        pkl.dump([allTrainingData, allTrainingAnswers], \
                     open(trainingDataPklFile,'wb'))
    return  allTrainingData, allTrainingAnswers
    

def filesToRecArray( files ):
    '''
    Take a list of fits files and append them together in 
    to one rec-array
    
    They MUST have the same dtypes
    '''
    
    for i, iFile in enumerate(files):
        
        data, iStarGalClass = matchStarGalaxiesToData( iFile )
        
        if i==0:
            allDataArray = rec2array( data )

            iDataNoNan, starGalNoNan = \
              removeNans( allDataArray, iStarGalClass )

            
            starGalClass = starGalNoNan
            allData = iDataNoNan
        else:
            
            iFileData = rec2array( data )
            iDataNoNan, starGalNoNan = \
              removeNans( iFileData, iStarGalClass )

            allData = np.vstack( (allData, iDataNoNan))
            
            starGalClass = np.append(starGalClass, starGalNoNan)
        
        
    return allData, starGalClass
    
def generateTestData():
    '''
    This will pick one of the datasets
    in trainingData to test on
    
    '''



    testGalaxyFiles = glob.glob('TestData/*uncor.cat')

    featureLabels = getFeatureLabels( testGalaxyFiles[0] )
    
    testFeatures, testAnswers = \
      generateTrainingData(allGalaxyFiles=testGalaxyFiles, \
                             trainingDataPklFile='testData.pkl')

    return featureLabels, testFeatures, testAnswers

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


    includeNames = \
      ['MAG_AUTO','gal_size','MU_MAX','MAG_ISO','RADIUS','FLUX_AUTO',\
           'xxxx','yyyy','xyyy','xxyy','xx','xy','yy','e1','e2','prob',\
      'ell','skymed','exp_time','skysd']
   
    newArray = np.zeros((len(recArray),len(includeNames)), float)

    for i, iField in enumerate(includeNames):
        newArray[:,i] = recArray[iField]

    return newArray
        
def getFeatureLabels( fitsFile ):
    includeNames = fits.open(fitsFile)[1].data.columns.names
    #remove all those with err in it
    namesNoErr = [ i for i in includeNames if not 'err' in i ]
    #includeNames.remove('skymed')
    #includeNames.remove('exp_time')
    #includeNames.remove('skysw')
    #includeNames.remove('skysd')

    print namesNoErr
    #namesNoErr = \
    #    ['MAG_AUTO','gal_size','MU_MAX','MAG_ISO','RADIUS',\
    #        'xxxx','yyyy','xyyy','xxyy','ell']


      #,'skymed']
    namesNoErr = \
      ['MAG_AUTO','gal_size','MU_MAX','MAG_ISO','RADIUS','FLUX_AUTO',\
           'xxxx','yyyy','xyyy','xxyy','xx','xy','yy','e1','e2','prob',\
      'ell','skymed','exp_time','skysd']
    print namesNoErr

    return np.array(namesNoErr)


def matchStarGalaxiesToData( iFile ):
    '''
    I need to add a column to 'data' that has the classifcation
    of stars (0), galaxies(1) and neither (-1)
    
    '''
    cluster=iFile.split('_')[0]

    data = fits.open(iFile)[1].data
    classification = np.zeros(len(data))-1

    newIDcol = \
      [fits.Column(name='ID', format='D', array=np.arange(len(data)))]


    dataCols = data.columns + fits.ColDefs(newIDcol)
    dataWithID =  fits.BinTableHDU.from_columns(dataCols)
    
    dataWithID.writeto('DataID.fits', clobber=True)
    matchedGalaxyData = at.run_match(cluster+'_galaxies.fits',\
                                         'DataID.fits')[1].data

    
    classification[matchedGalaxyData['ID'].astype(int)] = 1

    matchedStarData = at.run_match(cluster+'_stars.fits',\
                                       'DataID.fits')[1].data
                                    
    classification[matchedStarData['ID'].astype(int)] = 0

    

    return data, classification

def removeNans( newArray, starGal ):

    #remove nan
    #newArray[ np.isfinite(newArray) == False ] = -99
    nanCheck = np.isfinite(np.sum(newArray, axis=1))
    newArrayNansRemoved = newArray[nanCheck, :]
    Nremoved = newArray.shape[0] - newArrayNansRemoved.shape[0]
    
    print("%i/%i removed due to nans" % (Nremoved, newArray.shape[0]))

    nanCheckField = np.isfinite(np.sum(newArray, axis=0))


    return newArrayNansRemoved, starGal[nanCheck]
