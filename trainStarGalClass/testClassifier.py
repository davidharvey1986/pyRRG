
import pickle as pkl
import generateTrainingData as gtd
from matplotlib import pyplot as plt
import ipdb as pdb
import trainData as td
import numpy as np

def testClassifier(retrain=False):
    '''
    This will generate some test data and plot the stars and galaxies
    on the standard plot
    '''

        
    modelName = td.trainDataRF(retrain=retrain)
        
    featureLabels, testFeatures, testAnswers = gtd.generateTestData()

    mlClassifier = pkl.load(open(modelName))

    predictClassification = mlClassifier.predict(testFeatures)

    
    plt.figure(1)
    classifierScore( testAnswers, predictClassification)
    print("Classifier score is %0.2f" %mlClassifier.score(testFeatures,testAnswers))
    plt.plot( testFeatures[predictClassification==1,\
                               featureLabels=='MAG_AUTO'], \
                  testFeatures[predictClassification==1,\
                                   featureLabels=='MU_MAX'],'r.',\
                  label='Galaxies')

    featureLabels, testFeatures, testAnswers = \
      gtd.generateTestData()

    

    plt.plot( testFeatures[predictClassification==1,\
                               featureLabels=='MAG_AUTO'], \
                  testFeatures[predictClassification==1,\
                                   featureLabels=='MU_MAX'],'r.')
    plt.plot( testFeatures[predictClassification==0,\
                            featureLabels=='MAG_AUTO'], \
                  testFeatures[predictClassification==0,\
                                   featureLabels=='MU_MAX'],'y.',\
                  label='Stars')
    plt.plot( testFeatures[predictClassification==-1,\
                            featureLabels=='MAG_AUTO'], \
                  testFeatures[predictClassification==-1,\
                                   featureLabels=='MU_MAX'],'g.', \
                  label='Noise')
    plt.xlim(15,30)
    plt.xlabel('Magnitdue')
    plt.ylabel('Mu Max')
    plt.savefig('testRandomForest.pdf')
    plt.show()


def classifierScore( answer, predict):
    Nstars = len(answer[answer==0])
    Ngalaxies = len(answer[answer==1])
    Nnoise = len(answer[answer==-1])

    NstarCorrect = len(predict[answer==0][predict[answer==0] == 0])
    NgalCorrect = len(predict[answer==1][predict[answer==1] == 1])
    NnoiseCorrect = len(predict[answer==-1][predict[answer==-1] == -1])
    print("Number of correct stars: %i/%i" % (NstarCorrect,Nstars))
    print("Number of correct galaxies: %i/%i"% (NgalCorrect,Ngalaxies))
    print("Number of correct noise: %i/%i"% (NnoiseCorrect, Nnoise))


