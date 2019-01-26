
import pickle as pkl
import generateTrainingData as gtd
from matplotlib import pyplot as plt
import ipdb as pdb
import trainData as td
def testClassifier(retrain=False):
    '''
    This will generate some test data and plot the stars and galaxies
    on the standard plot
    '''
    if retrain:
        td.trainData()
    featureLabels, testFeatures, testAnswers = \
      gtd.generateTestData()

    mlClassifier = pkl.load(open('starGalaxyModel.pkl'))

    predictClassification = mlClassifier.predict(testFeatures)
    print mlClassifier.score(testFeatures, testAnswers)
    print("Number of objects in test data is %i"% len(testFeatures))

    plt.plot( testFeatures[predictClassification==1,\
                               featureLabels=='MAG_AUTO'], \
                  testFeatures[predictClassification==1,\
                                   featureLabels=='MU_MAX'],'r.')
    plt.plot( testFeatures[predictClassification==0,\
                            featureLabels=='MAG_AUTO'], \
                  testFeatures[predictClassification==0,\
                                   featureLabels=='MU_MAX'],'y.')
    plt.show()
