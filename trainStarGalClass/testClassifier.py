
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
    print  featureLabels.shape
    mlClassifier = pkl.load(open('starGalaxyModel.pkl'))
    print testFeatures.shape
    predictClassification = mlClassifier.predict(testFeatures)
    print plt.plot(predictClassification-testAnswers)
    plt.figure(2)
    print mlClassifier.score(testFeatures, testAnswers)
    plt.plot( testFeatures[predictClassification==1,\
                               featureLabels=='MAG_AUTO'], \
                  testFeatures[predictClassification==1,\
                                   featureLabels=='gal_size'],'r.')
    plt.plot( testFeatures[predictClassification==0,\
                            featureLabels=='MAG_AUTO'], \
                  testFeatures[predictClassification==0,\
                                   featureLabels=='gal_size'],'y.')
    plt.show()
