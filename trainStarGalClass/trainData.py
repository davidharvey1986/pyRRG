'''
This script will train the data at the moment using a very simple
logisitic regression

'''

from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
import generateTrainingData as gt
import pickle as pkl
import time
import ipdb as pdb
def trainData(nSamples=60000):

    trainingFeatures, trainingAnswers = \
      gt.generateTrainingData()
    nStart = 17040
    #12941 - 19041 seems to be an important sample
    iStart =  0 #nStart+nSamples
    
    #print trainingFeatures[nStart+nSamples+1,5:]
    #print trainingFeatures[nStart+nSamples,5:]
    #nFeat = 6
    
    #trainingFeatures = trainingFeatures[iStart:nSamples+iStart,:]
    #trainingAnswers = trainingAnswers[iStart:nSamples+iStart]
    standardGamma = 1./(trainingFeatures.shape[1] * trainingFeatures.std())
    print("Size of training data is %i with %i features" %\
              (trainingFeatures.shape[0], trainingFeatures.shape[1]))
              
    clf = SVC(C=100, cache_size=200, class_weight='balanced', coef0=0.0,
        decision_function_shape='ovr', degree=3, \
        gamma=0.0001, kernel='linear',
        max_iter=-1, probability=False, \
        random_state=None, shrinking=True,
        tol=0.001, verbose=False)
    print("Training started at %s" % time.ctime())
    startTime = time.time()
    fitClassifier = clf.fit( trainingFeatures, trainingAnswers)
    endTime = time.time()
    fitTime= endTime - startTime
    print("Time to fit classifier is %0.2f seconds" % fitTime)
    pkl.dump(fitClassifier, open('starGalaxyModel.pkl','wb'))

    
