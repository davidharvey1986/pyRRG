'''
This script will train the data at the moment using a very simple
logisitic regression

'''
import os as os
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier as RF

from sklearn.svm import SVC
import generateTrainingData as gt
import pickle as pkl
import time
import ipdb as pdb

def trainDataSVM(nSamples=60000):
    pickleFileName = 'starGalaxyModelSVM.pkl'
    
    trainingFeatures, trainingAnswers = \
      gt.generateTrainingData()
 
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

    pkl.dump(fitClassifier, open(pickleFileName,'wb'))
    return pickleFileName
    
def trainDataRF(nTrees=100, retrain=True):
    '''
    This uses a random forest to classify the data
    '''
    pickleFileName = 'starGalaxyModelRF.pkl'

    if (os.path.isfile(pickleFileName)) & (not retrain):
        return pickleFileName
    trainingFeatures, trainingAnswers = \
      gt.generateTrainingData()
  

    print("Size of training data is %i with %i features" %\
              (trainingFeatures.shape[0], trainingFeatures.shape[1]))
              

    
    clf = RF( n_estimators=nTrees, criterion="gini", n_jobs=4)

    print("Training started at %s with %i Trees" % (time.ctime(),nTrees))
    startTime = time.time()
    fitClassifier = clf.fit( trainingFeatures, trainingAnswers)
    endTime = time.time()
    fitTime= endTime - startTime
    print("Time to fit classifier is %0.2f seconds" % fitTime)
    pkl.dump(fitClassifier, open(pickleFileName,'wb'), 2)
    
    
    return pickleFileName

