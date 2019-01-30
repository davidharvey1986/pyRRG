'''
This script will train the data at the moment using a very simple
logisitic regression

'''

from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
import generateTrainingData as gt
import pickle as pkl

def trainData():

    trainingFeatures, trainingAnswers = \
      gt.generateTrainingData()
    standardGamma = 1./(trainingFeatures.shape[1] * trainingFeatures.std())
    print("Size of training data is %i with %i features" %\
              (trainingFeatures.shape[0], trainingFeatures.shape[1]))
    clf = SVC(C=100, cache_size=200, class_weight='balanced', coef0=0.0,
        decision_function_shape='ovr', degree=3, gamma=0.0001, kernel='linear',
        max_iter=-1, probability=False, random_state=None, shrinking=True,
        tol=0.001, verbose=False)
    
    fitClassifier = clf.fit( trainingFeatures, trainingAnswers)
    
    pkl.dump(fitClassifier, open('starGalaxyModel.pkl','wb'))

    
