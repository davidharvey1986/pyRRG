'''
Trying to find the best parameters for the SVM

Going to use this for now.
https://stackoverflow.com/questions/46330329/finding-the-values-of-c-and-gamma-to-optimise-svm

But i need to research some more


'''
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
import generateTrainingData as gt
import pickle as pkl
import os
def gridSearch():
    '''
    Using the grid search technique automatically with the 
    classifier

    '''
    trainingFeatures, trainingAnswers = \
      gt.generateTrainingData()
      
    params_grid = {'C': [0.001, 0.01, 0.1, 1, 10, 100],
          'gamma': [0.0001, 0.001, 0.01, 0.1],
          'kernel':['linear','rbf','poly'] }

    
    #Create the GridSearchCV object
    pickleFile = 'grid_clf.pkl'
    if not os.path.isfile( pickleFile ):
        grid_clf = GridSearchCV(SVC(class_weight='balanced'), \
                                params_grid,n_jobs=4)
        pkl.dump(grid_clf,open(pickleFile,'wb'))
    else:
        grid_clf = pkl.load(open(pickleFile,'rb'))
        

    #Fit the data with the best possible parameters
    print("Finding best fit parameters")
    grid_clf = grid_clf.fit(trainingFeatures, trainingAnswers)
    
    pkl.dump(grid_clf,open("bestSVMparams.pkl",'wb'))
    #Print the best estimator with it's parameter
    print grid_clf.best_estimator
