
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

    mlClassifier = pkl.load( open(modelName,'rb'), encoding='latin1')

    predictClassification = mlClassifier.predict(testFeatures)

    

    classifierScore( testAnswers, predictClassification)
    print(("Classifier score is %0.2f" %mlClassifier.score(testFeatures,testAnswers)))
    #plt.plot( testFeatures[predictClassification==1,\
    #                           featureLabels=='MAG_AUTO'], \
    #              testFeatures[predictClassification==1,\
    ##                               featureLabels=='MU_MAX'],'r.',\
     #             label='Galaxies')

    featureLabels, testFeatures, testAnswers = \
      gtd.generateTestData()

    


    fig = plt.figure( )

    axis = plt.gca()

    cleanFeatures, cleanImportance = \
      cleanLabelsAndImportance( featureLabels, \
                    mlClassifier.feature_importances_)

    
    axis.bar(np.arange(len(cleanFeatures)), \
            cleanImportance, edgecolor='black', color='grey')
    axis.set_xticks(np.arange(len(cleanFeatures)))
    axis.set_xticklabels(cleanFeatures, rotation='vertical')
    axis.set_ylabel('Relative Importance')
    axis.set_xlabel('Feature', labelpad=-20)

    left, bottom, width, height = [0.5, 0.5, 0.4, 0.4]
    ax2 = fig.add_axes([left, bottom, width, height])
    ax2.plot( testFeatures[predictClassification==1,\
                             featureLabels=='MAG_AUTO'], \
                  testFeatures[predictClassification==1,\
                                   featureLabels=='MU_MAX'],'r.', \
                  ms=1, label='Galaxies')
    ax2.plot( testFeatures[predictClassification==0,\
                            featureLabels=='MAG_AUTO'], \
                  testFeatures[predictClassification==0,\
                                   featureLabels=='MU_MAX'],'y.',\
                  label='Stars', ms=1)
    ax2.plot( testFeatures[predictClassification==-1,\
                            featureLabels=='MAG_AUTO'], \
                  testFeatures[predictClassification==-1,\
                                   featureLabels=='MU_MAX'],'g.', \
                  label='Noise', ms=1)
    ax2.set_xlim(15,30)
    ax2.set_xlabel('Magnitude')
    ax2.set_ylabel(r'$\mu_{\rm max}$')
    ax2.set_title('GAL-0364-52000-084')
    ax2.legend()

    fig.tight_layout()
    plt.savefig('classifierImportance.pdf')
    plt.show()

def cleanLabelsAndImportance( labels, importance):
    
    cleanLabels = []
    cleanImportance = []
    alreadyMag = False
    
    for iImportanceIndex, iLabel in enumerate(labels):
        
        if 'MAG' in iLabel:
            print(iLabel)
            if  alreadyMag ==False:
                cleanLabels.append('Magnitude')
                cleanImportance.append(importance[iImportanceIndex])
                alreadyMag = True
                continue
            else:
                continue
        if ('FLUX' in iLabel) | ('RADIUS' in iLabel) | \
          ('prob' in iLabel):
            continue
        if '_' in iLabel:
            cleanLabels.append(' '.join(iLabel.split('_')))
            cleanImportance.append(importance[iImportanceIndex])
        
        elif ('xx' in iLabel) | ('xy' in iLabel)| ('yy' in iLabel):
            cleanLabels.append(r'$J_{'+iLabel+'}$')
            cleanImportance.append(importance[iImportanceIndex])
        else:
            cleanLabels.append(iLabel)
            cleanImportance.append(importance[iImportanceIndex])
        
        

    return cleanLabels, cleanImportance
        
    
def classifierScore( answer, predict):
    Nstars = len(answer[answer==0])
    Ngalaxies = len(answer[answer==1])
    Nnoise = len(answer[answer==-1])

    NstarCorrect = len(predict[answer==0][predict[answer==0] == 0])
    NgalCorrect = len(predict[answer==1][predict[answer==1] == 1])
    NnoiseCorrect = len(predict[answer==-1][predict[answer==-1] == -1])
    print(("Number of correct stars: %0.2f" % (NstarCorrect /Nstars)))
    print(("Number of correct galaxies: %0.2f"% (NgalCorrect / Ngalaxies)))
    print(("Number of correct noise: %0.2f"% (NnoiseCorrect / Nnoise)))


if __name__ == '__main__':
    testClassifier()
