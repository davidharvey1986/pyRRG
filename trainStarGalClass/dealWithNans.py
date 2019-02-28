'''
What am i going to do with the nans in the data?
as they seem to have some okay galaxies,


'''
import pickle as pkl
from matplotlib import pyplot as plt
import numpy as np
#First I am going to look at what the nans come from.
def main():
    includeNames = \
      ['MAG_AUTO','gal_size','MU_MAX','MAG_ISO','RADIUS','FLUX_AUTO',\
           'xxxx','yyyy','xyyy','xxyy','xx','xy','yy','e1','e2','prob',\
      'ell','skymed','exp_time','skysd']
    
    allTrainingData, allTrainingAnswers =  \
      pkl.load(open('DataTrained.pkl','rb'))

    counter = np.ones(allTrainingData.shape[0])
    print("TOTAL COUNTER IS %i" % np.sum(counter))
    nanVector = [ np.sum(counter[(allTrainingData[:,i] == -99) |\
                             (allTrainingData[:,i] == 0)]) \
                      for i in xrange(allTrainingData.shape[1] )]
    plt.plot(np.arange(len(includeNames)), nanVector,'*')

    plt.xticks(np.arange(len(includeNames)), includeNames,rotation=90)
    plt.show()
