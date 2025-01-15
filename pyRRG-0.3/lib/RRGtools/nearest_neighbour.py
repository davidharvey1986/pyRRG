import numpy as np
import RRGtools as tools

def nearest_neighbour(  x, y):

    '''

    Take in a set of objects with their positions
    as described by x and y, and find the distance
    to the nearest neighbour

    INPTUS : X, Y : TWO VECTORS OF EQUAL LENGTH THAT
                    ARE THE POSITIONS OF EACH OBJECT
                    MUST BE A FLOAT


    OUTPUTS : A VECTOR OF LEN(X) LONG THAT GIVES THE DISTANCE
              (IN UNITS OF THE COORDINATES) TO THE NEAREST
              NEIGHBOUR

    '''

    x_arr = np.array(np.transpose(np.matrix(x))*np.ones(len(x)))
    y_arr = np.array(np.transpose(np.matrix(y))*np.ones(len(y)))

    x_arr1 = np.transpose(x_arr)
    y_arr1 = np.transpose(y_arr)

    radius = np.sqrt( (x_arr - x_arr1)**2 +(y_arr - y_arr1)**2)

   

    radius[ np.arange(len(x)), np.arange(len(x)) ] = np.max(radius)

    indexes = np.argmin(radius, axis=1)
    
    result = np.min(radius, axis=1)
        
    return result, indexes

    
def match_cats( catA_x, catA_y, catB_x, catB_y, cut=None, wcs=False):
    '''
    Simialr to nearest neighbour except instead
    of nearest to self, nearest to a nother cat

    keywords : cut : place a distane cut on the matching, if the distance is
                greater thant this to the nearest return -1 as th eindex

                wcs : if true the x and y are in wcs and therefore must use
                ra_separation

    '''

    x_arrA = np.array(np.transpose(np.matrix(catA_x))*np.ones(len(catB_x)))
    y_arrA = np.array(np.transpose(np.matrix(catA_y))*np.ones(len(catB_y)))

    x_arrB = np.transpose(np.array(np.transpose(np.matrix(catB_x))*np.ones(len(catA_x))))
    y_arrB = np.transpose(np.array(np.transpose(np.matrix(catB_y))*np.ones(len(catA_y))))

    if not wcs:
        radius = np.sqrt( (x_arrA - x_arrB)**2 +(y_arrA - y_arrB)**2)
    else:
        radius = tools.ra_separation( x_arrA, y_arrA, x_arrB, y_arrB, abs=True)

    indexes = np.argmin(radius, axis=1)
    
    result = np.min(radius, axis=1)

    if cut is not None:
        indexes[ result > cut ] = -1
        
    return result, indexes
