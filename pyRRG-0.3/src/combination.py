import numpy as np

def combination( n, k):

    '''
    ;PURPOSE : TO FIND THE BINOMIAL COEEFICIENT OF SET OF N 
    ;          ELEMENTS AND K COMINATIONS
    
    ;INPUTS : 
    ;     N : the number in the set
    ;     K : the number of combinations
    
    
    ;RETURNS : 
    ;    ( N )  =     N!
    ;    ( K )     ------
    ;              K!(N-K)!
    '''
    
    return np.math.factorial(n)/(np.math.factorial(k)*np.math.factorial(n-k))

