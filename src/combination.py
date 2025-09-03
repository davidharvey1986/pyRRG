import numpy as np
from scipy.special import factorial
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
    
    return factorial(n)/(factorial(k)*factorial(n-k))

