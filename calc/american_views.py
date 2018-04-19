from django.shortcuts import render
from django.http import HttpResponse

from math import *
from scipy.stats import norm
import numpy as np


def BinomialTree(S0,sigma,r,T,K,N,option_type):
    
        
    deltaT = float(T) / N
 
    u = np.exp(sigma * np.sqrt(deltaT))
    d = 1.0 / u
    
    vector =  np.asarray([0.0 for i in range(N+1)])
    vector1 = np.asarray([(S0 * u**j * d**(N - j)) for j in range(N+1)])
    vector2 =np.asarray( [float(K) for i in range(N+1)])
    
 
    num = np.exp(r * deltaT)
    p = (num - d)/ (u - d)
    q = 1.0 - p
 
   
    if option_type =="C":
        vector[:] = np.maximum(vector1-vector2, 0.0)
    elif option_type =="P":
        vector[:] = np.maximum(-vector1+vector2, 0.0)
    
   
    #calculate backward the option prices
    for i in range(N-1, -1, -1):
       vector[:-1]=np.exp(-r * deltaT) * (p * vector[1:] + q * vector[:-1])
       vector1[:]=vector1[:]*u
       if N>0: 
        if option_type =="C":
                vector[:]=np.maximum(vector[:],vector1[:]-vector2[:])
        elif option_type=="P":
                vector[:]=np.maximum(vector[:],-vector1[:]+vector2[:])
                
    # print fs
    return vector[0]

def american_index(request):
    return render(request, 'american.html')

def american(request,S, sigma, r, T, K, N, call_put):
    c = BinomialTree(int(S), float(sigma), float(r), float(T), float(K), int(N), str(call_put))
    r = HttpResponse(str(c))
    return r
