from django.shortcuts import render
from django.http import HttpResponse

from math import *
from scipy.stats import norm
import numpy as np
import random


def geo_asian_option(S0,sigma,r,T,K,step,option_type):
    N = float(step)
    
    sigmaHat = sigma * sqrt((N+1)*(2*N+1)/(6*pow(N,2)))
    muHat = (r-0.5*pow(sigma, 2))*(N+1)/(2*N)+0.5*pow(sigmaHat, 2)
    
    d1 = ( log( S0 / K )+ (muHat + 0.5 * pow(sigmaHat,2)) * (T)) / (sigmaHat * sqrt(T))
    d2 = d1 - sigmaHat*sqrt(T)

    if option_type == 'C':
        return round(exp(-r*(T))*(S0 * exp(muHat*(T))*norm.cdf(d1) - K*norm.cdf(d2)),4)
    elif option_type == 'P':
        return round(exp(-r*(T))*(-S0 * exp(muHat*(T))*norm.cdf(-d1) + K*norm.cdf(-d2)),4)

def g_asian_index(request):
    return render(request, 'g_asian.html')


def g_asian(request, S, sigma, r, T, K, step, call_put):
    c = geo_asian_option(int(S), float(sigma), float(r), float(T), float(K), float(step), str(call_put))
    r = HttpResponse( str(c))
    return r
