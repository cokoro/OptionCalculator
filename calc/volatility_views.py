from django.shortcuts import render
from django.http import HttpResponse

import numpy as np
from math import *
from scipy.stats import norm
from decimal import Decimal

def sigmacal(S0,r,q,T,K,premium,option_type):
    count = 0
    sigma = sqrt( 2 * abs((log( S0 / K ) + ( r - q ) * (T))/(T)))
    for i in range(1,1000):
        count = count+1
        d1 = ( log( S0 / K )+(r - q)*( T ))/( sigma * sqrt(T)) + 0.5 * sigma * sqrt(T)
        d2 = ( log( S0 / K )+(r - q)*( T ))/( sigma * sqrt(T)) - 0.5 * sigma * sqrt(T)
        vega = S0 * exp( -q * ( T )) * norm.pdf( d1 ) * sqrt( T )
        if option_type == 'C':
            calculatePrice = S0 * exp( -q * ( T )) * norm.cdf(d1) - K * exp( -r *(T))*norm.cdf(d2)
        else:
            calculatePrice = K * exp( -r * ( T )) * norm.cdf(-d2) - S0 * exp( -q *(T))*norm.cdf(-d1)
        
        #X(n+1) = X(n) - f(X(n))/f'(X(n))
        sigma = sigma - ( calculatePrice- premium ) / vega
    
        if abs(calculatePrice-premium)<1e-10:
            print(calculatePrice)
            break
    print(count)
    return round(sigma,4)


def volatility_index(request):
    return render(request, 'volatility.html')


def volatility(request, S, r, q, T, K, premium, call_put):
    c = sigmacal(int(S), float(r), float(q), float(T), float(K),float(premium), str(call_put))
    r = HttpResponse( str(c))
    return r
