from __future__ import division
import numpy as np
from math import sqrt,exp
import math
from scipy.stats import norm

from django.shortcuts import render
from django.http import HttpResponse

import scipy
number= 10000

def Geometric_Basket_Option(S1, S2, sigma_1, sigma_2, r, T, K, rau, option_type, flag):
    Bg_zero = sqrt(S1 * S2)
    Bg_sigma = sqrt(sigma_1 * sigma_1 + sigma_2 * sigma_2 + 2 * sigma_1 * sigma_2 * rau) / 2
    Bg_miu = r - 1/2 * (sigma_1 * sigma_1 + sigma_2 * sigma_2) / 2 + 1/2 * Bg_sigma*Bg_sigma
    d_1 = (math.log(Bg_zero / K) + (Bg_miu + 1/2 * Bg_sigma*Bg_sigma)*T) / (Bg_sigma*sqrt(T))
    d_2 = (math.log(Bg_zero / K) + (Bg_miu - 1/2 * Bg_sigma*Bg_sigma)*T) / (Bg_sigma*sqrt(T))
    if(flag=="Yes"): #closed form
        if (option_type == "C"):
            price = exp(-r * T) * (Bg_zero * exp(Bg_miu * T) * norm.cdf(d_1,0,1) - K * norm.cdf(d_2,0,1))
        else:
            price = exp(-r * T) *(K * norm.cdf(-d_2,0,1) - Bg_zero * exp(Bg_miu * T) * norm.cdf(-d_1,0,1))
        return round(price, 8)
    else: # Standard monte carlo
        Dt = 1e-2
        N = int(T / Dt)
        geoPayoff=np.zeros(number)
        drift_1 = exp((r - 0.5 * sigma_1 * sigma_1) * Dt)
        drift_2 = exp((r - 0.5 * sigma_2 * sigma_2) * Dt)
        Spath_1=np.zeros(N)
        Spath_2=np.zeros(N)
        scipy.random.seed(1)
        for i in range(0,number):
            Rand1 = scipy.random.randn(N-1)
            Rand2 = scipy.random.randn(N-1)
            Spath_1[0] = S1
            Spath_2[0] = S2
            for j in range(1,N):
                growthFactor_1 = drift_1 * exp(sigma_1 * sqrt(Dt) * Rand1[j-1])
                growthFactor_2 = drift_2 * exp(sigma_2 * sqrt(Dt) * (rau * Rand1[j-1] + sqrt(1 - rau * rau) * Rand2[j-1]))
                Spath_1[j] = Spath_1[j - 1] * growthFactor_1
                Spath_2[j] = Spath_2[j - 1] * growthFactor_2
            geomean = sqrt(Spath_1[-1] * Spath_2[-1])
            if (option_type == "C"):
                geoPayoff[i] = exp(-r * T) * max(geomean - K, 0)
            else:
                geoPayoff[i] = exp(-r * T) * max(K - geomean, 0)
        Pmean = np.mean(geoPayoff)
        Pstd = np.std(geoPayoff)
        confmc = [Pmean - 1.96 * Pstd / sqrt(number), Pmean + 1.96 * Pstd / sqrt(number)]
    return round(np.mean(confmc),4)

def g_basket_index(request):
    return render(request, 'g_basket.html')

def G_Basket(request,S1, S2, sigma_1, sigma_2, r, T, K, rau, call_put, flag):
    c = Geometric_Basket_Option(int(S1), int(S2), float(sigma_1), float(sigma_2), float(r), float(T), float(K), float(rau), str(call_put), str(flag))
    r = HttpResponse(str(c))
    return r
