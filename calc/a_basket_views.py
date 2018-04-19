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

def Arithmetic_Basket_Option(S1, S2, sigma_1, sigma_2, r, T, K, rau, option_type, flag):
    dt = T
    drift1 = exp((r - 0.5 * sigma_1 * sigma_1) * dt)
    drift2 = exp((r - 0.5 * sigma_2 * sigma_2) * dt)
    arithPayOff = np.empty(number, dtype=float)
    geoPayoff = np.empty(number, dtype=float)
    scipy.random.seed(10)

    for i in range(0, number, 1):
        Rand1 = scipy.random.randn(1)
        Rand2 = scipy.random.randn(1)
        #print(Rand1,Rand2)
        growthFactor1 = drift1 * exp(sigma_1 * sqrt(dt) * Rand1)
        S1next = S1 * growthFactor1
        growthFactor2 = drift2 * exp(sigma_2 * sqrt(dt) * (rau * Rand1 + sqrt(1 - rau * rau) * Rand2))
        S2next = S2 * growthFactor2

        # Arithmetic mean
        arithMean = 0.5 * (S1next + S2next)
        geoMean = sqrt(S1next * S2next)
        if (option_type == "C"):
            arithPayOff[i] = exp(-r * T) * max((arithMean - K), 0)
            geoPayoff[i] = exp(-r * T) * max((geoMean - K), 0)
        else:
            arithPayOff[i] = exp(-r * T) * max((K - arithMean), 0)
            geoPayoff[i] = exp(-r * T) * max((K - geoMean), 0)
    if(flag=="Yes"):    #control variate version
        geo = Geometric_Basket_Option(S1, S2, sigma_1, sigma_2, r, T, K, rau, option_type, 1)
        #print(geo)
        covXY = np.mean(arithPayOff * geoPayoff) - np.mean(arithPayOff) * np.mean(geoPayoff)
        theta = covXY / np.std(geoPayoff)
        Z = arithPayOff + theta * (geo - geoPayoff)
        Zmean = np.mean(Z)
        Zstd = np.std(Z)
        confcv = [Zmean - 1.96 * Zstd / sqrt(number), Zmean + 1.96 * Zstd / sqrt(number)]
        return round(np.mean(confcv),4)

    else: # Standard monte carlo
        Pmean = np.mean(arithPayOff)
        Pstd = np.std(arithPayOff)
        confmc = [Pmean - 1.96 * Pstd / sqrt(number), Pmean + 1.96 * Pstd / sqrt(number)]
        return round(np.mean(confmc),4)

def a_basket_index(request):
    return render(request, 'a_basket.html')

def a_basket(request,S1, S2, sigma_1, sigma_2, r, T, K, rau, call_put, flag):
    c = Arithmetic_Basket_Option(int(S1), int(S2), float(sigma_1), float(sigma_2), float(r), float(T), float(K), float(rau), str(call_put), str(flag))
    r = HttpResponse(str(c))
    return r
