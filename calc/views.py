from django.shortcuts import render
from django.http import HttpResponse

from math import *
import numpy as np
from scipy.stats import norm
import scipy
from math import sqrt,exp,log

def geo_asian_option(S0,sigma,r,T,K,option_type):
    dt = 1e-2
    N = int(T / dt)
    
    sigmaHat = sigma * sqrt((N+1)*(2*N+1)/(6*pow(N,2)))
    muHat = (r-0.5*pow(sigma, 2))*(N+1)/(2*N)+0.5*pow(sigmaHat, 2)
    
    d1 = ( log( S0 / K )+ (muHat + 0.5 * pow(sigmaHat,2)) * (T)) / (sigmaHat * sqrt(T))
    d2 = d1 - sigmaHat*sqrt(T)

    if option_type == 'C':
        return round(exp(-r*(T))*(S0 * exp(muHat*(T))*norm.cdf(d1) - K*norm.cdf(d2)),4)
    elif option_type == 'P':
        return round(exp(-r*(T))*(-S0 * exp(muHat*(T))*norm.cdf(-d1) + K*norm.cdf(-d2)),4)

def browianMotion(S0, dt, c1, c2, random):
    return S0 * np.exp(c1 * dt + c2 * random)


def arith_asian_option(S0,sigma,r,T,K,step,option_type,path,flag):
    dt = 1e-2
    N = int(T / dt)
    drift = exp((r - 0.5 * sigma * sigma) * dt)
    arithPayOff = np.empty(path, dtype=float)
    geoPayoff = np.empty(path, dtype=float)
    scipy.random.seed(5)
    for i in range(0,path):
        Rand = scipy.random.randn(N-1)
        Spath = np.zeros(N)
        Spath[0] = S0
        for j in range(1,N):
            growthFactor = drift * exp(sigma * sqrt(dt) * Rand[j-1])
            Spath[j] = Spath[j - 1] * growthFactor
        arithMean = np.mean(Spath)
        for k in range(0,N):
            Spath[k] = log(Spath[k])
        geoMean = exp((1 / N) * sum(Spath))

        if (option_type == "C"):
            arithPayOff[i] = exp(-r * T) * max(arithMean - K, 0)
            geoPayoff[i] = exp(-r * T) * max(geoMean - K, 0)
        else:
            arithPayOff[i] = exp(-r * T) * max(K - arithMean, 0)
            geoPayoff[i] = exp(-r * T) * max(K - geoMean, 0)

    #Standard Monte Carlo
    if flag == 'No':
        pmean=np.mean(arithPayOff)
        pstd=np.std(arithPayOff)
        confcv=[pmean-1.96*pstd/np.sqrt(path),pmean+1.96*pstd/np.sqrt(path)]
        return round(np.mean(confcv),4)

    #Control variates
    elif flag=='Yes':
        covXY = np.mean(arithPayOff * geoPayoff) - np.mean(arithPayOff) * np.mean(geoPayoff)
        theta = covXY / np.std(geoPayoff)
        geo = geo_asian_option(S0, sigma, r, T, K, option_type)
        Z = arithPayOff + theta * (geo - geoPayoff)
        Zmean = np.mean(Z)
        Zstd = np.std(Z)
        confcv = [Zmean - 1.96 * Zstd / sqrt(path), Zmean + 1.96 * Zstd / sqrt(path)]
        return round(np.mean(confcv), 4)

def index(request):
    return render(request, 'index.html')


def add(request, S, K, T, sigma, r,call_put,path,step,controlspecify):
    print "call = ",call_put
    #import pdb; pdb.set_trace()

    c = arith_asian_option(int(S), float(sigma), float(r), float(T), float(K), int(step), str(call_put), int(path), str(controlspecify))
    r = HttpResponse(str(c))
    return r
