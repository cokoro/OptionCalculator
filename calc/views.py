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


def browianMotion(S0, dt, c1, c2, random):
    return S0 * np.exp(c1 * dt + c2 * random)


def arith_asian_option(S0,sigma,r,T,K,step,option_type,path,flag):
    
    #This method is called when RandomState is initialized. It can be called again to re-seed the generator.
    np.random.seed(0)
    # get all paths
    paths = np.zeros((path, step))
    random = np.zeros((path, step))
    # generator random to each paths
    for i in range(0, path):
        random[i, :] = np.random.standard_normal(step)
    c1 = r - 0.5 * sigma * sigma
    dt = float(T) / step
    c2 = sigma * np.sqrt(dt)
    # set first column
    paths[:, 0] = browianMotion(S0, dt, c1, c2, random[:, 0])

    # set 1 to step-1 column
    for i in range(1, step):
        s = paths[:, i - 1]
        paths[:, i] = browianMotion(S0, dt, c1, c2, random[:, i])

    #mean(1) calculate mean of row
    #mean(0) calculate mean of column
    arithMean = paths.mean(1)
    getLogList = np.log(paths)
    geoMean = np.exp(1 / float(step) * getLogList.sum(1))

    if option_type == 'C':
        arith = np.maximum(arithMean - K, 0)*exp(-r*(T))
        geo = np.maximum(geoMean - K, 0)*np.exp(-r*(T))
    elif option_type == 'P':
        arith = np.maximum(K - arithMean, 0)*exp(-r*(T))
        geo = np.maximum(K - geoMean, 0)*exp(-r*(T))
   
    #Standard Monte Carlo
    if flag == 'No':
        pmean=np.mean(arith)
        pstd=np.std(arith)
        confcv=[pmean-1.96*pstd/np.sqrt(path),pmean+1.96*pstd/np.sqrt(path)]
        result=(pmean-1.96*pstd/np.sqrt(path)+pmean+1.96*pstd/np.sqrt(path))*0.5
        return round(result,4)

    #Control variates
    elif flag=='Yes':
        XY = arith * geo
        covXY = np.mean(XY) - (np.mean(geo) * np.mean(arith))
        theta = covXY/np.var(geo)
        geoMean = geo_asian_option(S0,sigma,r,T,K,step,option_type)
        z = arith + theta*(geoMean - geo)
        Zmean=np.mean(z)
        Zstd=np.std(z)
        confcv=[Zmean-1.96*Zstd/np.sqrt(path),Zmean+1.96*Zstd/np.sqrt(path)]
        result=(Zmean-1.96*Zstd/np.sqrt(path)+Zmean+1.96*Zstd/np.sqrt(path))*0.5
        return round(result,4)

def index(request):
    return render(request, 'index.html')


def add(request, S, K, T, sigma, r,call_put,path,step,controlspecify):
    print "call = ",call_put
    #import pdb; pdb.set_trace()

    c = arith_asian_option(int(S), float(sigma), float(r), float(T), float(K), int(step), str(call_put), int(path), str(controlspecify))
    r = HttpResponse(str(c))
    return r
