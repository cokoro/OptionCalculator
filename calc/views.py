from django.shortcuts import render
from django.http import HttpResponse

from math import *
import numpy as np
from scipy.stats import norm
import scipy
from math import sqrt,exp,log
import math
from scipy import stats

def arith_asian_option(S,sigma,r,T,E,N,Type,M,flag):
    Dt = T/N #Dt = 3/50=0.06
    n = N

    sigsqT = np.square(sigma)*T*(N+1)*(2*N+1)/(6*N*N)
    muT = 0.5*sigsqT+(r-0.5*np.square(sigma))*T*(N+1)/(2*N)
    
    d1 = (math.log(S/E)+(muT+0.5*sigsqT))/(math.sqrt(sigsqT))
    d2 = d1 - math.sqrt(sigsqT)
    
    drift = math.exp((r-0.5*np.square(sigma))*Dt)    
    

    if Type == "C":
        Spath = [0]*N
        arithPayoff = [0]*M
        geoPayoff = [0]*M
        #M unmber of samples
        geo = math.exp(-r*T)*(S*math.exp(muT)*stats.norm.cdf(d1)-E*stats.norm.cdf(d2)) #closed-form value
        for i in range(M):
            Spath = [0]* N
            growthFactor = drift*math.exp(sigma*math.sqrt(Dt)*np.random.normal())
            Spath[0] = S * growthFactor
            for j in range(1,N):
                growthFactor = drift*math.exp(sigma*math.sqrt(Dt)*np.random.normal())
                Spath[j] = Spath[j-1] * growthFactor
        
            #arithmetic mean
            arithMean = np.mean(Spath)
            arithPayoff[i] = math.exp(-r*T)*max(arithMean-E,0)  #payoffs

            
            #geometric mean
            SpathLog = [0]*N
            for n in range(len(Spath)):
                SpathLog[n]=math.log(Spath[n])
                
            geoMean = math.exp((1/N)*np.sum(SpathLog))
            geoPayoff[i] = math.exp(-r*T)*max(geoMean-E,0)
                
   
        if(flag =="Yes"):
			Pmean = np.mean(arithPayoff)
			Pstd = np.std(arithPayoff)
			confmc = [Pmean-1.96*Pstd/math.sqrt(M), Pmean+1.96*Pstd/math.sqrt(M)] 
			return round(Pmean,4) 
        else:
        	covXY = np.mean(np.multiply(arithPayoff, geoPayoff))-np.mean(arithPayoff)*np.mean(geoPayoff)
        	theta = covXY/np.var(geoPayoff)
        	Z = np.array(arithPayoff) + theta * (geo - np.array(geoPayoff))
        	Zmean = np.mean(Z)
        	Zstd = np.std(Z)
        	confcv = [Zmean-1.96*Zstd/math.sqrt(M), Zmean+1.96*Zstd/math.sqrt(M)]
        	return round(Zmean,4)
    else:
        Spath = [0]*N
        arithPayoff = [0]*M
        geoPayoff = [0]*M
        geo = math.exp(-r*T)*(E*stats.norm.cdf(-d2)-S*math.exp(muT)*stats.norm.cdf(-d1)) #closed-form value
        for i in range(M):
            Spath = [0]* N
            growthFactor = drift*math.exp(sigma*math.sqrt(Dt)*np.random.normal())
            Spath[0] = S * growthFactor
            for j in range(1,N):
                growthFactor = drift*math.exp(sigma*math.sqrt(Dt)*np.random.normal())
                Spath[j] = Spath[j-1] * growthFactor
        
            #arithmetic mean
            arithMean = np.mean(Spath)
            arithPayoff[i] = math.exp(-r*T)*max(E-arithMean,0)  #payoffs

            
            #geometric mean
            SpathLog = [0]*N
            for n in range(len(Spath)):
                SpathLog[n]=math.log(Spath[n])
                
            geoMean = math.exp((1/N)*np.sum(SpathLog))
            geoPayoff[i] = math.exp(-r*T)*max(E-geoMean,0)
                
        if(flag =="Yes"):
            Pmean = np.mean(arithPayoff)
            Pstd = np.std(arithPayoff)
            confmc = [Pmean-1.96*Pstd/math.sqrt(M), Pmean+1.96*Pstd/math.sqrt(M)]   
            return round(Pmean,4)
        else:
			covXY = np.mean(np.multiply(arithPayoff, geoPayoff))-np.mean(arithPayoff)*np.mean(geoPayoff)
			theta = covXY/np.var(geoPayoff)
			Z = np.array(arithPayoff) + theta * (geo - np.array(geoPayoff))
			Zmean = np.mean(Z)
			Zstd = np.std(Z)
			confcv = [Zmean-1.96*Zstd/math.sqrt(M), Zmean+1.96*Zstd/math.sqrt(M)]
			return round(Zmean,4)

def index(request):
    return render(request, 'index.html')


def add(request, S, K, T, sigma, r,call_put,path,step,controlspecify):
    print "call = ",call_put
    #import pdb; pdb.set_trace()

    c = arith_asian_option(int(S), float(sigma), float(r), float(T), float(K), int(step), str(call_put), int(path), str(controlspecify))
    print "res = ", c
    r = HttpResponse(str(c))
    return r
