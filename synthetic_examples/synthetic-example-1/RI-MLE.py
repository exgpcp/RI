#!/usr/bin/env python
from scipy.stats import expon
from scipy.stats import uniform
from scipy.stats import norm
from scipy.stats import multivariate_normal
from numpy.random import multinomial
from numpy.random import uniform
import numpy as np
import matplotlib.pyplot as plt
import scipy 
from scipy.stats import truncnorm
import math
from scipy.stats import mvn
np.set_printoptions(suppress=True)
import pyreadr
import time
import sys

data = pyreadr.read_r('syndata/data_'+sys.argv[1]+'.rda') #args[1] range from 1 to 300
points_inhomo = np.array(data["dataa"]).squeeze()
hyperpara=(pyreadr.read_r('MLE.rda')["groundhypest"]).squeeze()
theta0=hyperpara[0][int(sys.argv[1])-1]
theta1=hyperpara[1][int(sys.argv[1])-1]
T=50
bin_num=1000
x=np.linspace(T/bin_num/2,T-T/bin_num/2,bin_num)
c=math.ceil(int(sys.argv[1])/100)
intensity=c*(2*np.exp(-x/15)+np.exp(-((x-25)/10)**2))

def expo_quad_kernel(theta0,theta1,xn,xm): # 1,0.1
    return theta0*np.exp(-theta1/2*np.sum((xn - xm)**2))

def expo_quad_kernel2(theta0,theta1,xn,T): # 1,0.1
    return np.sqrt(np.pi/2/theta1)*theta0*(math.erf(np.sqrt(theta1/2)*(T-xn))+math.erf(np.sqrt(theta1/2)*xn))

def expo_quad_kernel3(theta0,theta1,T): # 1,0.1
    return 2*theta0/theta1*(np.sqrt(np.pi*theta1/2)*T*math.erf(np.sqrt(theta1/2)*T)+ np.exp(-theta1/2*(T**2)) -1)

def chol_sample(mean, cov_chol):
    return mean + cov_chol @ np.random.standard_normal(mean.size)

def log_lik(f,ns):
    if np.prod(f>0)==0:
        return float('-inf') 
    else:
        return np.sum(np.log(f[:,Ngrid:Ngrid+ns]))-f[:,Ngrid+ns]
    
def elliptical_slice(initial_theta,prior,lnpdf,pdf_params=(),
                     cur_lnpdf=None,angle_range=None):
    """
    NAME:
       elliptical_slice
    PURPOSE:
       Markov chain update for a distribution with a Gaussian "prior" factored out
    INPUT:
       initial_theta - initial vector
       prior - cholesky decomposition of the covariance matrix 
               (like what numpy.linalg.cholesky returns), 
               or a sample from the prior
       lnpdf - function evaluating the log of the pdf to be sampled
       pdf_params= parameters to pass to the pdf
       cur_lnpdf= value of lnpdf at initial_theta (optional)
       angle_range= Default 0: explore whole ellipse with break point at
                    first rejection. Set in (0,2*pi] to explore a bracket of
                    the specified width centred uniformly at random.
    OUTPUT:
       new_theta, new_lnpdf
    HISTORY:
       Originally written in matlab by Iain Murray (http://homepages.inf.ed.ac.uk/imurray2/pub/10ess/elliptical_slice.m)
       2012-02-24 - Written - Bovy (IAS)
    """
    D= len(initial_theta)
    if cur_lnpdf is None:
        cur_lnpdf= lnpdf(initial_theta,*pdf_params)

    # Set up the ellipse and the slice threshold
    if len(prior.shape) == 1: #prior = prior sample
        nu= prior
    else: #prior = cholesky decomp
        if not prior.shape[0] == D or not prior.shape[1] == D:
            raise IOError("Prior must be given by a D-element sample or DxD chol(Sigma)")
        nu= np.dot(prior,np.random.normal(size=D))
    hh = math.log(np.random.uniform()) + cur_lnpdf

    # Set up a bracket of angles and pick a first proposal.
    # "phi = (theta'-theta)" is a change in angle.
    if angle_range is None or angle_range == 0.:
        # Bracket whole ellipse with both edges at first proposed point
        phi= np.random.uniform()*2.*math.pi
        phi_min= phi-2.*math.pi
        phi_max= phi
    else:
        # Randomly center bracket on current point
        phi_min= -angle_range*np.random.uniform()
        phi_max= phi_min + angle_range
        phi= np.random.uniform()*(phi_max-phi_min)+phi_min

    # Slice sampling loop
    while True:
        # Compute xx for proposed angle difference and check if it's on the slice
        xx_prop = initial_theta*math.cos(phi) + nu*math.sin(phi)
        cur_lnpdf = lnpdf(xx_prop,*pdf_params)
        if cur_lnpdf > hh:
            # New point is on slice, ** EXIT LOOP **
            break
        # Shrink slice to rejected point
        if phi > 0:
            phi_max = phi
        elif phi < 0:
            phi_min = phi
        else:
            raise RuntimeError('BUG DETECTED: Shrunk to current position and still not acceptable.')
        # Propose new angle difference
        phi = np.random.uniform()*(phi_max - phi_min) + phi_min
    return (xx_prop,cur_lnpdf)

K=len(points_inhomo)
N=K+1
Ngrid=100
xxx=np.linspace(0,T,Ngrid+1)[:Ngrid] 
g_mk=3*c+np.zeros((Ngrid+K+1))## 0--(Ngrid-1):function values at grids, Ngrid--(Ngrid+K-1):K observations, last:integral term.
g_mk[Ngrid+K]=150*c
Nfinal=Ngrid+N
x4=np.concatenate([xxx, points_inhomo])
cov_K=np.zeros((Nfinal,Nfinal))
noise_var=1e-13
nsim1=10000
nsim2=50000
for i in range(Nfinal):
        for j in range(i,Nfinal):
            if (j<(Nfinal-1)) and (i<(Nfinal-1)):
                cov_K[i][j]=expo_quad_kernel(theta0,theta1,x4[i],x4[j])
            if (j==(Nfinal-1)) and (i<(Nfinal-1)):
                cov_K[i][j]=expo_quad_kernel2(theta0,theta1,x4[i],T)
            if (j==(Nfinal-1)) and (i==(Nfinal-1)):
                cov_K[i][j]=expo_quad_kernel3(theta0,theta1,T)  
            if j!=i:
                cov_K[j][i]=cov_K[i][j]
cov_K_noise=cov_K
min_eig=np.min(np.real(np.linalg.eigvals(cov_K_noise)))
while(min_eig<1e-10):
    cov_K_noise += np.eye(Nfinal)*noise_var
    min_eig=np.min(np.real(np.linalg.eigvals(cov_K_noise)))
cov_K_chol=np.linalg.cholesky(cov_K_noise)

for ite in range(nsim1):
    prior=chol_sample(mean=np.zeros(Nfinal), cov_chol=cov_K_chol)
    g_mk,curloglike=elliptical_slice(g_mk.reshape(1,-1),prior,log_lik,pdf_params=[K],cur_lnpdf=None,angle_range=None)
 
g_mk_list=[]
g_mk_list2=[]
g_mk_list3=[]
start_time2 = time.time()
for ite in range(nsim2):
    prior=chol_sample(mean=np.zeros(Nfinal), cov_chol=cov_K_chol)
    g_mk,curloglike=elliptical_slice(g_mk.reshape(1,-1),prior,log_lik,pdf_params=[K],cur_lnpdf=None,angle_range=None)
    g_mk_list.append(g_mk[0][Ngrid:Ngrid+K])
    g_mk_list2.append(g_mk[0][0:Ngrid])
    g_mk_list3.append(g_mk[0][Ngrid+K])
timerun2=time.time()-start_time2  

low=np.quantile(g_mk_list2, 0.025, axis=0)
me=np.quantile(g_mk_list2, 0.5, axis=0)
high=np.quantile(g_mk_list2, 0.975, axis=0)

#about the integral values
comp1=30*(1-np.exp(-10/3))
comp2=10*np.sqrt(np.pi)*math.erf(5/2)
integral_truth=(comp1+comp2)*c
integral_mean=np.mean(g_mk_list3)
integral_sd=np.sqrt(np.cov(g_mk_list3))
np.savez('/output/simulation/mymethodfinal3/ESSsyn1_truthmle'+sys.argv[1]+'.npz', aaa=g_mk_list3,aa=g_mk_list2,a=g_mk_list,c=points_inhomo,d=xxx,b=x,e=intensity,g=theta0,h=theta1,i=noise_var,m=integral_truth,n=integral_mean,o=integral_sd,q=timerun2)
