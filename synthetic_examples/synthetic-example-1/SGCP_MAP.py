#!/usr/bin/env python
# coding: utf-8

from scipy.stats import expon
from scipy.stats import uniform
from scipy.stats import norm
from scipy.stats import multivariate_normal
from numpy.random import multinomial
from scipy.stats import uniform
import numpy as np
import matplotlib.pyplot as plt
import scipy 
from scipy.stats import truncnorm
import math
import pyreadr
import time
import sys

def expo_quad_kernel(theta0,theta1,xn,xm): # 1,0.1
    return theta0*np.exp(-theta1/2*np.sum((xn - xm)**2))

def sigmoid(x): #"Numerically-stable sigmoid function."
    result=np.array([])
    for x_i in x:
        if x_i >= 0:
            z = np.exp(-x_i)
            result=np.append(result, 1 / (1 + z))
        else:
            z = np.exp(x_i)
            result=np.append(result, z / (1 + z))
    return result

def GP_regression(xi,yi,theta0,theta1,noise_var,rang,num_points):
    N=len(xi)
    cov_K=np.zeros((N,N))
    for i in range(N):
        for j in range(i,N):
            cov_K[i][j]=expo_quad_kernel(theta0,theta1,xi[i],xi[j])
            cov_K[j][i]=cov_K[i][j]
    min_eig=np.min(np.real(np.linalg.eigvals(cov_K))) # numerical float truncation error refine
    while(min_eig<0):
        cov_K += -10*min_eig*np.eye(N)
        min_eig=np.min(np.real(np.linalg.eigvals(cov_K)))
    cov_K_noise=cov_K+np.eye(N)*noise_var
    x1=np.linspace(0,rang,num_points+1)      # prediction points, integer is to make it easy
    M=len(x1)-1
    mean=np.zeros((1,M))[0]
    posterior_cov=np.zeros((M,M))
    k_matrix=np.zeros((M,N))
    k_matrix_pre=np.zeros((M,M))
    for i in range(M):
        for j in range(N):
            k_matrix[i][j]=expo_quad_kernel(theta0,theta1,x1[i],xi[j])
    k_C=np.dot(k_matrix,np.linalg.inv(cov_K_noise))
    mean=np.dot(k_C,yi)
    for i in range(M):
        for j in range(i,M):
            k_matrix_pre[i][j]=expo_quad_kernel(theta0,theta1,x1[i],x1[j])
            k_matrix_pre[j][i]=k_matrix_pre[i][j]
    posterior_cov=k_matrix_pre-np.dot(k_C,k_matrix.T)+np.eye(M)*noise_var
    min_eig=np.min(np.real(np.linalg.eigvals(posterior_cov))) # numerical float truncation error refine
    while(min_eig<0):
        posterior_cov += -10*min_eig*np.eye(posterior_cov.shape[0])
        min_eig=np.min(np.real(np.linalg.eigvals(posterior_cov)))
    return x1[:num_points],mean, posterior_cov

def GP_regression_one_pred(xi,yi,theta0,theta1,noise_var,x_pred):
    N=len(xi)
    cov_K=np.zeros((N,N))
    for i in range(N):
        for j in range(i,N):
            cov_K[i][j]=expo_quad_kernel(theta0,theta1,xi[i],xi[j])
            if j!=i:
                cov_K[j][i]=cov_K[i][j]

    cov_K_noise=cov_K+np.eye(N)*noise_var
    min_eig=np.min(np.real(np.linalg.eigvals(cov_K_noise))) # numerical float truncation error refine
    i=0
    while(min_eig<0):
        cov_K_noise += -10*min_eig*np.eye(N)
        i+=1
        min_eig=np.min(np.real(np.linalg.eigvals(cov_K_noise)))
        if i>10: break
    
    # prediction point
    k_vec=np.array([])
    for i in range(N):
        k_vec=np.append(k_vec,expo_quad_kernel(theta0,theta1,x_pred,xi[i]))
    k_C=np.dot(k_vec,np.linalg.inv(cov_K_noise))
    mean=np.dot(k_C,yi)
    var=theta0+noise_var-np.dot(k_C,k_vec)
    if var<0:var=0.0001
    std_dev=np.sqrt(var)
    
    return mean,std_dev






data = pyreadr.read_r('/syndata/data_'+sys.argv[1]+'.rda')
points_inhomo = list(np.array(data["dataa"]).squeeze())

hyperpara=(pyreadr.read_r('/MAP/func1_priorest'+sys.argv[1]+'.rda')["b"]).squeeze()
theta0=hyperpara[0]
theta1=hyperpara[1]
c=math.ceil(int(sys.argv[1])/100)
measure_sup=3*c
noise_var=1e-4
T=50

def sampling_M(M,s_m,g_mk,measure_sup,T,theta0,theta1,noise_var):
    temp=uniform.rvs(0,1)
    if temp<=0.5: ## M-->M+1
        s_prime=scipy.stats.uniform.rvs(0,T) # propose a s'
        mean_s_prime,std_s_prime=GP_regression_one_pred(np.array(points_inhomo+list(s_m)),g_mk,theta0,theta1,noise_var,s_prime)
        g_s_prime=norm.rvs(mean_s_prime,std_s_prime) # propose a g(s')
        accept_rate=T*measure_sup/(M+1)/(1+np.exp(g_s_prime))
        temp_1=scipy.stats.uniform.rvs(0,1)
        if temp_1<accept_rate:
            M=M+1
            s_m=np.append(s_m,s_prime)
            g_mk=np.append(g_mk,g_s_prime)
            assert len(s_m)==M
            assert len(g_mk)==M+K
    else:   ## M-->M-1
        index=np.where(multinomial(1,[1/M]*M)==1)[0][0]
        s_prime=s_m[index]
        g_s_prime=g_mk[K+index]
        accept_rate=M*(1+np.exp(g_s_prime))/T/measure_sup
        temp_1=scipy.stats.uniform.rvs(0,1)
        if temp_1<accept_rate:
            M=M-1
            s_m=np.delete(s_m,index)
            g_mk=np.delete(g_mk,K+index)
            assert len(s_m)==M
            assert len(g_mk)==M+K
    return M, s_m, g_mk

def sampling_s_m(M,s_m,g_mk,measure_sup,T,theta0,theta1,noise_var):
    for i in range(M):
        s_prime=norm.rvs(s_m[i],1)
        while (s_prime>T or s_prime<0):
            s_prime=norm.rvs(s_m[i],1)
        s_m_1=np.delete(s_m,i)
        g_mk_1=np.delete(g_mk,K+i)
        mean_s_prime,std_s_prime=GP_regression_one_pred(np.array(points_inhomo+list(s_m_1)),g_mk_1,theta0,theta1,noise_var,s_prime)
        g_s_prime=norm.rvs(mean_s_prime,std_s_prime) # propose a g(s')
        accept_rate=(1+np.exp(g_mk[K+i]))/(1+np.exp(g_s_prime))
        temp_1=uniform.rvs(0,1)
        if temp_1<accept_rate:
            s_m=np.insert(s_m_1,i,s_prime)
            g_mk=np.insert(g_mk_1,K+i,g_s_prime)
            assert len(s_m)==M
            assert len(g_mk)==M+K
    return M, s_m, g_mk

M=20    ## initial
K=len(points_inhomo)
s_m=np.linspace(1,T-1,M)
g_mk=np.zeros((M+K))## 0--(K-1) observe K--last unobserve

def chol_sample(mean, cov,jit):
    return mean + np.linalg.cholesky(cov+np.diag(np.full(cov.shape[0], jit))) @ np.random.standard_normal(mean.size)

def log_lik_adam(f,ns,cov_inv):
    #ns is len(points_inhomo)
    #f[0:ns] is observations
    #f[ns:] is thinned
    #print(-0.5*f@cov_inv@f.transpose()-np.sum(np.log(1+np.exp(-f[0][0:ns])))-np.sum(np.log(1+np.exp(f[0][ns:]))))
    #print(-0.5*f@cov_inv@f.transpose()-np.sum(np.log(1+np.exp(-f[:,0:ns])))-np.sum(np.log(1+np.exp(f[:,ns:]))))
    #comp1=-0.5*f@cov_inv@f.transpose()
    return -np.sum(np.log(1+np.exp(-f[:,0:ns])))-np.sum(np.log(1+np.exp(f[:,ns:])))
    
def elliptical_slice_adam(initial_theta,prior,lnpdf,pdf_params=(),
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

# In[105]:

nsim1=10000
nsim2=10000
start_time1=time.time()
for k in range(nsim1):
    M,s_m,g_mk=sampling_M(M,s_m,g_mk.squeeze(),measure_sup,T,theta0,theta1,noise_var)
    M,s_m,g_mk=sampling_s_m(M,s_m,g_mk.squeeze(),measure_sup,T,theta0,theta1,noise_var) # sampling s_m
    loc=np.array(points_inhomo+list(s_m))
    cov_K=np.zeros((K+M,K+M))
    for i in range(K+M):
        for j in range(i,K+M):
            cov_K[i][j]=expo_quad_kernel(theta0,theta1,loc[i],loc[j])
            cov_K[j][i]=cov_K[i][j]
    cov_K_noise=cov_K+np.eye(K+M)*noise_var
    cov_K_noise_inv=np.linalg.inv(cov_K_noise)
    prior=chol_sample(mean=np.zeros(K+M), cov=cov_K, jit=noise_var)#nu
    g_mk,curloglike=elliptical_slice_adam(g_mk.reshape(1,-1),prior,log_lik_adam,pdf_params=[K,cov_K_noise_inv],cur_lnpdf=None,angle_range=None)
timerun1=time.time()-start_time1

M_list=[] 
s_m_list=[]
g_mk_list3=[]
g_mk_list2=[]
g_mk_list=[]
Ngrid=100
start_time2=time.time()
for k in range(nsim2):
    M,s_m,g_mk=sampling_M(M,s_m,g_mk.squeeze(),measure_sup,T,theta0,theta1,noise_var)
    M,s_m,g_mk=sampling_s_m(M,s_m,g_mk.squeeze(),measure_sup,T,theta0,theta1,noise_var) # sampling s_m
    loc=np.array(points_inhomo+list(s_m))
    cov_K=np.zeros((K+M,K+M))
    for i in range(K+M):
        for j in range(i,K+M):
            cov_K[i][j]=expo_quad_kernel(theta0,theta1,loc[i],loc[j])
            cov_K[j][i]=cov_K[i][j]
    cov_K_noise=cov_K+np.eye(K+M)*noise_var
    cov_K_noise_inv=np.linalg.inv(cov_K_noise)
    prior=chol_sample(mean=np.zeros(K+M), cov=cov_K, jit=noise_var)#nu
    g_mk,curloglike=elliptical_slice_adam(g_mk.reshape(1,-1),prior,log_lik_adam,pdf_params=[K,cov_K_noise_inv],cur_lnpdf=None,angle_range=None)    
    xxx,mean,cov=GP_regression(np.array(points_inhomo+list(s_m)),g_mk.squeeze(),theta0,theta1,noise_var,T,Ngrid)
    line=multivariate_normal.rvs(mean,cov)
    M_list.append(M)
    s_m_list.append(s_m)
    g_mk_list3.append(g_mk)
    g_mk_list2.append(line)
    g_mk_list.append(g_mk[0][0:K])
timerun2=time.time()-start_time2    

inten_est=measure_sup/(1+np.exp(-np.array(g_mk_list)))
low=np.quantile(inten_est, 0.025, axis=0)
med=np.quantile(inten_est, 0.5, axis=0)
high=np.quantile(inten_est, 0.975, axis=0)
truth=2*np.exp(-np.array(points_inhomo)/15)+np.exp(-((np.array(points_inhomo)-25)/10)**2)
truth=c*truth
l2_dist1=sum((np.array(med).squeeze()-truth)**2)
coverage1=np.sum((truth>=low.squeeze()) * (truth<=high.squeeze()))/K
width1=sum(high-low)/K

inten_est2=measure_sup/(1+np.exp(-np.array(g_mk_list2)))
low=np.quantile(inten_est2, 0.025, axis=0)
med=np.quantile(inten_est2, 0.5, axis=0)
high=np.quantile(inten_est2, 0.975, axis=0)

truth=2*np.exp(-np.array(xxx)/15)+np.exp(-((np.array(xxx)-25)/10)**2)
truth=c*truth
l2_dist2=sum((np.array(med).squeeze()-truth)**2)
coverage2=np.sum((truth>=low.squeeze()) * (truth<=high.squeeze()))/Ngrid
width2=sum(high-low)/Ngrid


np.savez('/output/simulation/adam1_1/priormle/syn1/syn1_priormle_'+sys.argv[1]+'.npz', aaa=M_list,aa=g_mk_list2,a=g_mk_list,c=points_inhomo,d=xxx,
    e=theta0,f=theta1,g=measure_sup,h=noise_var,i=coverage1,j=coverage2,k=l2_dist1,l=l2_dist2,m=width1,n=width2,o=timerun1,p=timerun2)


import pickle
with open("/output/simulation/adam1_1/priormle/syn1/syn1_priormle1_"+sys.argv[1]+".bin", "wb") as output:
    pickle.dump(g_mk_list, output)

with open("/output/simulation/adam1_1/priormle/syn1/syn1_priormle2_"+sys.argv[1]+".bin", "wb") as output:
    pickle.dump(s_m_list, output)
