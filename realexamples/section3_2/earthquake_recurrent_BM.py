#!/usr/bin/env python
# coding: utf-8

# In[5]:


#!/usr/bin/env python
# coding: utf-8

# In[2]:


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


data = np.load('/earthquakes/earthquakes.npz')
points_inhomo=list(data['a'])


T=365

def expo_quad_kernel(xn,xm): # 1,0.1
    return min(xn,xm)

def expo_quad_kernel2(xn,T): # 1,0.1
    return (T*xn-0.5*xn**2)

def expo_quad_kernel3(T): # 1,0.1
    return (T**3/3)



# ## MCMC inference

# In[7]:

def chol_sample(mean, cov_chol):
    #return mean + np.linalg.cholesky(cov) @ np.random.standard_normal(mean.size)
    return mean + cov_chol @ np.random.standard_normal(mean.size)


def log_lik(f,ns):
    #x=pdf_params[0]
    #ns=pdf_params[1]
    #pr=scipy.stats.norm.logpdf(x[:,1],mu_y+rho0*math.sqrt(sigmasq_y)/math.sqrt(sigmasq_x)*(x[:,0]-mu_x) ,math.sqrt(sigmasq_y*(1-rho0**2)))

    if np.prod(f>0)==0:
        return float('-inf') 
    else:
        #print(f[:,0:ns].shape)
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
xxx=np.linspace(0,T,Ngrid+1)[1:Ngrid+1] 
g_mk=3+np.zeros((Ngrid+K+1))## 0--(Ngrid-1) Ngrid--(Ngrid+K-1) observe K--last unobserve +integral term
g_mk[Ngrid+K]=3*T
Nfinal=Ngrid+N
x4=np.concatenate([xxx, points_inhomo])

cov_K=np.zeros((Nfinal,Nfinal))
noise_var=1e-7
theta=np.exp(4)
rem=1
alpha=.1
beta=.1
for i in range(Nfinal):
    for j in range(i,Nfinal):
        if (j<(Nfinal-1)) and (i<(Nfinal-1)):
            cov_K[i][j]=expo_quad_kernel(x4[i],x4[j])
        if (j==(Nfinal-1)) and (i<(Nfinal-1)):
            cov_K[i][j]=expo_quad_kernel2(x4[i],T)
        if (j==(Nfinal-1)) and (i==(Nfinal-1)):
            cov_K[i][j]=expo_quad_kernel3(T)  
        if j!=i:
            cov_K[j][i]=cov_K[i][j]



cov_K_inv=np.linalg.inv(cov_K+np.eye(Nfinal)*noise_var)

ones_vec=np.ones(Nfinal)
ones_vec[Nfinal-1]=T
ones_vec=ones_vec.reshape(Nfinal,1)
cov_K_mod_inv=cov_K_inv-cov_K_inv@ones_vec@ones_vec.transpose()@cov_K_inv/(ones_vec.transpose()@cov_K_inv@ones_vec)


cov_K_mod_inv_add=cov_K_mod_inv+np.eye(Nfinal)*noise_var

cov_K_mod=np.linalg.inv(cov_K_mod_inv_add)

cov_K_mod_chol=np.linalg.cholesky(cov_K_mod) 

nsim1=10000
nsim2=50000
#L_inv=np.linalg.cholesky(cov_K_inv_add)
#cov_K_chol=np.linalg.inv(L_inv)
for ite in range(nsim1):
    if (ite%10000==0):
        print(ite)#cov_K is the covariance matrix
    cov_K_chol_final=cov_K_mod_chol/np.sqrt(theta)
    prior=chol_sample(mean=np.zeros(Nfinal), cov_chol=cov_K_chol_final)#nu
    g_mk,curloglike=elliptical_slice(g_mk.reshape(1,-1),prior,log_lik,pdf_params=[K],cur_lnpdf=None,angle_range=None)
   
    alpha_pos=alpha+Nfinal/2
    beta_pos=beta+1/2*g_mk[0]@cov_K_mod_inv@g_mk[0]
    theta=np.random.gamma(alpha_pos, 1/beta_pos,1)
   
g_mk_list=[]
g_mk_list2=[]
g_mk_list3=[]
theta_list=[]
#theta=0.01

for ite in range(nsim2):
    if (ite%10000==0):
        print(ite)#cov_K is the covariance matrix
    cov_K_chol_final=cov_K_mod_chol/np.sqrt(theta)
    prior=chol_sample(mean=np.zeros(Nfinal), cov_chol=cov_K_chol_final)#nu
    g_mk,curloglike=elliptical_slice(g_mk.reshape(1,-1),prior,log_lik,pdf_params=[K],cur_lnpdf=None,angle_range=None)
    g_mk_list.append(g_mk[0][Ngrid:Ngrid+K])
    g_mk_list2.append(g_mk[0][0:Ngrid])
    g_mk_list3.append(g_mk[0][Ngrid+K])
    alpha_pos=alpha+Nfinal/2
    beta_pos=beta+1/2*g_mk[0]@cov_K_mod_inv@g_mk[0]
    theta=np.random.gamma(alpha_pos, 1/beta_pos,1)
    theta_list.append(theta)    #print(beta_pos)
        #print(theta)



np.savez('/output/graphs/real2/BM.npz', aaa=g_mk_list3,aa=g_mk_list2,a=g_mk_list,c=points_inhomo,d=xxx,i=noise_var,p=timerun,q=width1,r=width2,s=theta_list)

#1e-7
plt.figure(1)                # the first figure
plt.subplot(111)             # the first subplot in the first figure
low=np.quantile(g_mk_list2, 0.025, axis=0)
plt.plot(xxx,low,'--',lw=1,alpha=0.6,c='tab:grey',label='2.5%-50%-97.5% Quantiles')
me=np.quantile(g_mk_list2, 0.5, axis=0)
plt.plot(xxx,me,'--',lw=1,alpha=0.6,c='tab:grey')
high=np.quantile(g_mk_list2, 0.975, axis=0)
plt.plot(xxx,high,'--',lw=1,alpha=0.6,c='tab:grey')
plt.plot(points_inhomo,np.zeros(len(points_inhomo)),linestyle='None', marker='|', color='k', markersize=10,label='events')
plt.ylim(0,6)
plt.ylabel(r'Poisson Intensity $\lambda(s)$')
plt.title('Earthquake dataset with RI-BM fit',y=-0.2)
plt.legend()
plt.savefig("/output/graphs/real2/BM_median.pdf",bbox_inches='tight')
plt.show()


plt.figure()
plt.hist(g_mk_list3,bins=20)
plt.xlim(700,1200)
plt.title('Histogram of $\int_0^T\lambda(s)ds$', y=-0.2)
plt.savefig("/output/graphs/real2/BM_histogram.pdf",bbox_inches='tight')
plt.show()

plt.figure()
plt.plot(np.array(g_mk_list2)[:,50])

plt.title('latent GP at midpoint',y=-0.2)
plt.savefig("/output/graphs/real2/BM_traceplot1.pdf",bbox_inches='tight')
plt.show()



