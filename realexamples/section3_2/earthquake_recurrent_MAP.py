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

theta0= 4.135945
theta1=4.756541e-05

def expo_quad_kernel(theta0,theta1,xn,xm): # 1,0.1
    return theta0*np.exp(-theta1/2*np.sum((xn - xm)**2))

def expo_quad_kernel2(theta0,theta1,xn,T): # 1,0.1
    return np.sqrt(np.pi/2/theta1)*theta0*(math.erf(np.sqrt(theta1/2)*(T-xn))+math.erf(np.sqrt(theta1/2)*xn))

def expo_quad_kernel3(theta0,theta1,T): # 1,0.1
    return 2*theta0/theta1*(np.sqrt(np.pi*theta1/2)*T*math.erf(np.sqrt(theta1/2)*T)+ np.exp(-theta1/2*(T**2)) -1)

T=365

def chol_sample(mean, cov_chol):
    #return mean + np.linalg.cholesky(cov) @ np.random.standard_normal(mean.size)
    return mean + cov_chol @ np.random.standard_normal(mean.size)

def log_lik(f,ns):
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

def GP_regression_ess(xi,yi,theta0,theta1,rang,num_points,cov):
    N=len(yi)
   
     #cov is sigma_22
    
    x1=np.linspace(0,rang,num_points+1)      # prediction points, integer is to make it easy
    M=len(x1)-1
    mean=np.zeros((1,M))[0]
    posterior_cov=np.zeros((M,M)) #final covariance returned
    k_matrix=np.zeros((M,N)) #Sigma_12
    
    k_matrix_pre=np.zeros((M,M)) #Sigma_11
    for i in range(M):
        for j in range(N):
            if j<(N-1):
                k_matrix[i][j]=expo_quad_kernel(theta0,theta1,x1[i],xi[j])
            else:
                k_matrix[i][j]=expo_quad_kernel2(theta0,theta1,x1[i],T)
              
    k_C=np.dot(k_matrix,np.linalg.inv(cov))
   
    mean=np.dot(k_C,yi)
    for i in range(M):
        for j in range(i,M):
            k_matrix_pre[i][j]=expo_quad_kernel(theta0,theta1,x1[i],x1[j])
            k_matrix_pre[j][i]=k_matrix_pre[i][j]
    posterior_cov=k_matrix_pre-np.dot(k_C,k_matrix.T)
    
    return x1[:num_points],mean, posterior_cov


K=len(points_inhomo)
N=K+1
Ngrid=100
xxx=np.linspace(0,T,Ngrid+1)[:Ngrid] 
g_mk=3+np.zeros((Ngrid+K+1))## 0--(Ngrid-1) Ngrid--(Ngrid+K-1) observe K--last unobserve +integral term
g_mk[Ngrid+K]=3*T
 
Nfinal=Ngrid+N
x4=np.concatenate([xxx, points_inhomo])
cov_K=np.zeros((Nfinal,Nfinal))
noise_var=1e-11
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
    prior=chol_sample(mean=np.zeros(Nfinal), cov_chol=cov_K_chol)#nu
    g_mk,curloglike=elliptical_slice(g_mk.reshape(1,-1),prior,log_lik,pdf_params=[K],cur_lnpdf=None,angle_range=None)

g_mk_list=[]
g_mk_list2=[]
g_mk_list3=[]
for ite in range(nsim2):
    prior=chol_sample(mean=np.zeros(Nfinal), cov_chol=cov_K_chol)
    g_mk,curloglike=elliptical_slice(g_mk.reshape(1,-1),prior,log_lik,pdf_params=[K],cur_lnpdf=None,angle_range=None)
    g_mk_list.append(g_mk[0][Ngrid:Ngrid+K])
    g_mk_list2.append(g_mk[0][0:Ngrid])
    g_mk_list3.append(g_mk[0][Ngrid+K])

np.savez('/output/graphs2/real2/real2_priormle.npz', aaa=g_mk_list3,aa=g_mk_list2,a=g_mk_list,c=points_inhomo,d=xxx,g=theta0,h=theta1,i=noise_var)

plt.figure()                
low=np.quantile(g_mk_list2, 0.025, axis=0)
plt.plot(xxx,low,'--',lw=1,alpha=0.6,c='tab:grey',label='2.5%-50%-97.5% Quantiles')
me=np.quantile(g_mk_list2, 0.5, axis=0)
plt.plot(xxx,me,'--',lw=1,alpha=0.6,c='tab:grey')
high=np.quantile(g_mk_list2, 0.975, axis=0)
plt.plot(xxx,high,'--',lw=1,alpha=0.6,c='tab:grey')
plt.plot(points_inhomo,np.zeros(len(points_inhomo)),linestyle='None', marker='|', color='k', markersize=10,label='events')
plt.ylim(0,6)
plt.ylabel(r'Poisson Intensity $\lambda(s)$')
plt.title('Earthquake dataset with RI-MAP fit',y=-0.2)
plt.legend()
plt.savefig("/output/graphs/real2/Priormle_median.pdf",bbox_inches='tight')
plt.show()

plt.figure()
plt.hist(g_mk_list3,bins=20)
plt.xlim(700,1200)
plt.title('Histogram of $\int_0^T\lambda(s)ds$', y=-0.2)
plt.savefig("/output/graphs/real2/Priormle_histogram.pdf",bbox_inches='tight')
plt.show()

plt.figure()
plt.plot(np.array(g_mk_list2)[:,50])
plt.title('latent GP at midpoint',y=-0.2)
plt.savefig("/output/graphs/real2/Priormle_traceplot1.pdf",bbox_inches='tight')
plt.show()

