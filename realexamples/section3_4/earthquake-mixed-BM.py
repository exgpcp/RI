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
points_inhomo=np.array(points_inhomo)
counts=[]
counts.append(((281 <= points_inhomo) & (points_inhomo < 288)).sum())
counts.append(((288 <= points_inhomo) & (points_inhomo < 295)).sum())
counts.append(((295 <= points_inhomo) & (points_inhomo < 302)).sum())
counts.append(((302 <= points_inhomo) & (points_inhomo < 309)).sum())
counts.append(((309 <= points_inhomo) & (points_inhomo < 316)).sum())
counts.append(((316 <= points_inhomo) & (points_inhomo < 323)).sum())
counts.append(((323 <= points_inhomo) & (points_inhomo < 330)).sum())
counts.append(((330 <= points_inhomo) & (points_inhomo < 337)).sum())
counts.append(((337 <= points_inhomo) & (points_inhomo < 344)).sum())
counts.append(((344 <= points_inhomo) & (points_inhomo < 351)).sum())
counts.append(((351 <= points_inhomo) & (points_inhomo < 358)).sum())
counts.append(((358 <= points_inhomo) & (points_inhomo < 365)).sum())

T=365
days=[0,281,288,295,302,309,316,323,330,337,344,351,358,365]
recurrent=points_inhomo[points_inhomo<281]

def expo_quad_kernel(xn,xm): # 1,0.1
    return min(xn,xm)

def expo_quad_kernel2(xn,t1,t2): # t1<t2   
    if (xn<=t1):
        return xn*(t2-t1)   
    elif (xn>t1 and xn<t2):
        return -xn**2/2-t1**2/2+xn*t2
    elif (xn>=t2):
        return (t2**2-t1**2)/2

def expo_quad_kernel3(t1,t2,t3,t4): # 1,0.1 #t1<t2, t3<t4
    if (t1==t3 and t2==t4 ):
        return t2**3/3+2*t1**3/3-t1**2*t2/2-t2*t1**2/2
    else:
        return (t4-t3)*(t2**2-t1**2)/2
        
def chol_sample(mean, cov_chol):
    return mean + cov_chol @ np.random.standard_normal(mean.size)

def log_lik(f,ns):
    #ns is # of Ngrids
    if np.prod(f>0)==0:
        return float('-inf') 
    else:
        return np.sum(np.log(f[:,Ngrid:Ngrid+ns])) -np.sum(f[:,(Ngrid+ns):Nfinal])+np.sum(  counts*np.log( f[:,(Ngrid+ns+1):Nfinal] ) )

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

K=len(counts)
N=len(recurrent)+1
Ngrid=100
x4=np.linspace(0,T,Ngrid+1)[1:Ngrid+1] 
xxx=np.concatenate([x4, points_inhomo])
Nfinal=Ngrid+K+N
g_mk=3+np.zeros(Nfinal)
g_mk[Ngrid+N-1]=3*281
nsim1=10000
nsim2=50000
cov_K=np.zeros((Nfinal,Nfinal))
theta=np.exp(4)
rem=1
alpha=.1
beta=.1

for i in range(Nfinal):
        for j in range(i,Nfinal):
            if (j<Ngrid+N-1) and (i<Ngrid+N-1):
                cov_K[i][j]=expo_quad_kernel(xxx[i],xxx[j])
            if (j>(Ngrid+N-2)) and (i<(Ngrid+N-1)):
                cov_K[i][j]=expo_quad_kernel2(xxx[i],    days[  j-Ngrid-N+1  ], days[  j-Ngrid-N+2           ]    )
            if  (i>(Ngrid+N-2)):
                cov_K[i][j]=expo_quad_kernel3(days[i-Ngrid-N+1],days[i-Ngrid-N+2],  days[j-Ngrid-N+1],days[j-Ngrid-N+2 ]          )  
            if j!=i:
                cov_K[j][i]=cov_K[i][j]

noise_var=1e-7
cov_K_inv=np.linalg.inv(cov_K+np.eye(Nfinal)*noise_var)
ones_vec=np.ones(Nfinal)
ones_vec[Ngrid+N-1]=281
ones_vec[(Ngrid+N):Nfinal]=7
ones_vec=ones_vec.reshape(Nfinal,1)
cov_K_mod_inv=cov_K_inv-cov_K_inv@ones_vec@ones_vec.transpose()@cov_K_inv/(ones_vec.transpose()@cov_K_inv@ones_vec)
cov_K_mod_inv_add=cov_K_mod_inv+np.eye(Nfinal)*noise_var
cov_K_mod=np.linalg.inv(cov_K_mod_inv_add)
cov_K_mod_chol=np.linalg.cholesky(cov_K_mod) 

for ite in range(nsim1):
    cov_K_chol_final=cov_K_mod_chol/np.sqrt(theta)
    prior=chol_sample(mean=np.zeros(Nfinal), cov_chol=cov_K_chol_final)#nu
    g_mk,curloglike=elliptical_slice(g_mk.reshape(1,-1),prior,log_lik,pdf_params=[N-1],cur_lnpdf=None,angle_range=None)
    alpha_pos=alpha+Nfinal/2
    beta_pos=beta+1/2*g_mk[0]@cov_K_mod_inv@g_mk[0]
    theta=np.random.gamma(alpha_pos, 1/beta_pos,1)
      
g_mk_list2=[]
g_mk_list3=[]
theta_list=[]
for ite in range(nsim2):
    cov_K_chol_final=cov_K_mod_chol/np.sqrt(theta)
    prior=chol_sample(mean=np.zeros(Nfinal), cov_chol=cov_K_chol_final)#nu
    g_mk,curloglike=elliptical_slice(g_mk.reshape(1,-1),prior,log_lik,pdf_params=[N-1],cur_lnpdf=None,angle_range=None)
    g_mk_list2.append(g_mk[0][0:Ngrid])
    g_mk_list3.append(g_mk[0][Ngrid:Nfinal])
    alpha_pos=alpha+Nfinal/2
    beta_pos=beta+1/2*g_mk[0]@cov_K_mod_inv@g_mk[0]
    theta=np.random.gamma(alpha_pos, 1/beta_pos,1)
    theta_list.append(theta)   

data=np.load('/output/graphs2/real2/BM.npz')
g_mk_list3_rec=data['aaa']
g_mk_list2_rec=data['aa']
g_mk_list_rec=data['a']
xxx_rec=data['d']
theta_list_rec=data['s']
low_rec=np.quantile(g_mk_list2_rec, 0.025, axis=0)
high_rec=np.quantile(g_mk_list2_rec, 0.975, axis=0)
me_rec=np.quantile(g_mk_list2_rec, 0.5, axis=0)

low=np.quantile(g_mk_list2, 0.025, axis=0)
high=np.quantile(g_mk_list2, 0.975, axis=0)
me=np.quantile(g_mk_list2, 0.5, axis=0)

fig, ax = plt.subplots()
ax.plot(x4,me,'-',lw=1,alpha=0.6,c='r')
ax.plot(x4,low,'-',lw=1,alpha=0.6,c='r',label='Mixed: 2.5%-50%-97.5% Quantiles')
ax.plot(x4,high,'-',lw=1,alpha=0.6,c='r')
ax.fill_between(x4, low, high, color='r', alpha=.05)
ax.plot(xxx_rec,me_rec,'--',lw=1,alpha=0.6,c='b')
ax.plot(xxx_rec,low_rec,'--',lw=1,alpha=0.6,c='b',label='Recurrent: 2.5%-50%-97.5% Quantiles')
ax.plot(xxx_rec,high_rec,'--',lw=1,alpha=0.6,c='b')
ax.fill_between(xxx_rec, low_rec, high_rec, color='b', alpha=.05)

plt.plot(points_inhomo,np.zeros(len(points_inhomo)),linestyle='None', marker='|', color='k', markersize=10,label='events')
plt.ylim(0,6)
plt.xlabel('Time (days)')
plt.ylabel(r'Poisson Intensity $\lambda(s)$')
plt.title('Mixed and Recurrent Earthquake Datasets with RI-BM Fit')
plt.bar(days[1:13], height = np.array(counts)/7, width = 6 , label='average weekly counts')
plt.legend(loc=9)
plt.savefig("/output/graphs2/real5/BM_median.pdf",bbox_inches='tight')
plt.show()

np.savez('/output/graphs2/real5/BM.npz', aaa=g_mk_list3,aa=g_mk_list2,c=points_inhomo,d=x4,s=theta_list)

