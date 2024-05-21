#!/usr/bin/env python
# coding: utf-8

# In[1]:


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

theta0=37370.420361   
theta1=64.850410

loc1=[0.931481481481481,0.938888888888889,0.935185185185185,0.97962962962963,0.787037037037037,0.842592592592593,0.938888888888889,0.735185185185185,0.731481481481482,0.831481481481481,0.833333333333333,0.931481481481481,0.938888888888889,0.598148148148148,0.783333333333333,0.907407407407407,0.725925925925926,0.594444444444444,0.661111111111111,0.874074074074074,0.888888888888889,0.661111111111111,0.814814814814815,0.587037037037037,0.97962962962963,0.966666666666667,0.725925925925926,0.788888888888889,0.653703703703704,0.598148148148148,0.942592592592593,0.725925925925926,0.785185185185185,0.546296296296296,0.709259259259259,0.903703703703704,0.77962962962963,0.777777777777778,0.462962962962963,0.47037037037037,0.748148148148148,0.888888888888889,0.57962962962963,0.951851851851852,0.657407407407407,0.674074074074074,0.744444444444444,0.72962962962963,0.72037037037037,0.67962962962963,0.544444444444445,0.911111111111111,0.844444444444444,0.42962962962963,0.435185185185185,0.559259259259259,0.909259259259259,0.966666666666667,0.485185185185185,0.827777777777778,0.711111111111111,0.888888888888889,0.624074074074074,0.887037037037037,0.761111111111111,0.898148148148148,0.411111111111111,0.972222222222222,0.357407407407407,0.953703703703704,0.459259259259259,0.827777777777778,0.831481481481481,0.848148148148148,0.711111111111111,0.866666666666667,0.711111111111111,0.718518518518519,0.0388888888888889,0.716666666666667,0.8,0.612962962962963,0.711111111111111,0.764814814814815,0.740740740740741,0.792592592592593,0.162962962962963,0.707407407407407,0.755555555555556,0.237037037037037,0.768518518518519,0.674074074074074,0.655555555555556,0.711111111111111,0.140740740740741,0.694444444444444,0.681481481481481,0.138888888888889,0.657407407407407,0.677777777777778,0.155555555555556,0.390740740740741,0.474074074074074,0.398148148148148,0.218518518518519,0.122222222222222,0.205555555555556,0.264814814814815,0.296296296296296,0.303703703703704,0.312962962962963,0.25,0.637037037037037,0.272222222222222,0.731481481481482,0.625925925925926,0.311111111111111,0.062962962962963,0.594444444444444,0.618518518518519,0.601851851851852,0.222222222222222,0.627777777777778,0.235185185185185,0.598148148148148,0.281481481481481,0.218518518518519,0.248148148148148,0.231481481481481,0.475925925925926,0.538888888888889,0.198148148148148,0.277777777777778,0.531481481481482,0.205555555555556,0.52037037037037,0.185185185185185,0.177777777777778,0.403703703703704,0.475925925925926,0.407407407407407,0.233333333333333,0.187037037037037,0.411111111111111,0.431481481481482,0.22037037037037,0.437037037037037,0.242592592592593,0.383333333333333,0.237037037037037,0.37037037037037,0.42962962962963,0.383333333333333,0.407407407407407,0.427777777777778,0.331481481481481,0.342592592592593,0.335185185185185,0.316666666666667,0.325925925925926,0.292592592592593,0.303703703703704,0.316666666666667,0.0851851851851852,0.0796296296296296,0.0925925925925926,0.114814814814815,0.116666666666667,0.0944444444444444,0.1,0.0888888888888889,0.087037037037037,0.248148148148148,0.0592592592592593,0.0444444444444444,0.0759259259259259,0.281481481481481,0.244444444444444,0.0555555555555556,0.268518518518519,0.251851851851852,0.231481481481481,0.125925925925926,0.237037037037037,0.214814814814815,0.231481481481481,0.101851851851852,0.207407407407407,0.22962962962963,0.101851851851852,0.161111111111111,0.174074074074074,0.125925925925926,0.138888888888889,0.190740740740741]
loc2=[0.81767955801105,0.76427255985267,0.721915285451197,0.664825046040516,0.661141804788214,0.644567219152855,0.622467771639042,0.611418047882136,0.596685082872928,0.556169429097606,0.543278084714549,0.574585635359116,0.524861878453039,0.499079189686924,0.488029465930018,0.458563535911602,0.449355432780847,0.447513812154696,0.445672191528545,0.441988950276243,0.440147329650092,0.432780847145488,0.397790055248619,0.392265193370166,0.392265193370166,0.386740331491713,0.373848987108656,0.373848987108656,0.342541436464088,0.337016574585635,0.331491712707182,0.300184162062615,0.300184162062615,0.287292817679558,0.283609576427256,0.276243093922652,0.263351749539595,0.244935543278085,0.244935543278085,0.233885819521179,0.224677716390424,0.219152854511971,0.21731123388582,0.204419889502762,0.20073664825046,0.187845303867403,0.187845303867403,0.18232044198895,0.173112338858195,0.169429097605893,0.165745856353591,0.151012891344383,0.136279926335175,0.116022099447514,0.104972375690608,0.101289134438306,0.0994475138121547,0.0957642725598527,0.0902394106813996,0.0810313075506446,0.0791896869244936,0.0791896869244936,0.0736648250460405,0.0662983425414365,0.0644567219152855,0.0589318600368324,0.0534069981583794,0.0515653775322284,0.0497237569060774,0.0405156537753223,0.0294659300184162,0.0276243093922652,0.988950276243094,0.987108655616943,0.983425414364641,0.979742173112339,0.970534069981584,0.955801104972376,0.946593001841621,0.942909760589319,0.930018416206262,0.922651933701658,0.920810313075506,0.918968692449355,0.913443830570902,0.913443830570902,0.913443830570902,0.907918968692449,0.904235727440147,0.902394106813996,0.898710865561694,0.89134438305709,0.887661141804788,0.887661141804788,0.887661141804788,0.885819521178637,0.878453038674033,0.874769797421731,0.87292817679558,0.867403314917127,0.850828729281768,0.847145488029466,0.826887661141805,0.791896869244936,0.779005524861878,0.775322283609576,0.771639042357274,0.769797421731123,0.760589318600368,0.753222836095764,0.744014732965009,0.744014732965009,0.734806629834254,0.734806629834254,0.731123388581952,0.721915285451197,0.716390423572744,0.712707182320442,0.712707182320442,0.70902394106814,0.707182320441989,0.701657458563536,0.697974217311234,0.696132596685083,0.696132596685083,0.692449355432781,0.692449355432781,0.69060773480663,0.685082872928177,0.679558011049724,0.651933701657459,0.626151012891344,0.616942909760589,0.604051565377532,0.60036832412523,0.596685082872928,0.589318600368324,0.569060773480663,0.556169429097606,0.513812154696133,0.456721915285451,0.449355432780847,0.447513812154696,0.447513812154696,0.447513812154696,0.441988950276243,0.43646408839779,0.432780847145488,0.427255985267035,0.421731123388582,0.419889502762431,0.414364640883978,0.410681399631676,0.39963167587477,0.395948434622468,0.373848987108656,0.366482504604052,0.357274401473297,0.348066298342541,0.331491712707182,0.329650092081031,0.324125230202578,0.320441988950276,0.305709023941068,0.294659300184162,0.29097605893186,0.270718232044199,0.255985267034991,0.239410681399632,0.226519337016575,0.219152854511971,0.209944751381215,0.204419889502762,0.202578268876611,0.20073664825046,0.20073664825046,0.197053406998158,0.191528545119705,0.187845303867403,0.180478821362799,0.176795580110497,0.160220994475138,0.136279926335175,0.121546961325967,0.121546961325967,0.110497237569061,0.0994475138121547,0.0957642725598527,0.0920810313075506,0.0883977900552486,0.0755064456721915,0.0755064456721915,0.0736648250460405,0.0662983425414365,0.0662983425414365]


points_inhomo=np.concatenate((np.array(loc1).reshape(len(loc1),1),np.array(loc2).reshape(len(loc2),1)),axis=1)

T1=1
T2=1

def expo_quad_kernel(theta0,theta1,xn,xm): # 1,0.1
    return theta0*np.exp(-theta1/2*np.sum((xn - xm)**2))

def expo_quad_kernel2(theta0,theta1,xn,T1,T2): # 1,0.1
    return theta0*np.pi/2/theta1*(math.erf(np.sqrt(theta1/2)*(T1-xn[0]))+math.erf(np.sqrt(theta1/2)*xn[0]))*(math.erf(np.sqrt(theta1/2)*(T2-xn[1]))+math.erf(np.sqrt(theta1/2)*xn[1]))

def expo_quad_kernel3(theta0,theta1,T1,T2): # 1,0.1
    return theta0*2/theta1*(np.sqrt(np.pi*theta1/2)*T1*math.erf(np.sqrt(theta1/2)*T1)+ np.exp(-theta1/2*(T1**2))-1)*2/theta1*(np.sqrt(np.pi*theta1/2)*T2*math.erf(np.sqrt(theta1/2)*T2)+ np.exp(-theta1/2*(T2**2)) -1)

# ## simulation

# In[6]:


T1U=1
T2U=1
T1l=0
T2l=0
X=np.arange(T1l, T1U+0.05,0.05)
Y = np.arange(T2l, T2U+0.05,0.05)
X, Y = np.meshgrid(X, Y)
X, Y = X.flatten(),Y.flatten()
XY=np.column_stack([X,Y])

x4=np.concatenate((XY,points_inhomo))
# ## MCMC inference


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

# In[17]:
Ngrid=XY.shape[0]
K=len(points_inhomo)
N=K+1


#g_mk=chol_sample(mean=np.zeros(N), cov=cov_K, jit=noise_var)
g_mk=300+np.zeros((Ngrid+K+1))## 0--(Ngrid-1) Ngrid--(Ngrid+K-1) observe K--last unobserve +integral term

g_mk[Ngrid+K]=300*(T1U-T1l)*(T2U-T2l)
#print(g_mk)

Nfinal=Ngrid+N
cov_K=np.zeros((Nfinal,Nfinal))



nsim1=100000
nsim2=100000
for i in range(Nfinal):
        for j in range(i,Nfinal):
            if (j<(Nfinal-1)) and (i<(Nfinal-1)):
                cov_K[i][j]=expo_quad_kernel(theta0,theta1,x4[i],x4[j])
            if (j==(Nfinal-1)) and (i<(Nfinal-1)):
                cov_K[i][j]=expo_quad_kernel2(theta0,theta1,x4[i],T1,T2)
            if (j==(Nfinal-1)) and (i==(Nfinal-1)):
                cov_K[i][j]=expo_quad_kernel3(theta0,theta1,T1,T2)  
            if j!=i:
                cov_K[j][i]=cov_K[i][j]
noise_var=1e-8
cov_K_noise=cov_K+np.eye(Nfinal)*noise_var
min_eig=np.min(np.real(np.linalg.eigvals(cov_K_noise)))
print(min_eig)
cov_K_chol= np.linalg.cholesky(cov_K_noise)


for ite in range(nsim1):
    prior=chol_sample(mean=np.zeros(Nfinal), cov_chol=cov_K_chol)#nu
    g_mk,curloglike=elliptical_slice(g_mk.reshape(1,-1),prior,log_lik,pdf_params=[K],cur_lnpdf=None,angle_range=None)

        

g_mk_list=[]
g_mk_list2=[]
g_mk_list3=[]
for ite in range(nsim2):
    prior=chol_sample(mean=np.zeros(Nfinal), cov_chol=cov_K_chol)#nu
    g_mk,curloglike=elliptical_slice(g_mk.reshape(1,-1),prior,log_lik,pdf_params=[K],cur_lnpdf=None,angle_range=None)
    g_mk_list.append(g_mk[0][Ngrid:Ngrid+K])
    g_mk_list2.append(g_mk[0][0:Ngrid])
    g_mk_list3.append(g_mk[0][Ngrid+K])



low=np.quantile(g_mk_list2, 0.025, axis=0)
me=np.quantile(g_mk_list2, 0.5, axis=0)
high=np.quantile(g_mk_list2, 0.975, axis=0)


integral_mean=np.mean(g_mk_list3)
integral_sd=np.sqrt(np.cov(g_mk_list3))


plt.figure()
plt.plot(np.array(g_mk_list2)[:,220])
plt.title('latent GP at midpoint',y=-0.2)
plt.savefig("/output/graphs/real3/priormle_traceplot1.pdf",bbox_inches='tight')
plt.show()

plt.figure()
plt.hist(np.array(g_mk_list3))
plt.xlim(140,280)
plt.title('Histogram of $\int_S \lambda(s)ds$', y=-0.2)
plt.savefig("/output/graphs/real3/priormle_histogram.pdf",bbox_inches='tight')
plt.show()


X=np.arange(T1l, T1U+0.05,0.05)
Y = np.arange(T2l, T2U+0.05,0.05)
nn1=Y.shape[0]
nn2=X.shape[0]
X, Y = np.meshgrid(X, Y)
X, Y = X.flatten(),Y.flatten()
Xnew=X.reshape(nn1,nn2).T
Ynew=Y.reshape(nn1,nn2).T

pdf_pos_new2=me.reshape(nn1,nn2).T

fig1, ax2 = plt.subplots()
CS=ax2.contourf(Xnew, Ynew, pdf_pos_new2,cmap='Blues',alpha=0.5)
plt.contour(Xnew, Ynew, pdf_pos_new2,cmap='Blues',linewidths=3)
plt.plot(points_inhomo[:,0], points_inhomo[:,1],'.',color='black',alpha=0.5)

cbar = fig1.colorbar(CS)
cbar.ax.set_ylabel(r'Poisson Intensity $\lambda(s)$')



plt.title('Redwoods Data with RI-MAP Fit')
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
plt.savefig("/output/graphs2/real3/priomle_median.pdf",bbox_inches='tight')

np.savez('/output/graphs/real3/real3_priormle.npz', aaa=g_mk_list3,aa=g_mk_list2,a=g_mk_list,c=points_inhomo,n=integral_mean,o=integral_sd,s=theta0,t=theta1)






