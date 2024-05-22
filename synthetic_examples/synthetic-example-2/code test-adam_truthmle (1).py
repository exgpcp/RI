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

import sys
#mvn.mvnun(np.ones(D)*float("-Inf"),np.zeros(D),np.zeros(D),cov_K_noise)[0]


# In[15]:


l2_dist1_1=np.zeros((100,1))
l2_dist1_2=np.zeros((100,1))
width1_1=np.zeros((100,1))
width1_2=np.zeros((100,1))
coverage1_1=np.zeros((100,1))
coverage1_2=np.zeros((100,1))
timerun1_1=np.zeros((100,1))
timerun1_2=np.zeros((100,1))


# In[16]:


for jjj in range(100):
    print(jjj)
    data = np.load('/scratch/groups/juliapr/output_Bingjing/simulation/adam1_1/truthmle/syn1/syn1_truthmle_'+str(jjj+1)+'.npz')
    c=math.ceil((jjj+1)/100)
    
    #,aa=g_mk_list2,a=g_mk_list,c=points_inhomo,d=xxx,
    #e=theta0,f=theta1,g=measure_sup,h=noise_var,i=coverage1,j=coverage2,k=l2_dist1,l=l2_dist2,m=width1,n=width2,o=timerun1,p=timerun2)
    
    
    
    
    g_mk_list2=data['aa']
    g_mk_list=data['a']
    points_inhomo=data['c']
    xxx=data['d']
    
    
    coverage1=data['i']
    theta0=data['e']
    theta1=data['f']
    measure_sup=data['g']
    noise_var=data['h']
    coverage2=data['j']
    l2_dist1=data['k']
    l2_dist2=data['l']
    width1=data['m']
    width2=data['n']
    
    timerun1=data['o']
    timerun2=data['p']


    
    l2_dist1_1[jjj]=l2_dist1
    coverage1_1[jjj]=coverage1
    width1_1[jjj]=width1


    
    l2_dist1_2[jjj]=l2_dist2
    coverage1_2[jjj]=coverage2
    width1_2[jjj]=width2
    timerun1_1[jjj]=timerun1
    timerun1_2[jjj]=timerun2
    
   


# In[5]:


len(g_mk_list)


# In[17]:


k=0
print([np.round(np.quantile(l2_dist1_1[k:(k+100)], 0.5, axis=0),2)[0],
      np.round(np.quantile(l2_dist1_1[k:(k+100)], 0.025, axis=0),2)[0],
       np.round(np.quantile(l2_dist1_1[k:(k+100)], 0.25, axis=0),2)[0],
       np.round(np.quantile(l2_dist1_1[k:(k+100)], 0.75, axis=0),2)[0],
       np.round(np.quantile(l2_dist1_1[k:(k+100)], 0.975, axis=0),2)[0]])
       
       

print([np.round(np.quantile(100*coverage1_1[k:(k+100)], 0.5, axis=0),0)[0],
       np.round(np.quantile(100*coverage1_1[k:(k+100)], 0.025, axis=0),0)[0],
       np.round(np.quantile(100*coverage1_1[k:(k+100)], 0.25, axis=0),0)[0],
       np.round(np.quantile(100*coverage1_1[k:(k+100)], 0.75, axis=0),0)[0],
       np.round(np.quantile(100*coverage1_1[k:(k+100)], 0.975, axis=0),0)[0]])
      
    
print([np.round(np.quantile(width1_1[k:(k+100)], 0.5, axis=0),2)[0],
      np.round(np.quantile(width1_1[k:(k+100)], 0.025, axis=0),2)[0],
       np.round(np.quantile(width1_1[k:(k+100)], 0.25, axis=0),2)[0],
       np.round(np.quantile(width1_1[k:(k+100)], 0.75, axis=0),2)[0],
       np.round(np.quantile(width1_1[k:(k+100)], 0.975, axis=0),2)[0]])
       
      


print([np.round(np.quantile(l2_dist1_2[k:(k+100)], 0.5, axis=0),2)[0],
      np.round(np.quantile(l2_dist1_2[k:(k+100)], 0.025, axis=0),2)[0],
      np.round(np.quantile(l2_dist1_2[k:(k+100)], 0.25, axis=0),2)[0],
      np.round(np.quantile(l2_dist1_2[k:(k+100)], 0.75, axis=0),2)[0],
      np.round(np.quantile(l2_dist1_2[k:(k+100)], 0.975, axis=0),2)[0]])
      
      
      
      
      
    
print([np.round(np.quantile(100*coverage1_2[k:(k+100)], 0.5, axis=0),0)[0],
       np.round(np.quantile(100*coverage1_2[k:(k+100)], 0.025, axis=0),0)[0],
       np.round(np.quantile(100*coverage1_2[k:(k+100)], 0.25, axis=0),0)[0],
       np.round(np.quantile(100*coverage1_2[k:(k+100)], 0.75, axis=0),0)[0],
       np.round(np.quantile(100*coverage1_2[k:(k+100)], 0.975, axis=0),0)[0]])
      
print([np.round(np.quantile(width1_2[k:(k+100)], 0.5, axis=0),2)[0],
      np.round(np.quantile(width1_2[k:(k+100)], 0.025, axis=0),2)[0],
       np.round(np.quantile(width1_2[k:(k+100)], 0.25, axis=0),2)[0],
       np.round(np.quantile(width1_2[k:(k+100)], 0.75, axis=0),2)[0],
       np.round(np.quantile(width1_2[k:(k+100)], 0.975, axis=0),2)[0]])
       

print( [np.round(np.mean(timerun1_2[k:(k+100)]),2),  np.round(np.std(timerun1_2[k:(k+100)]) ,2) ])


# In[10]:


l2_dist2_1=np.zeros((300,1))
l2_dist2_2=np.zeros((300,1))
coverage2_1=np.zeros((300,1))
coverage2_2=np.zeros((300,1))
timerun2_1=np.zeros((300,1))
timerun2_2=np.zeros((300,1))
width2_1=np.zeros((300,1))
width2_2=np.zeros((300,1))


# In[13]:


for jjj in range(100):
    print(jjj)
    data = np.load('/scratch/groups/juliapr/output_Bingjing/simulation/adam1_1/truthmle/syn2/syn2_truthmle_'+str(jjj+301)+'.npz')
    c=math.ceil((jjj+1)/100)
    
    #,aa=g_mk_list2,a=g_mk_list,c=points_inhomo,d=xxx,
    #e=theta0,f=theta1,g=measure_sup,h=noise_var,i=coverage1,j=coverage2,k=l2_dist1,l=l2_dist2,m=width1,n=width2,o=timerun1,p=timerun2)
    
    
    
    
    g_mk_list2=data['aa']
    g_mk_list=data['a']
    points_inhomo=data['c']
    xxx=data['d']
    
    
    coverage1=data['i']
    theta0=data['e']
    theta1=data['f']
    measure_sup=data['g']
    noise_var=data['h']
    coverage2=data['j']
    l2_dist1=data['k']
    l2_dist2=data['l']
    width1=data['m']
    width2=data['n']
    
    timerun1=data['o']
    timerun2=data['p']


    
    l2_dist2_1[jjj]=l2_dist1
    coverage2_1[jjj]=coverage1
    width2_1[jjj]=width1


    
    l2_dist2_2[jjj]=l2_dist2
    coverage2_2[jjj]=coverage2
    width2_2[jjj]=width2
    timerun2_1[jjj]=timerun1
    timerun2_2[jjj]=timerun2
    
   


# In[14]:


k=0
print([np.round(np.quantile(l2_dist2_1[k:(k+100)], 0.5, axis=0),2)[0],
       np.round(np.quantile(l2_dist2_1[k:(k+100)], 0.025, axis=0),2)[0],
      np.round(np.quantile(l2_dist2_1[k:(k+100)], 0.25, axis=0),2)[0],
       np.round(np.quantile(l2_dist2_1[k:(k+100)], 0.75, axis=0),2)[0],
       np.round(np.quantile(l2_dist2_1[k:(k+100)], 0.975, axis=0),2)[0]
      ])


print([np.round(np.quantile(coverage2_1[k:(k+100)], 0.5, axis=0),2)[0],
     np.round(np.quantile(coverage2_1[k:(k+100)], 0.025, axis=0),2)[0],
      np.round(np.quantile(coverage2_1[k:(k+100)], 0.25, axis=0),2)[0],
      np.round(np.quantile(coverage2_1[k:(k+100)], 0.75, axis=0),2)[0],
      np.round(np.quantile(coverage2_1[k:(k+100)], 0.975, axis=0),2)[0]])


print([np.round(np.quantile(width2_1[k:(k+100)], 0.5, axis=0),2)[0],
      np.round(np.quantile(width2_1[k:(k+100)], 0.025, axis=0),2)[0],
       np.round(np.quantile(width2_1[k:(k+100)], 0.25, axis=0),2)[0],
       np.round(np.quantile(width2_1[k:(k+100)], 0.75, axis=0),2)[0],
       np.round(np.quantile(width2_1[k:(k+100)], 0.975, axis=0),2)[0]])
       

print([np.round(np.quantile(l2_dist2_2[k:(k+100)], 0.5, axis=0),2)[0],
      np.round(np.quantile(l2_dist2_2[k:(k+100)], 0.025, axis=0),2)[0],
      np.round(np.quantile(l2_dist2_2[k:(k+100)], 0.25, axis=0),2)[0],
      np.round(np.quantile(l2_dist2_2[k:(k+100)], 0.75, axis=0),2)[0],
      np.round(np.quantile(l2_dist2_2[k:(k+100)], 0.975, axis=0),2)[0]])


print([np.round(np.quantile(coverage2_2[k:(k+100)], 0.5, axis=0),2)[0],
      np.round(np.quantile(coverage2_2[k:(k+100)], 0.025, axis=0),2)[0],
    np.round(np.quantile(coverage2_2[k:(k+100)], 0.25, axis=0),2)[0],
      np.round(np.quantile(coverage2_2[k:(k+100)], 0.75, axis=0),2)[0],
np.round(np.quantile(coverage2_2[k:(k+100)], 0.975, axis=0),2)[0]])


print([np.round(np.quantile(width2_2[k:(k+100)], 0.5, axis=0),2)[0],
      np.round(np.quantile(width2_2[k:(k+100)], 0.025, axis=0),2)[0],
       np.round(np.quantile(width2_2[k:(k+100)], 0.25, axis=0),2)[0],
       np.round(np.quantile(width2_2[k:(k+100)], 0.75, axis=0),2)[0],
       np.round(np.quantile(width2_2[k:(k+100)], 0.975, axis=0),2)[0]])
       

print( [np.round(np.mean(timerun2_2[k:(k+100)]),2),  np.round(np.std(timerun2_2[k:(k+100)]) ,2) ])

#print(np.round(np.mean(timerun2_2[k:(k+100)]/5),2)[0],
      
      
#print(np.round(np.quantile(timerun2_2[k:(k+100)]/5, 0.975, axis=0),2)[0])


# In[8]:


def inten3(t):
    if t<=25:
        return 2+1/25*t
    if t<=50:
        return 5-2/25*t
    if t<=75:
        return 3/50*t-2
    return 1/50*t+1
inten3v = np.vectorize(inten3)


# In[9]:


l2_dist3_1=np.zeros((300,1))
l2_dist3_2=np.zeros((300,1))
coverage3_1=np.zeros((300,1))
coverage3_2=np.zeros((300,1))
timerun3_1=np.zeros((300,1))
timerun3_2=np.zeros((300,1))
width3_1=np.zeros((300,1))
width3_2=np.zeros((300,1))


# In[10]:


for jjj in range(300):
    print(jjj)
    data = np.load('/scratch/groups/juliapr/output_Bingjing/simulation/mymethodfinal3/ESSsyn3_priormle'+str(jjj+601)+'.npz')
    c=math.ceil((jjj+1)/100)
    g_mk_list3=data['a']
    g_mk_list3=data['aaa']
    g_mk_list2=data['aa']
    g_mk_list=data['a']
    points_inhomo=data['c']
    xxx=data['d']
    x=data['b']
    intensity=data['e']
    coverage=data['f']
    theta0=data['g']
    theta1=data['h']
    noise_var=data['i']
    l_2norm=data['j']
    integral_truth=data['m']
    integral_mean=data['n']
    integral_sd=data['o']
    timerun1=data['p']
    timerun2=data['q']


    low=np.quantile(g_mk_list, 0.025, axis=0)
    high=np.quantile(g_mk_list, 0.975, axis=0)
    med=np.quantile(g_mk_list, 0.5, axis=0)
    truth=inten3v(points_inhomo)
    truth=c*truth
    l2_dist3_1[jjj]=sum((np.array(med).squeeze()-truth)**2)
    width3_1[jjj]=sum(high-low)/len(points_inhomo)
    coverage3_1[jjj]=np.sum((truth>=low.squeeze()) * (truth<=high.squeeze()))/len(points_inhomo)


    low=np.quantile(g_mk_list2, 0.025, axis=0)
    high=np.quantile(g_mk_list2, 0.975, axis=0)
    med=np.quantile(g_mk_list2, 0.5, axis=0)
    truth=inten3v(xxx)
    truth=c*truth
    l2_dist3_2[jjj]=sum((np.array(med).squeeze()-truth)**2)
    width3_2[jjj]=sum(high-low)/len(xxx)
    coverage3_2[jjj]=np.sum((truth>=low.squeeze()) * (truth<=high.squeeze()))/len(xxx)
    timerun3_1[jjj]=timerun1
    timerun3_2[jjj]=timerun2
    
    


# In[39]:


k=0

print([np.round(np.quantile(l2_dist3_1[k:(k+100)], 0.5, axis=0),2)[0],
       np.round(np.quantile(l2_dist3_1[k:(k+100)], 0.025, axis=0),2)[0],
      np.round(np.quantile(l2_dist3_1[k:(k+100)], 0.25, axis=0),2)[0],
       np.round(np.quantile(l2_dist3_1[k:(k+100)], 0.75, axis=0),2)[0],
       np.round(np.quantile(l2_dist3_1[k:(k+100)], 0.975, axis=0),2)[0]
      ])


print([np.round(np.quantile(100*coverage3_1[k:(k+100)], 0.5, axis=0),0)[0],
     np.round(np.quantile(100*coverage3_1[k:(k+100)], 0.025, axis=0),0)[0],
      np.round(np.quantile(100*coverage3_1[k:(k+100)], 0.25, axis=0),0)[0],
      np.round(np.quantile(100*coverage3_1[k:(k+100)], 0.75, axis=0),0)[0],
      np.round(np.quantile(100*coverage3_1[k:(k+100)], 0.975, axis=0),0)[0]])


print([np.round(np.quantile(width3_1[k:(k+100)], 0.5, axis=0),2)[0],
      np.round(np.quantile(width3_1[k:(k+100)], 0.025, axis=0),2)[0],
       np.round(np.quantile(width3_1[k:(k+100)], 0.25, axis=0),2)[0],
       np.round(np.quantile(width3_1[k:(k+100)], 0.75, axis=0),2)[0],
       np.round(np.quantile(width3_1[k:(k+100)], 0.975, axis=0),2)[0]])
       


print([np.round(np.quantile(l2_dist3_2[k:(k+100)], 0.5, axis=0),2)[0],
      np.round(np.quantile(l2_dist3_2[k:(k+100)], 0.025, axis=0),2)[0],
      np.round(np.quantile(l2_dist3_2[k:(k+100)], 0.25, axis=0),2)[0],
      np.round(np.quantile(l2_dist3_2[k:(k+100)], 0.75, axis=0),2)[0],
      np.round(np.quantile(l2_dist3_2[k:(k+100)], 0.975, axis=0),2)[0]])


print([np.round(np.quantile(100*coverage3_2[k:(k+100)], 0.5, axis=0),0)[0],
      np.round(np.quantile(100*coverage3_2[k:(k+100)], 0.025, axis=0),0)[0],
    np.round(np.quantile(100*coverage3_2[k:(k+100)], 0.25, axis=0),0)[0],
      np.round(np.quantile(100*coverage3_2[k:(k+100)], 0.75, axis=0),0)[0],
np.round(np.quantile(100*coverage3_2[k:(k+100)], 0.975, axis=0),0)[0]])

print([np.round(np.quantile(width3_2[k:(k+100)], 0.5, axis=0),2)[0],
      np.round(np.quantile(width3_2[k:(k+100)], 0.025, axis=0),2)[0],
       np.round(np.quantile(width3_2[k:(k+100)], 0.25, axis=0),2)[0],
       np.round(np.quantile(width3_2[k:(k+100)], 0.75, axis=0),2)[0],
       np.round(np.quantile(width3_2[k:(k+100)], 0.975, axis=0),2)[0]])
       
print( [np.round(np.mean(timerun3_2[k:(k+100)]/5),2),  np.round(np.std(timerun3_2[k:(k+100)]/5) ,2) ])


# In[46]:


np.savez('/scratch/groups/juliapr/output_Bingjing/simulation/mymethod2/analysis.npz', a=l2_dist1_1,b=l2_dist1_2,c=coverage1_1,d=
coverage1_2,
e=timerun1_1,
f=timerun1_2,
aa=l2_dist2_1,bb=l2_dist2_2,cc=coverage2_1,dd=
coverage2_2,
ee=timerun2_1,
ff=timerun2_2,
aaa=l2_dist3_1,bbb=l2_dist3_2,ccc=coverage3_1,ddd=
coverage3_2,
eee=timerun3_1,
fff=timerun3_2)
        
        
        
        
        #aaa=g_mk_list3,aa=g_mk_list2,a=g_mk_list,c=points_inhomo,d=xxx,b=x,e=intensity,f=coverage,g=theta0,h=theta1,i=noise_var,j=l_2norm,k=lp,l=lpdist,m=integral_truth,n=integral_mean,o=integral_sd)

