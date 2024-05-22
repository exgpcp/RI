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


# In[2]:


l2_dist1_1=np.zeros((300,1))
l2_dist1_2=np.zeros((300,1))
coverage1_1=np.zeros((300,1))
coverage1_2=np.zeros((300,1))
timerun1=np.zeros((300,1))
width1_1=np.zeros((300,1))
width1_2=np.zeros((300,1))


# In[3]:


for jjj in range(300):
    print(jjj)
    data = np.load('/scratch/groups/juliapr/output_Bingjing/simulation/mymethod_BM1_5_4/ESSsyn1_BM1_'+str(jjj+1)+'.npz')
    c=math.ceil((jjj+1)/100)
    g_mk_list3=data['aaa']
    g_mk_list2=data['aa']
    g_mk_list=data['a']
    points_inhomo=data['c']
    xxx=data['d']
    x=data['b']
    intensity=data['e']
    coverage1=data['f']
    coverage2=data['ff']
    noise_var=data['i']
    l2_dist1=data['j']
    l2_dist2=data['jj']
    integral_truth=data['m']
    integral_mean=data['n']
    integral_sd=data['o']
    timerun=data['p']
    width1=data['q']
    width2=data['r']
    theta_list=data['s']
    l2_dist1_1[jjj]=l2_dist1
    coverage1_1[jjj]=coverage1
    width1_1[jjj]=width1


    
    l2_dist1_2[jjj]=l2_dist2
    coverage1_2[jjj]=coverage2
    timerun1[jjj]=timerun
    width1_2[jjj]=width2
   


# In[38]:


len(g_mk_list)


# In[36]:


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
      
      

print( [np.round(np.mean(timerun1[k:(k+100)]/6),2),  np.round(np.std(timerun1[k:(k+100)]/6) ,2) ])

print([np.round(np.quantile(width1_1[k:(k+100)], 0.5, axis=0),2)[0],
      np.round(np.quantile(width1_1[k:(k+100)], 0.025, axis=0),2)[0],
       np.round(np.quantile(width1_1[k:(k+100)], 0.25, axis=0),2)[0],
       np.round(np.quantile(width1_1[k:(k+100)], 0.75, axis=0),2)[0],
       np.round(np.quantile(width1_1[k:(k+100)], 0.975, axis=0),2)[0]])     

print([np.round(np.quantile(width1_2[k:(k+100)], 0.5, axis=0),2)[0],
      np.round(np.quantile(width1_2[k:(k+100)], 0.025, axis=0),2)[0],
       np.round(np.quantile(width1_2[k:(k+100)], 0.25, axis=0),2)[0],
       np.round(np.quantile(width1_2[k:(k+100)], 0.75, axis=0),2)[0],
       np.round(np.quantile(width1_2[k:(k+100)], 0.975, axis=0),2)[0]])     
      


# In[7]:


#########syn2

l2_dist2_1=np.zeros((300,1))
l2_dist2_2=np.zeros((300,1))
coverage2_1=np.zeros((300,1))
coverage2_2=np.zeros((300,1))

timerun2=np.zeros((300,1))
width2_1=np.zeros((300,1))
width2_2=np.zeros((300,1))
for jjj in range(300):
    print(jjj)
    data = np.load('/scratch/groups/juliapr/output_Bingjing/simulation/mymethod_BM1_5_4/ESSsyn2_BM1_'+str(jjj+301)+'.npz')
    c=math.ceil((jjj+1)/100)
    g_mk_list3=data['aaa']
    g_mk_list2=data['aa']
    g_mk_list=data['a']
    points_inhomo=data['c']
    xxx=data['d']
    x=data['b']
    intensity=data['e']
    coverage1=data['f']
    coverage2=data['ff']
    noise_var=data['i']
    l2_dist1=data['j']
    l2_dist2=data['jj']
    integral_truth=data['m']
    integral_mean=data['n']
    integral_sd=data['o']

    
    
    
    timerun=data['p']
    width1=data['q']
    width2=data['r']
    theta_list_2=data['s']
    l2_dist2_1[jjj]=l2_dist1
    coverage2_1[jjj]=coverage1
    width2_1[jjj]=width1


    
    l2_dist2_2[jjj]=l2_dist2
    coverage2_2[jjj]=coverage2
    timerun2[jjj]=timerun
    width2_2[jjj]=width2


# In[41]:


k=200
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

print([np.round(np.quantile(l2_dist2_2[k:(k+100)], 0.5, axis=0),2)[0],
      np.round(np.quantile(l2_dist2_2[k:(k+100)], 0.025, axis=0),2)[0],
      np.round(np.quantile(l2_dist2_2[k:(k+100)], 0.25, axis=0),2)[0],
      np.round(np.quantile(l2_dist2_2[k:(k+100)], 0.75, axis=0),2)[0],
      np.round(np.quantile(l2_dist2_2[k:(k+100)], 0.95, axis=0),2)[0]])


print([np.round(np.quantile(coverage2_2[k:(k+100)], 0.5, axis=0),2)[0],
      np.round(np.quantile(coverage2_2[k:(k+100)], 0.025, axis=0),2)[0],
    np.round(np.quantile(coverage2_2[k:(k+100)], 0.25, axis=0),2)[0],
      np.round(np.quantile(coverage2_2[k:(k+100)], 0.75, axis=0),2)[0],
np.round(np.quantile(coverage2_2[k:(k+100)], 0.95, axis=0),2)[0]])


print( [np.round(np.mean(timerun2[k:(k+100)]/6),2),  np.round(np.std(timerun2[k:(k+100)]/6) ,2) ])

print([np.round(np.quantile(width2_1[k:(k+100)], 0.5, axis=0),2)[0],
      np.round(np.quantile(width2_1[k:(k+100)], 0.025, axis=0),2)[0],
       np.round(np.quantile(width2_1[k:(k+100)], 0.25, axis=0),2)[0],
       np.round(np.quantile(width2_1[k:(k+100)], 0.75, axis=0),2)[0],
       np.round(np.quantile(width2_1[k:(k+100)], 0.975, axis=0),2)[0]])     

print([np.round(np.quantile(width2_2[k:(k+100)], 0.5, axis=0),2)[0],
      np.round(np.quantile(width2_2[k:(k+100)], 0.025, axis=0),2)[0],
       np.round(np.quantile(width2_2[k:(k+100)], 0.25, axis=0),2)[0],
       np.round(np.quantile(width2_2[k:(k+100)], 0.75, axis=0),2)[0],
       np.round(np.quantile(width2_2[k:(k+100)], 0.975, axis=0),2)[0]])     
      


#print(np.round(np.mean(timerun2_2[k:(k+100)]/5),2)[0],
      
      
#print(np.round(np.quantile(timerun2_2[k:(k+100)]/5, 0.975, axis=0),2)[0])


# In[25]:


plt.plot(theta_list_2[10000:])


# In[26]:


np.where(coverage2_1[100:200]<0.8)


# In[12]:


#######syn3


l2_dist3_1=np.zeros((300,1))
l2_dist3_2=np.zeros((300,1))
coverage3_1=np.zeros((300,1))
coverage3_2=np.zeros((300,1))



timerun3=np.zeros((300,1))
width3_1=np.zeros((300,1))
width3_2=np.zeros((300,1))


# In[13]:


for jjj in range(300):
    print(jjj)
    #data = np.load('/scratch/groups/juliapr/output_Bingjing/simulation/mymethod_BM3/ESSsyn3_BM'+str(jjj+601)+'.npz')
    data = np.load('/scratch/groups/juliapr/output_Bingjing/simulation/mymethod_BM1_5_4/ESSsyn3_BM1'+str(jjj+601)+'.npz')
    c=math.ceil((jjj+1)/100)
  
    
    g_mk_list3=data['aaa']
    g_mk_list2=data['aa']
    g_mk_list=data['a']
    points_inhomo=data['c']
    xxx=data['d']
    x=data['b']
    intensity=data['e']
    coverage1=data['f']
    coverage2=data['ff']
    noise_var=data['i']
    l2_dist1=data['j']
    l2_dist2=data['jj']
    integral_truth=data['m']
    integral_mean=data['n']
    integral_sd=data['o']
    
    
    timerun=data['p']
    width1=data['q']
    width2=data['r']
    theta_list_3=data['s']
    l2_dist3_1[jjj]=l2_dist1
    coverage3_1[jjj]=coverage1
    width3_1[jjj]=width1


    
    l2_dist3_2[jjj]=l2_dist2
    coverage3_2[jjj]=coverage2
    timerun3[jjj]=timerun
    width3_2[jjj]=width2
    


# In[26]:


np.where(l2_dist3_1[200:300]>800)


# In[29]:


np.where(l2_dist3_2[200:300]>118)


# In[26]:


print(np.where(coverage3_1[200:300]<0.88))
print(np.where(coverage3_2[200:300]<0.88))


# In[25]:


print(l2_dist3_1[130])
print( l2_dist3_2[130] )
print(coverage3_1[130])
print(coverage3_2[130])


# In[44]:


k=200

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


print( [np.round(np.mean(timerun3[k:(k+100)]/6),2),  np.round(np.std(timerun3[k:(k+100)]/6) ,2) ])


print([np.round(np.quantile(width3_1[k:(k+100)], 0.5, axis=0),2)[0],
      np.round(np.quantile(width3_1[k:(k+100)], 0.025, axis=0),2)[0],
       np.round(np.quantile(width3_1[k:(k+100)], 0.25, axis=0),2)[0],
       np.round(np.quantile(width3_1[k:(k+100)], 0.75, axis=0),2)[0],
       np.round(np.quantile(width3_1[k:(k+100)], 0.975, axis=0),2)[0]])     

print([np.round(np.quantile(width3_2[k:(k+100)], 0.5, axis=0),2)[0],
      np.round(np.quantile(width3_2[k:(k+100)], 0.025, axis=0),2)[0],
       np.round(np.quantile(width3_2[k:(k+100)], 0.25, axis=0),2)[0],
       np.round(np.quantile(width3_2[k:(k+100)], 0.75, axis=0),2)[0],
       np.round(np.quantile(width3_2[k:(k+100)], 0.975, axis=0),2)[0]])     


# In[ ]:





# In[21]:


jjj=200

data = np.load('/scratch/groups/juliapr/output_Bingjing/simulation/mymethod2/ESSsyn3_priormle'+str(jjj+601)+'.npz')
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
truth=inten3v(xxx)
truth=c*truth


low=np.quantile(g_mk_list2, 0.025, axis=0)
high=np.quantile(g_mk_list2, 0.975, axis=0)
med=np.quantile(g_mk_list2, 0.5, axis=0)
truth=inten3v(xxx)
truth=c*truth


    

plt.figure(1)                # the first figure
plt.subplot(111)             # the first subplot in the first figure

#plt.scatter(points_inhomo,np.average(g_mk_list_ess,axis=0))
plt.plot(xxx,low,'--',lw=1,alpha=0.6,c='tab:grey',label='2.5%-50%-97.5% Quantiles')

plt.plot(xxx,med,'--',lw=1,alpha=0.6,c='tab:grey')

# Calculate middle values

#plt.plot(xxx,np.average(g_mk_list2,axis=0),'g-',lw=1,alpha=0.6)
#plt.plot(points_inhomo,sigmoid(np.average(MM,axis=0)),'g-')

plt.plot(xxx,high,'--',lw=1,alpha=0.6,c='tab:grey')
plt.plot(xxx,truth,'k-',lw=1,alpha=0.6,label='Truth')
plt.plot(points_inhomo,np.zeros(len(points_inhomo)),linestyle='None', marker='|', color='k', markersize=10,label='events')
#plt.xlabel("t",y=-0.3)
plt.ylabel(r'Poisson Intensity $\lambda(s)$')
plt.title('The first synthetic data and our proposed model fit',y=-0.15)
plt.legend()


# In[22]:


np.savez('/scratch/groups/juliapr/output_Bingjing/simulation/mymethod2truth/analysis.npz', a=l2_dist1_1,b=l2_dist1_2,c=coverage1_1,d=
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
        
        


# In[31]:


k=200
print(np.round(np.quantile(l2_dist1_1[k:(k+100)], 0.5, axis=0),2)[0])
print(np.round(np.quantile(l2_dist1_1[k:(k+100)], 0.025, axis=0),2)[0])
print(np.round(np.quantile(l2_dist1_1[k:(k+100)], 0.975, axis=0),2)[0])

print(np.round(np.quantile(coverage1_1[k:(k+100)], 0.5, axis=0),2)[0])
print(np.round(np.quantile(coverage1_1[k:(k+100)], 0.025, axis=0),2)[0])
print(np.round(np.quantile(coverage1_1[k:(k+100)], 0.975, axis=0),2)[0])

print(np.round(np.quantile(l2_dist1_2[k:(k+100)], 0.5, axis=0),2)[0])
print(np.round(np.quantile(l2_dist1_2[k:(k+100)], 0.025, axis=0),2)[0])
print(np.round(np.quantile(l2_dist1_2[k:(k+100)], 0.975, axis=0),2)[0])


print(np.round(np.quantile(coverage1_2[k:(k+100)], 0.5, axis=0),2)[0])
print(np.round(np.quantile(coverage1_2[k:(k+100)], 0.025, axis=0),2)[0])
print(np.round(np.quantile(coverage1_2[k:(k+100)], 0.975, axis=0),2)[0])

print(np.round(np.quantile(timerun1_2[k:(k+100)], 0.5, axis=0),2)[0])
print(np.round(np.quantile(timerun1_2[k:(k+100)], 0.025, axis=0),2)[0])
print(np.round(np.quantile(timerun1_2[k:(k+100)], 0.975, axis=0),2)[0])

