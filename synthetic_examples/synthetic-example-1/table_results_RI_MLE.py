#!/usr/bin/env python
# coding: utf-8
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

l2_dist1_1=np.zeros((300,1))
l2_dist1_2=np.zeros((300,1))
width1_1=np.zeros((300,1))
width1_2=np.zeros((300,1))
coverage1_1=np.zeros((300,1))
coverage1_2=np.zeros((300,1))
timerun1_2=np.zeros((300,1))

for jjj in range(300):
    print(jjj)
    data = np.load('/output/simulation/mymethodfinal3/ESSsyn1_truthmle'+str(jjj+1)+'.npz')
    c=math.ceil((jjj+1)/100)
    g_mk_list3=data['aaa']
    g_mk_list2=data['aa']
    g_mk_list=data['a']
    points_inhomo=data['c']
    xxx=data['d']
    x=data['b']
    intensity=data['e']
    theta0=data['g']
    theta1=data['h']
    noise_var=data['i']
    integral_truth=data['m']
    integral_mean=data['n']
    integral_sd=data['o']
    timerun2=data['q']

    low=np.quantile(g_mk_list, 0.025, axis=0)
    high=np.quantile(g_mk_list, 0.975, axis=0)
    med=np.quantile(g_mk_list, 0.5, axis=0)
    truth=2*np.exp(-np.array(points_inhomo)/15)+np.exp(-((np.array(points_inhomo)-25)/10)**2)
    truth=c*truth
    l2_dist1_1[jjj]=sum((np.array(med).squeeze()-truth)**2)
    coverage1_1[jjj]=np.sum((truth>=low.squeeze()) * (truth<=high.squeeze()))/len(points_inhomo)
    width1_1[jjj]=sum(high-low)/len(points_inhomo)

    low=np.quantile(g_mk_list2, 0.025, axis=0)
    high=np.quantile(g_mk_list2, 0.975, axis=0)
    med=np.quantile(g_mk_list2, 0.5, axis=0)
    truth=2*np.exp(-np.array(xxx)/15)+np.exp(-((np.array(xxx)-25)/10)**2)
    truth=c*truth
    l2_dist1_2[jjj]=sum((np.array(med).squeeze()-truth)**2)
    coverage1_2[jjj]=np.sum((truth>=low.squeeze()) * (truth<=high.squeeze()))/len(xxx)
    width1_2[jjj]=sum(high-low)/len(xxx)
    timerun1_2[jjj]=timerun2
    
for k in ([0,100,200]):
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
    
    print( [np.round(np.mean(timerun1_2[k:(k+100)]/5),2),  np.round(np.std(timerun1_2[k:(k+100)]/5) ,2) ])
