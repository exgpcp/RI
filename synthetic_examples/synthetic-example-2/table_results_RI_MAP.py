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

l2_dist2_1=np.zeros((300,1))
l2_dist2_2=np.zeros((300,1))
coverage2_1=np.zeros((300,1))
coverage2_2=np.zeros((300,1))
timerun2_2=np.zeros((300,1))
width2_1=np.zeros((300,1))
width2_2=np.zeros((300,1))

for jjj in range(300):
    print(jjj)
    data = np.load('/output/simulation/mymethodfinal3/ESSsyn2_priormle'+str(jjj+301)+'.npz')
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
    timerun2=data['q']

    low=np.quantile(g_mk_list, 0.025, axis=0)
    high=np.quantile(g_mk_list, 0.975, axis=0)
    med=np.quantile(g_mk_list, 0.5, axis=0)
    truth=10+points_inhomo-points_inhomo
    truth=c*truth
    l2_dist2_1[jjj]=sum((np.array(med).squeeze()-truth)**2)
    width2_1[jjj]=sum(high-low)/len(points_inhomo)
    coverage2_1[jjj]=np.sum((truth>=low.squeeze()) * (truth<=high.squeeze()))/len(points_inhomo)
    low=np.quantile(g_mk_list2, 0.025, axis=0)
    high=np.quantile(g_mk_list2, 0.975, axis=0)
    med=np.quantile(g_mk_list2, 0.5, axis=0)
    truth=10+xxx-xxx
    truth=c*truth
    l2_dist2_2[jjj]=sum((np.array(med).squeeze()-truth)**2)
    width2_2[jjj]=sum(high-low)/len(xxx)
    coverage2_2[jjj]=np.sum((truth>=low.squeeze()) * (truth<=high.squeeze()))/len(xxx)
    timerun2_2[jjj]=timerun2

for k in ([0,100,200]):
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
           
    print( [np.round(np.mean(timerun2_2[k:(k+100)]/5),2),  np.round(np.std(timerun2_2[k:(k+100)]/5) ,2) ])
