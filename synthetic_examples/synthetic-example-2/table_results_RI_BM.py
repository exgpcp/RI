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
timerun2=np.zeros((300,1))
width2_1=np.zeros((300,1))
width2_2=np.zeros((300,1))
for jjj in range(300):
    print(jjj)
    data = np.load('/output/simulation/mymethod_BM1_5_4/ESSsyn2_BM1_'+str(jjj+301)+'.npz')
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
      




