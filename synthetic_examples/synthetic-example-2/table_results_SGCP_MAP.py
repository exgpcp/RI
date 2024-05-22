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

l2_dist2_1=np.zeros((100,1))
l2_dist2_2=np.zeros((100,1))
coverage2_1=np.zeros((100,1))
coverage2_2=np.zeros((100,1))
timerun2_1=np.zeros((100,1))
timerun2_2=np.zeros((100,1))
width2_1=np.zeros((100,1))
width2_2=np.zeros((100,1))

for jjj in range(100):
    print(jjj)
    data = np.load('/output/simulation/adam1_1/priormle/syn2/syn2_priormle_'+str(jjj+301)+'.npz')    
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
    timerun2=data['p']

    l2_dist2_1[jjj]=l2_dist1
    coverage2_1[jjj]=coverage1
    width2_1[jjj]=width1
    l2_dist2_2[jjj]=l2_dist2
    coverage2_2[jjj]=coverage2
    width2_2[jjj]=width2
    timerun2_2[jjj]=timerun2

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
