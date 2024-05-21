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

l2_dist1_1=np.zeros((200,1))
l2_dist1_2=np.zeros((200,1))
width1_1=np.zeros((200,1))
width1_2=np.zeros((200,1))
coverage1_1=np.zeros((200,1))
coverage1_2=np.zeros((200,1))
timerun1_2=np.zeros((200,1))

for jjj in range(200):
    print(jjj)
    data = np.load('/output/simulation/adam1_1/BM/syn1/syn1_BM_'+str(jjj+1)+'.npz')
    c=math.ceil((jjj+1)/100)
    g_mk_list2=data['aa']
    g_mk_list=data['a']
    points_inhomo=data['c']
    xxx=data['d']
    coverage1=data['i']
    tau_list=data['e']
    measure_sup=data['g']
    noise_var=data['h']
    coverage2=data['j']
    l2_dist1=data['k']
    l2_dist2=data['l']
    width1=data['m']
    width2=data['n']
    timerun2=data['p']

    l2_dist1_1[jjj]=l2_dist1
    coverage1_1[jjj]=coverage1
    width1_1[jjj]=width1

    l2_dist1_2[jjj]=l2_dist2
    coverage1_2[jjj]=coverage2
    width1_2[jjj]=width2
    timerun1_2[jjj]=timerun2
   
for k in ([0,100]):
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
