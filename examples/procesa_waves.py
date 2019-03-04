#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 10:39:23 2018

@author: casierraa
"""
from __future__ import division
from os import sys
sys.path.append('../pyseries/')
import lectura as lec
import signals as sig
import numpy as np
import matplotlib.pyplot as plt

y , U , V = lec.plot_txt_waves('Dcolina' , 0.0005)

N = len(U)
T_t = N*0.0005
#
FS = 10.0
x1 , Samag1 , A1 , nfs1 = sig.Ftrans(U , N , T_t/N , FS)
lec.grafFourier(Samag1 , x1 , nfs1 , 'result' , 0.0 , 16.0 , 0.0 , 1000.0 , 2)

#N = 3999
#T_t = 2.0
#tc  = 1.0
#fc  = 4.0
#t_b = np.sqrt(6)/np.pi/fc
#print ("Breadth =") , t_b , ('s')
#Rick , time = sig.ricker(N , T_t , tc, fc)
#lec.grafsignalG(Rick , 'Pulse', 'Disp','L' ,-1.0,+1.0, T_t/N , 3)
#
#x2 , Samag2 , A2 , nfs2 = sig.Ftrans(Rick , len(Rick) , T_t/N , FS)
#lec.grafFourier(Samag2 , x2 , nfs2 , 'result' , 0.0 , 16.0 , 0.0 , 750.0 , 4)
#
#plt.figure(5)
#plt.plot(x2 , Samag2/Samag1)
