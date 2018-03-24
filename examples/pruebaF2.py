#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 10:39:23 2018

@author: casierraa
"""

from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from sympy import *
import fourier as fou
from sympy import init_printing
init_printing()
#
# Playing with the Fourier transform
#
def harmonic(n , x , sn): 
    an = 8.0/(np.pi*(2*n-1))
    sn[:] = an*np.sin((2.0*n-1.0)*x[:])
    plt.plot(x , sn)
    
    return sn


x = np.arange(-5.0 , 5.0 , 0.05)
n = len(x)
sum = np.zeros(n)
sn  = np.zeros(n)
# S0
sum[:] = 2.0
#
nn = 20
an  = np.zeros(nn)
fn  = np.zeros(nn)
for n in range(1 ,nn):
    sn = harmonic(n , x , sn)
    sum[:] = sum[:] + sn[:]
    sn[:] = 0.0
plt.figure(nn+1)
plt.plot(x , sum)
plt.show()

plt.figure(nn+2)
plt.plot(fn , an)
plt.show()
#
# PLaying with the Fourier transform
#
N   = 2048
T_t = 4.0
x = np.arange(-2.0 , 2.0 ,T_t/N)
n = len(x)
FX = np.zeros((n), dtype = complex)
xm = x 
pot = -np.pi*xm**2
FX = np.cos(6.0*np.pi*xm)*np.exp(pot)
fou.grafsignalG(FX , T_t/N , 0)

FS = 10.0
x , Samag , A , nfs = fou.Ftrans(FX , N , T_t/N , FS)
fou.grafFourier(Samag , x , nfs , 1)
#
# Now let us proceed numerically
#
ZZ = np.fft.fft(FX)
XX = np.fft.ifft(ZZ)
#
# Let us plot the FAS
#
plt.figure(2)
plt.plot(abs(ZZ))
plt.show()
#
# And now again the original function
#
plt.figure(3)
plt.plot(xm , XX)
plt.show()
###
