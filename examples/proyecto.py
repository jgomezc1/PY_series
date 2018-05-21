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
from sympy import init_printing
from os import sys
sys.path.append('/Users/juan/Dropbox/REPOS/PY_series/pyseries/')
import lectura as lec
import signals as sig

y = lec.loadcsv('manyanet' , 'armenia')
lec.grafsignalA(y , 'armenia' , -12.0 , 12.0 , 0.005 , 0)
nt = len(y)
T_max = 5.0
Sa , Ta = sig.respesp(y , nt , 0.005 , T_max, 50)
plt.figure(1)
plt.grid()
plt.xlabel('Period (s)')
plt.ylabel('Sa (gales)')
plt.xlim(0.001 , T_max )
plt.plot(Ta , Sa)