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
sys.path.append('/Users/casierraa/dropbox/REPOS/PY_series/pyseries/')
import lectura as lec
import signals as sig

lec.plot_txt_waves('vhill' , 0.0005)