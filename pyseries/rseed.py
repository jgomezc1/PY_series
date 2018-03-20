# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 14:04:25 2016

@author: casierraa
"""

from obspy import read
import numpy as np


path = '../DATOS/'
var = 'SS.S175..EHE.D.2016.333'

st = read(var)
print st

tr = st[0]

print tr

print tr.stats

#tr.stats.station
#
#tr.stats.gse2.datatype
#
#tr.data

np.savetxt('Este.txt', tr)

st.plot()
























