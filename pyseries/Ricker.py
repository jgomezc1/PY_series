# -*- coding: utf-8 -*-
"""
Created on Tuesday Nov  1 11:25:26 2016

@author: jvergar2
@email: jvergar2@gmail.com
"""

import numpy as np
import RDCsignals as rdc
import RDClectura as lec

#
Tt=18.0
Tc=1.
fc=2.5
Nt=18001
dt = Tt/(Nt-1)

Rick, T=rdc.ricker(Nt, Tt, Tc, fc)
np.savetxt("pulso.txt", Rick)
signal = lec.readsignal_VEC('pulso' , 1.0)
ndats = len(signal)
lec.grafsignalA(signal, 'pulso' , -1.0 , 1.0 , dt , 1)
x, Sas , nfs = rdc.Ftrans(signal , ndats , dt , 10.0)
rdc.grafFourier(Sas , x , nfs , 'spectro' , 0.0 , 10.0 , 0.0 , 200.0 , 2)