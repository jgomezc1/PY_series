# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 20:25:46 2016

@author: Juan Carlos Vergara
"""
import numpy as np
import matplotlib.pyplot as pl


def Sa(Tmax, Nper, Aa, Av, Perfil, Uso):


#    Tmax=   10.0    #Máximo periodo estructural
#    Nper= 1000      #Número de periodos
    Dper= Tmax/Nper #Delta de periodo estructural
#    
#    Aa= 0.15        #Aceleración pico efectiva
#    Av= 0.20        #Velocidad pico efectiva
#    
#    Perfil= 4       #Perfil DE SUELO, A:1, B:2, C:3, D:4, E:5
#    Uso= 1          #Grupo de Uso
    
    
    
    
    
    #Fa   COEFICIENTE PARA AMPLifICAR EL ESPECTRO EN ROCA PARA TENER EN CUENTA LOS EFECTOS DE SITIO EN EL RANGO DE PERIODOS CORTOS
    #Fv   COEFICIENTE PARA AMPLifICAR EL ESPECTRO EN ROCA PARA TENER EN CUENTA LOS EFECTOS DE SITIO EN EL RANGO DE PERIODOS INTERMEDIOS
    #Impor    COEFICIENTE DE ImporTANCIA
    #T0       PERIODO DE VIBRACION LIMITE INFERIOR
    #TC       PERIODO DE VIBRACION INTERMEDIO
    #TL       PERIODO DE VIBRACION LIMITE SUPERIOR
    
    if (Perfil==1):
        Fa= 0.80
        Fv= 0.80
    
    if (Perfil==2):
        Fa=1.0
        Fv=1.0
    
    if (Perfil==3):
    
        if (Aa <= 0.20):
            Fa=1.20
    
        if (Aa > 0.20 and Aa <= 0.40):
            Fa=1.40-Aa
    
        if (Aa > 0.40):
            Fa=1.00
    
        if (Av <= 0.10):
            Fv=1.70
    
        if (Av > 0.10 and Av<= 0.50):
            Fv=1.80-Av
    
        if (Av > 0.50):
            Fv=1.30
    
    if (Perfil==4):
    
        if (Aa <= 0.10):
            Fa=1.60
    
        if (Aa > 0.10 and Aa <= 0.30):
            Fa=1.80-2.0*Aa
    
        if (Aa > 0.30 and Aa <= 0.50 ):
            Fa=1.50-Aa
    
        if (Aa > 0.50):
    
            Fa=1.00
    
        if (Av <= 0.10):
    
            Fv=2.40
    
    
        if (Av > 0.10 and Av<= 0.20):
    
            Fv=2.80-4.0*Av
    
        if (Av > 0.20 and Av <= 0.40):
    
            Fv=2.40-2.0*Av
    
    
    
        if ( Av > 0.40 and Av <= 0.50):
    
            Fv=2.00-Av
    
    
    
        if (Av > 0.50):
    
            Fv=1.50
    
    if (Perfil==5):
    
        if (Aa <= 0.10):
    
            Fa=2.0
    
        if (Aa > 0.10 and Aa <= 0.20):
    
            Fa=3.30-8.0*Aa
    
        if (Aa > 0.20 and Aa <= 0.30):
    
            Fa=2.70-5.0*Aa
    
        if (Aa > 0.30 and Aa <= 0.40):
    
            Fa=2.10-3.0*Aa
    
        if (Av <= 0.10):
    
            Fv=3.50
    
        if (Av > 0.10 and Av <= 0.20):
    
            Fv=3.80-3.0*Av
    
        if (Av > 0.20 and Av <= 0.40):
    
            Fv=4.0-4.0*Av
    
        if (Av > 0.40):
    
            Fv=2.40
    
    
    if (Uso==1):
        Impor=1.00
    
    if (Uso==2):
        Impor=1.10
    
    
    if (Uso==3):
        Impor=1.25
    
    
    if (Uso==4):
        Impor=1.50
    
    
    
    T0=0.10*Av*Fv/(Aa*Fa)
    
    TC=0.48*Av*Fv/(Aa*Fa)
    
    TL=2.40*Fv
    
    Dper=Tmax/Nper
    
    Sd=np.zeros([Nper,2])
    Sv=np.zeros([Nper,2])
    Sa=np.zeros([Nper,2])
    
    
    for i in range(Nper):
        T=Dper*(i+1)
    
        if (T <= T0):

            Sd[i,1]= 0.62*Aa*Fa*Impor*T**2*(0.40+0.60*T/T0)
            Sv[i,1]= 3.90*Aa*Fa*Impor*T*(0.40+0.60*T/T0)
            Sa[i,1]= 2.50*Aa*Fa*Impor*(0.40+0.60*T/T0)


        if (T > T0 and T <= TC):
    
            Sd[i,1]= 0.62*Aa*Fa*Impor*T**2
            Sv[i,1]= 3.90*Aa*Fa*Impor*T
            Sa[i,1]= 2.50*Aa*Fa*Impor
    
        if (T > TC and T <= TL):
    
            Sd[i,1]= 0.30*Av*Fv*Impor*T
            Sv[i,1]= 1.87*Av*Fv*Impor
            Sa[i,1]= 1.20*Av*Fv*Impor/T
    
        if (T > TL):
    
            Sd[i,1]= 0.30*Av*Fv*Impor*TL
            Sv[i,1]= 1.87*Av*Fv*Impor*TL/T
            Sa[i,1]= 1.20*Av*Fv*TL*Impor/T**2
    
    
        Sd[i,0]= T
        Sv[i,0]= T
        Sa[i,0]= T

    return Sa
    #pl.plot(Sa[:,0], Sa[:,1], color='black')




Tmax=   5.0    #Máximo periodo estructural
Nper= 2000      #Número de periodos
Dper= Tmax/Nper #Delta de periodo estructural

Aa= 0.15        #Aceleración pico efectiva
Av= 0.20        #Velocidad pico efectiva

Perfil= 4       #Perfil DE SUELO, A:1, B:2, C:3, D:4, E:5
Uso= 1          #Grupo de Uso


Sa=Sa(Tmax, Nper, Aa, Av, Perfil, Uso)


####################
valsy=np.linspace(0., .6, 7, endpoint=True)
valsx=np.linspace(0., 4., 9, endpoint=True)

pl.figure(figsize=(15./2.54, 4.5/2.54), facecolor='white')

pl.plot(Sa[:,0],Sa[:,1], color='black', linewidth=1.5)
pl.legend(loc='lower right', fontsize=7)
pl.xticks(valsx, fontsize=10)
pl.yticks(valsy, fontsize=10)
pl.ylabel(r'$\mathbf{S_a \ \left( \% g \right)}$',fontsize=10)
pl.xlabel(u'Periodo (s)',fontsize=10)
# pl.title(u'Promedio Relaciones Espectrales',fontsize=12)
pl.xlim(0.01, 4.)
pl.ylim(0., .6)
pl.grid()

pl.savefig("TipoB_NSR-10.png", dpi=600, facecolor='w', edgecolor='w',
           figsize=(30./2.54, 18./2.54), orientation='landscape',
           bbox_inches='tight')
#################
