"""

"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter
from scipy.signal import freqz
import lectura as lec
#
def linebaseace(a):
    """
     Conducts baseline correction of the acceleration time signal a[].
     """
    FacTaper = 5.
    ndats = len(a)
    ndatst = int(FacTaper / 100 * ndats)
    Tapper = np.zeros([ndats], dtype=float)

    ac = a - np.mean(a)
#
    for pos in range(ndats):
        if pos < ndatst:
            Tapper[pos] = 0.5 * (1. - np.cos(np.pi * pos / ndatst))
        elif pos < (ndats - ndatst):
            Tapper[pos] = 1.0
        else:
            Tapper[pos] = 0.5 * (1. - np.cos(np.pi * (ndats - pos) / ndatst))

    actap = ac * Tapper
#
    return(actap)
#
def linebasevel(ao,v,time):

    # ****Statement of variables *****
    #   v               : initial acceleration
    #   ao              : value mean for adjust line base
    #   vc              : adjusted acceleration
    # ***** End *****

    vlb = ao * time

    vc = v - vlb
    return(vc)
#
def Ftrans(datos, ndats, dt, fs):
    """
    Compute the Fourier spectra of datos[] of
    length ndats and sampled at dt.
    Returns the result in Samag after smoothing by the
    smoothing factor fs.
    """
    nfr = int(ndats/2)
    df = 1.0/(ndats*dt)
    x = np.arange(df,nfr*df, df)
    A = np.fft.fft(datos)
    Aa = np.abs(A)

    # Smooth the spectrum.
    Sa = Aa[1:nfr]
    Samag = smooth(Sa , x , fs) 
    nfs = nfr-1
    return x , Samag , A , nfs

#
def FtransV(datos , ndats , dt):
#
# Intgrates an acceleration history into a velocity history
# proceeding in the frequency domain.
#
	x=np.fft.fftfreq(ndats, dt)
	A=np.fft.fft(datos)
	Aomega=A
	
	for i in range(ndats-1):
		Aomega[i+1] = A[i+1]/(1j*x[i+1]*2.0*np.pi)

	B=np.fft.ifft(Aomega)
	Bb = np.real(B)
#
	return Bb    
#
#
def IFtrans(datos , ndats , dt):
#
	B=np.fft.ifft(datos)
#
	return np.real(B)
#
def smooth(Sa, Freq , fftfs):
#
    #  **** Input :
    #  Sa: Original spectrum
    #  Freq: Frequency
    #  fftfs: Smoothing factor
#
    Sas  = np.zeros([len(Sa)],dtype=float)
    fia = 1
    fma = 1
    suma = Sa[0] * Sa[0]
    pot = 1./2./fftfs
    fsexpi = 2**(-pot)
    fsexpm = 2**( pot)
    Sas[0] = Sa[0]
    NNfft = len(Sa)
    for i in range(1,NNfft):
    # #    for i=2:NNfft
        fi = int((i + 1) * fsexpi)
        fm = int((i + 1) * fsexpm)
        if fi < 1:
            fi = 1
        if fm > NNfft:
            fm = NNfft

        for Num in range(fia - 1, fi - 1):
        #         #for j=fia:fi-1:
            suma = suma - Sa[Num] * Sa[Num]

        for Num in range(fma, fm):
        #         #for j=fma+1:fm:
            suma = suma + Sa[Num]*Sa[Num]

        Nf = fm - fi + 1
        fia = fi
        fma = fm
        Sas[i]=np.sqrt(suma/Nf)
#
    return (Sas)
#
def cutsignal(senal ,  t0 , tf , dt):
    ndats = int((tf-t0)/dt)
    nsenal  = np.zeros(ndats,dtype=float)
    NI = int(t0/dt)
    for i in range(ndats):
        nsenal[i] = senal[NI+i]
#    
    return nsenal   
#
def rqma(senal):
    ndats  = len(senal)
    asq = 0.0
    for i in range(ndats):
        asq = asq + senal[i]*senal[i]
    qma = asq / ndats
    rms = np.sqrt(qma)
#
    return rms    
#
def filtro(senal , lowcut , highcut , fs , filtype , order , NGra):
    """
     senal  : time history.
     lowcut : lowcut frequency
     highcut: Highcut frequency
     fs     : sampling (200 samples per second in the CUSP)
     filtype: highpass, lowpass, bandpass.
     order  : integer. Speed of the decreasing branch.
    """
    b, a = butter_bandpass(lowcut, highcut, fs, filtype  , order=order)
    w, h = freqz(b, a , worN=2000)
    plt.figure(NGra)
    plt.plot((fs * 0.5 / np.pi) * w, abs(h), label="order = %d" % 4)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Gain')
    y = butter_bandpass_filter(senal , lowcut, highcut, filtype , fs, order=6)

    return y    
#
def butter_bandpass(lowcut, highcut, fs, filtype  , order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    if (filtype == 'highpass'):
        b, a = butter(order, low , btype = filtype)
    else:
        b, a = butter(order, [low , high] , btype=filtype)
#    
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, filtype , fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, filtype , order=order)
    y = lfilter(b, a, data)
    return y
#
def qfspectra(datos1,datos2,ndats):   
#
    TF=np.zeros([ndats],dtype=float) 
    TF=datos1/datos2
#
    return TF
#
def ricker(nt, Tt, tc, fc):
	
    Rick=np.zeros(nt)
    T=np.zeros(nt)
    dt=Tt/(nt-1)
	
    for i in range(nt):
        tao=np.pi*fc*(dt*i-tc)
        Rick[i]=(2.*tao**2-1.)*np.exp(-tao**2)
        T[i]= i*dt

    return Rick, T
#
def Sa(Tmax, Nper, Aa, Av, Perfil, Uso):


#    Tmax=   10.0    #Maximo periodo estructural
#    Nper= 1000      #Numero de periodos
    Dper= Tmax/Nper #Delta de periodo estructural
#    
#    Aa= 0.15        #Aceleracion pico efectiva
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


#Rutina para calcular el espectro de repuesta a un acelerograma

###
# acel1: Historia de aceleraciones
# nt: numero de tiempos
# dt: delta de tiempo
# tmax: periodo maximo hasta el cual se va a evaluar
# nper: numero de periodos a evaluar
# Todos los espectros se van a calcular con un amortiguamiento del 5 porciento, pero
# es posible cambiar dicho valor dentro de la rutina en la variable si
###
def respesp(acel1, nt, dt, tmax, nper):
    """
     acel1: Historia de aceleraciones
     nt: numero de tiempos
     dt: delta de tiempo
     tmax: periodo maximo hasta el cual se va a evaluar
     nper: numero de periodos a evaluar
     Todos los espectros se van a calcular con un amortiguamiento del 5 porciento, pero
     es posible cambiar dicho valor dentro de la rutina en la variable si
    """	
    dper=tmax/nper
    nt=nt+1
    acel=np.zeros(nt)
    Ta=np.zeros(nper)
    acel[:nt-1]=acel1
    del acel1
    si=0.02		#Amortiguamiento

    dper=tmax/nper
    for i in range(nper):
        Ta[i]=dper*i
    sa=np.zeros(nper)
    for i in range(nper):
        w=2.*np.pi/(i+1)/dper
        wd=w*np.sqrt(1.-si**2)
	
        xt0=0.
        vt0=0.
        sd0=0.
        for j in range(nt-1):
            a0=acel[j]
            pa=(acel[j+1]-acel[j])/dt
            d=-pa/w**2
            cc=(-a0+2.*si*pa/w)/w**2
            g1=xt0-cc
            g2=(vt0+si*w*g1+pa/w**2)/wd
            xt0=np.exp(-si*w*dt)*(g1*np.cos(wd*dt)+g2*np.sin(wd*dt))+cc+d*dt
		
            vt0=np.exp(-si*w*dt)*(-g1*wd*np.sin(wd*dt)+g2*wd*np.cos(wd*dt))
            vt0=vt0-si*w*np.exp(-si*w*dt)*(g1*np.cos(wd*dt)+g2*np.sin(wd*dt))
            vt0=vt0+d

            if (sd0<np.abs(xt0)):
                sd0=np.abs(xt0)
			
        sa[i]=sd0*w**2
    return sa , Ta

#def respesp(acel1, nt, dt, tmax, nper):
#    """Computes the acceleration response spectra for an
#    acceleration time history using excitation interpolation
#    based method (See Chopra Section 5.2).
#
#    Parameters
#    ----------
#    acel1 : ndarray (float)
#        Array with acceleration time history:
#    nt : scalar (float)
#        Number of data points.
#    dt : scalar (float)
#        Time step of the signal.
#    tmax : scalar (float)
#        Maximum period.
#    nper : scalar (float)
#       Number of periods to evaluate.
#    Todos los espectros se van a calcular con un amortiguamiento del 5 porciento, pero
#    es posible cambiar dicho valor dentro de la rutina en la variable si
#    """	
#    dper=tmax/nper
#    nt=nt+1
#    acel=np.zeros(nt)
#    Ta=np.zeros(nper)
#    acel[:nt-1]=acel1
#    del acel1
## Change damping if desired
#    si=0.05		
#
#    dper=tmax/nper
#    for i in range(nper):
#        Ta[i]=dper*i
#    sa=np.zeros(nper)
#    for i in range(nper):
#        w=2.*np.pi/(i+1)/dper
#		wd=w*np.sqrt(1.-si**2)
#	
#		xt0=0.
#		vt0=0.
#		sd0=0.
#	
#		for j in range(nt-1):
#			a0=acel[j]
#			pa=(acel[j+1]-acel[j])/dt
#			d=-pa/w**2
#			cc=(-a0+2.*si*pa/w)/w**2
#			g1=xt0-cc
#			g2=(vt0+si*w*g1+pa/w**2)/wd
#			xt0=np.exp(-si*w*dt)*(g1*np.cos(wd*dt)+g2*np.sin(wd*dt))+cc+d*dt
#		
#			vt0=np.exp(-si*w*dt)*(-g1*wd*np.sin(wd*dt)+g2*wd*np.cos(wd*dt))
#			vt0=vt0-si*w*np.exp(-si*w*dt)*(g1*np.cos(wd*dt)+g2*np.sin(wd*dt))
#			vt0=vt0+d
#
#			if (sd0<np.abs(xt0)):
#				sd0=np.abs(xt0)
#			
#		sa[i]=sd0*w**2
#    return sa , Ta

#
def transfunction(datos1,datos2,ndats,var1,dt,Nro):
    nfr=ndats/2
    df=1.0/(ndats*dt)
    x=np.zeros([nfr-1], dtype=float)
    x=np.arange(df,nfr*df,df)
    A1   = np.zeros([ndats],dtype=float)
    Aa1  = np.zeros([ndats],dtype=float) 
    A1=np.fft.fft(datos1)
    Aa1=np.abs(A1)
#
    A2   = np.zeros([ndats],dtype=float)
    Aa2  = np.zeros([ndats],dtype=float) 
    A2=np.fft.fft(datos2)
    Aa2=np.abs(A2)   
#
    TF=np.zeros([ndats],dtype=float) 
    TF=Aa1/Aa2
#
    plt.figure(Nro)
    plt.plot(x,TF[nfr+1:])
#    plt.plot(TF[1:ndats/2])
    plt.grid()  
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude')
    plt.legend(['Transfer function'])
    plt.xlim(0.0,10)
    plt.ylim(0.0,30)
#    plt.xscale('log')
#    plt.yscale('log')
    plt.savefig(var1)
#
def cut_tail(signal, nnew):
    """
     Given an initial signal generates a new one
     of nnew samples.
    """
    signalnew = np.zeros([nnew], dtype=float)
    for i in range(0 , nnew):
        signalnew[i] = signal[i]
    return signalnew
#
def splitsignal(signal,ndats,nr):
    """
     Given a single signal of ndats samples
     generate nr signals of nt = ndats/nr samples.
     
     signal : signal ndats
     nr     : number of windows
     nt     : 
     mati   : array with signals
    """
    nt =int(ndats/nr)
    mati = np.zeros([nt , nr], dtype=float)
    li = 0
    for i in range(0,nr):
        ls = li + (nt-1)
        j = 0
        for k in range(li , ls):
            mati[j , i] = signal[k]
            j = j + 1
        li = ls + 1#
    return nt , mati 
#
def hvlist(signals1 , signals2 , signals3 , nt , nvr , Rlist):
    """
     From a list of signals select the correct ones
     and compute H to V ratio.
     signals1, signals2 : horizonal signals 
     signals3           : vertical signals
    """
    xmin = 1.0; xmax = 100; ymin=0.1; ymax = 10
    for i in range(0 , nvr):
        idr = Rlist[i]
        x, Samag1 , A1 , nfs = Ftrans(signals1[: , idr] , nt , 0.01 , 10.0)
        x, Samag2 , A2 , nfs = Ftrans(signals2[: , idr] , nt , 0.01 , 10.0)
        u1s =  Samag1**2
        u2s =  Samag2**2
        SamagC = np.sqrt(u1s +u2s)
        x, Samag3 , A3 , nfs = Ftrans(signals3[: , idr] , nt , 0.01 , 10.0)
        FT = SamagC/Samag3 
        lec.grafFourier(FT , x , nfs , 'FAS', xmin , xmax , ymin , ymax , 20)
    return
#
#
def SaT(Per , Nper, Aa, Av, Perfil, Uso):


#    Tmax=   10.0    #Maximo periodo estructural
#    Nper= 1000      #Numero de periodos
#    Dper= Tmax/Nper #Delta de periodo estructural
#    
#    Aa= 0.15        #Aceleracion pico efectiva
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
    
#    Dper=Tmax/Nper
    
    Sd=np.zeros([Nper,2])
    Sv=np.zeros([Nper,2])
    Sa=np.zeros([Nper,2])
    
    
    for i in range(Nper):
#        T=Dper*(i+1)
        T = Per[i]
    
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




































