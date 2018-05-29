"""
@author: Risk and Design Consulting
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import signals as sig
#
def readsignal_MAT(var, fa):
    """
     Reads the time signal stored in the file var.dat and
     written in 8-column format. Returns the signal into the
     single vector signal.
     fa is an instrumental amplification factor
    """
    path = '../data/'
    channel = np.loadtxt(path + var + '.dat')
    nrows  = len(channel[:])
    ndats  = nrows*8
    signal   = np.zeros([ndats], dtype=float)
    k=0
    for i in range(0,nrows):
        for j in range(0,8):
            signal[k]=channel[i,j]*fa
            k=k+1
#
    return ndats , signal
#
def readsignal_ncol(var , ncol , fa):
    """
     Reads the time signal stored in the file var (variable extension)
     and written in ncol-column format. Returns the signal into the
     single vector signal.
     fa is a units convesrion factor
    """
    path = '../data/'
    channel = np.loadtxt(path + var)
    nrows  = len(channel[:])
    ndats  = nrows*ncol
    signal   = np.zeros([ndats], dtype=float)
    k=0
    for i in range(0,nrows):
        for j in range(0,ncol):
            signal[k]=channel[i,j]*fa
            k=k+1
#
    return ndats , signal
#
def readsignal_SAF(var , fac):
    """
     Reads the time signal stored in the file var.saf and
     written in 3-column format corresponding to the
     vertical, north and east components. Returns
     each component in independent vectors like
     SIGNALV = COL1*fac where fac is a unit conversion factor.
    """
    path = '../data/'
    channel = np.loadtxt(path + var + '.saf' , skiprows = 14)
    ndats  = len(channel[:])
    signalV   = np.zeros([ndats], dtype=float)
    signalN   = np.zeros([ndats], dtype=float)
    signalE   = np.zeros([ndats], dtype=float)
    for i in range(0,ndats):
            signalV[i]=channel[i,0]*fac
            signalN[i]=channel[i,1]*fac
            signalE[i]=channel[i,2]*fac
#
    return ndats , signalV , signalN, signalE
#
def readsignal_VEC(name , fa):
    """
     Reads the time signal stored in the file var.txt and
     written in a single column format. Returns the signal into the
     single vector signal.
     fa is an instrumental amplification factor
    """
    path = '../data/'
    channel = np.loadtxt(path + name + '.txt')
    ndats  = len(channel)
    signal   = np.zeros([ndats], dtype=float)
    for i in range(ndats):
        signal[i]=channel[i]*fa
#
    return ndats , signal
#
def readsignal_OPENSEES(name , fa):
    """
     Reads the OPENSEES generated time signal stored in
     the file name.out and written in a 3 column format
     storing time step, u-field and v-field.
     Returns the signal into the single vector signal.
    """
    path = '../data/'
    channel = np.loadtxt(path + name + '.out')
    ndats  = len(channel)
    signal   = np.zeros([ndats , 3], dtype=float)
    for i in range(ndats):
        signal[i , :]=channel[i, :]*fa
#
    return ndats , signal


#
def grafsignalA(U, var , ymin , ymax , dt , Ngra):
    """
     Plots the acceleration history U[ndats] into
     Ngra. The plot is also stored into var.pdf
    """
    path = '../results/'
    ndats  = len(U)
    x=np.zeros([ndats], dtype=float)
    x=np.arange(0,ndats*dt,dt)
    var1=path + var + '.pdf'
    plt.figure(Ngra)
    plt.plot(x,U)
    plt.grid()
    plt.xlabel('Tiempo (s)')
    plt.ylabel('Aceleracion (gales)')
    plt.ylim(ymin,ymax)
    plt.savefig(var1)
#
    return
#
def grafsignalV(U, var , dt , Ngra):
    """
     Plots the velocity history U[ndats] into
     Ngra. The plot is also stored into var.pdf
    """
    path = '../results/'
    ndats  = len(U)
    x=np.zeros([ndats], dtype=float)
    x=np.arange(0,ndats*dt,dt)
    var1=path + var + '.pdf'
    plt.figure(Ngra)
    plt.plot(x,U)
    plt.grid()
    plt.xlabel('Tiempo (s)')
    plt.ylabel('Velocidad (cm/s)')
    plt.savefig(var1)
#
    return
#
def grafsignalG(A , var1 , label1 , units , ymin , ymax , dt , Ngra):
    """
     Plots the generalized time signal A[ndats] into Ngra
     The plot is also stored into var.pdf
    """
    path = '../results/'
    ndats  = len(A)
    x=np.zeros([ndats], dtype=float)
    x=np.arange(0,ndats*dt,dt)
    var1=path + var1 + '.pdf'
    plt.figure(Ngra)
    plt.plot(x,A)
    plt.grid()
    plt.xlabel('Tiempo (s)' + str(Ngra))
    plt.ylabel(label1 + ' ' + units)
#    plt.ylim(ymin,ymax)
    plt.savefig(var1)
#
    return
#
def grafFourier(Sas , x , nfr , var, xmin , xmax , ymin , ymax , Nfig):
    """
     Plots the Fourier spectral amplitude Sas into Nfig.
     Sas : Spectrum
     x   : frecuency
     xmin,xmax,ymin,ymax
    """
    path = '../results/'
    plt.figure(Nfig)
    plt.plot(x,Sas)
    var1= path + var + '.pdf'
    plt.grid()  
    plt.xlabel('Frecuencia (Hz)')
    plt.ylabel('Amplitud')
    #plt.legend(['Fourier spectral amplitude'])
#    plt.xlim(xmin,xmax); plt.ylim(ymin,ymax)
#    plt.xscale('log')
#    plt.yscale('log')
    plt.savefig(var1)
#
    return

def grafFourierDouble(Sas1 , Sas2, x1 , x2, nfr , var, xmin , xmax , ymin , ymax,Nfig):
    """
     Plots the Fourier spectral amplitude Sas into Nfig.
     Sas : Spectrum
     x   : frecuency
     xmin,xmax,ymin,ymax
    """
    path = '../results/'
    plt.figure(Nfig)
    plt.plot(x1,Sas1)
    plt.plot(x2,Sas2)
    var1= path + var + '.pdf'
    plt.grid()  
    plt.xlabel('Frecuencia (Hz)')
    plt.ylabel('Amplitud')
    #plt.legend(['Fourier spectral amplitude'])
    #plt.legend(['Sensor Sixaola', 'Sensor Goebox'])
    plt.xlim(xmin,xmax)
    #plt.ylim(ymin,ymax)
    plt.xscale('log')
    #plt.yscale('log')
    plt.savefig(var1)
#
    return    
# 
def grafFourierG(Sas , x , nfr , var,  Nfig):
    """
     Plots the Fourier spectral amplitude Sas into
     Nfig.
    """
    path = '../results/'
    plt.figure(Nfig)
    plt.plot(x,Sas)
    var1= path + var + '.pdf'
    plt.grid()  
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude')
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig(var1)
#
    return
#
def plot_csv(filename, xlabel, ylabel,var):
#
    path1 = '../data/' 
    path2 = '../results/'
    DATOS=np.genfromtxt(path1 + filename +'.csv',delimiter=',')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    X = DATOS[:,0]
    Y = DATOS[:,1]
    plt.plot(X,Y)
#
    var1= path2 + var + '.pdf'
    plt.grid()  
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
#    plt.xscale('log')
#    plt.yscale('log')
    plt.savefig(var1)
#
#    plt.xscale('log')
    plt.grid()
    plt.figure()
#
    return DATOS
#
def plot_csvf(filename):
#
    path = '../data/'   
    DATOS=np.genfromtxt(path + filename +'.csv',delimiter=',')
    X = DATOS[:,0]
    Y = DATOS[:,1]
    plt.plot(X,Y)
#    plt.xscale('log')
#    plt.yscale('log')
    plt.grid()
    plt.figure()
#
    return  X , Y  
#
def plot_txt(filename, xlabel, ylabel):
# 
    path = '../data/'    
    DATOS=np.loadtxt(path + filename +'.txt' )
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    X = DATOS[:,0]
    Y = DATOS[:,1]
    plt.plot(X,Y)
    plt.grid()
    plt.figure()
    
    return
#
def plot_txtf(filename):
#
    path = '../data/'     
    DATOS=np.loadtxt(path + filename +'.txt' )
    X = DATOS[:,0]
    Y = DATOS[:,1]
    plt.plot(X,Y)
    plt.grid()
    plt.figure()
    
    return 
#
def plot_csvf2(filename1,filename2,label1,label2):
#
    path = '../data/'
    DATOS1=np.genfromtxt(path + filename1 +'.csv',delimiter=',')
    DATOS2=np.genfromtxt(path + filename2 +'.csv',delimiter=',')
    X = DATOS1[:,0]
    Y = DATOS1[:,1]
#
    x = DATOS2[:,0]
    y = DATOS2[:,1]    
#
    plt.plot(X,Y, 'k--')
    plt.plot(x,y, 'g')
    plt.legend( (label1, label2), loc=2, fontsize='small' )
#    plt.xscale('log')
    plt.grid()
    plt.figure()
#
    return
#
def plot_csvf3(filename1,filename2,filename3,label1,label2,label3):
#
    path = '../data/'
    DATOS1=np.genfromtxt(path + filename1 +'.csv',delimiter=',')
    DATOS2=np.genfromtxt(path + filename2 +'.csv',delimiter=',')
    DATOS3=np.genfromtxt(path + filename3 +'.csv',delimiter=',')
    X1 = DATOS1[:,0]
    Y1 = DATOS1[:,1]
#
    X2 = DATOS2[:,0]
    Y2 = DATOS2[:,1]    
#
    X3 = DATOS3[:,0]
    Y3 = DATOS3[:,1]    
#
    plt.plot(X1,Y1, 'k--')
    plt.plot(X2,Y2, 'g-.')
    plt.plot(X3,Y3, 'r')
    plt.legend( (label1, label2, label3), loc=2, fontsize='small' )
    plt.xscale('log')
    plt.grid()
    plt.figure()
#
    return    
#
def plot_txtf_log(filename):
#
    path = '../data/'    
    DATOS=np.loadtxt(path + filename +'.txt' )
    X = DATOS[:,0]
    Y = DATOS[:,1]
    plt.plot(X,Y)
    plt.xscale('log')
    plt.grid()
    plt.figure()
    
    return
#
def read_csvf(filename , icol , xlabel , ylabel , var , title):
    """
     Reads and plot an ncol csv file contained in
     filename and plots the column corresponding to ncol-icol
    """
#
    path1 = '../data/'
    path2 = '../results/' 
    DATOS=np.genfromtxt(path1 + filename +'.csv',delimiter=',')
    nrow , ncol = np.shape(DATOS)
    X = DATOS[:,0]
    Y = DATOS[:,ncol-icol]
    plt.grid()
    plt.plot(X,Y)  
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim(0.0, 5.0)
    var2=path2 + var + '.pdf'
    plt.savefig(var2)    
#
    return X , Y
#
def read_multiple_csvf(filename , xlabel , ylabel , var , title):
    """
     Reads and plot an ncol filename.csv file
    """
#  
    path1 = '../data/'
    path2 = '../results/'
    var2=path2 + var + '.pdf'
    DATOS=np.genfromtxt(path1+ filename +'.csv',delimiter=',')
    nrow , ncol = np.shape(DATOS)
    plt.figure()
    X = DATOS[:,0]
    for i in range(1 , ncol):
        Y = DATOS[:,i]
        plt.plot(X,Y)    
    plt.grid()
#    plt.xscale('log')
#    plt.yscale('log')
    plt.title(title)
    plt.xlim(0.0, 5.0)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(var2)
    return
#
def read_multiple_csvf_log(filename , xlabel , ylabel , title):
    """
     Reads and plot an ncol CSV file contained in
     filename and plots the column corresponding to ncol-icol
    """
#  
    path = '../data/'  
    DATOS=np.genfromtxt(path + filename +'.csv',delimiter=',')
    nrow , ncol = np.shape(DATOS)
    plt.figure()
    X = DATOS[:,0]
    for i in range(1 , ncol):
        Y = DATOS[:,i]
        plt.plot(X,Y)
    
    plt.grid()
    plt.xscale('log')
#    plt.yscale('log')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    return
#
def read_single_csvf(filename):
    """
     Reads and plot an ncol CSV file contained in
     filename and plots the column corresponding to ncol-icol
    """
#
    path = '../data/'    
    DATOS=np.genfromtxt(path + filename +'.csv',delimiter=',')
    nrow , ncol = np.shape(DATOS)
    X = DATOS[:,0]
    Y = DATOS[:,1]
    return X , Y
#
def loadcsv(filename , label):
#
    path1 = '../data/'
    path2 = '../results/'    
    DATOS=np.genfromtxt(path1 + filename +'.csv',delimiter=',')
    np.savetxt(path2 + label , DATOS)
#
    return DATOS
#
def grafdata(datay, datax , var , labelx , labely , title , Ngra):
    """
     Plots the acceleration history U[ndats] into
     Ngra. The plot is also stored into var.pdf
    """
    path = '../results/'
    var1=path + var + '.pdf'
    plt.figure(Ngra)
    plt.plot(datax , datay)
    plt.grid()
    plt.title(title)
    plt.xlabel(labelx)
    plt.ylabel(labely)
    plt.savefig(var1)
#
    return
#
def grafmany(signals , nt , nr , ifig):
    """
     Given nr signals each one of nt samples
     plot each signal and its FAS.
    """
    for jj in range(0 , nr):
        grafsignalG(signals[:, jj] , 'test' , 'Velocidad' , 'm/s' , -0.003 , 0.003 , 0.01 , ifig - 3)
        ifig =ifig + 1
#        x, Samag , A , nfs = rdc.Ftrans(signals[:, jj] , nt , 0.01 , 10.0)
#        grafFourier(Samag , x , nfs , 'FAS', 0.0 , 0.0 , 0.0 , 0.0 , ifig)
#        ifig =ifig + 1

def plot_txt_waves(filename , dt):
#
    path = '../data/'     
    DATOS=np.loadtxt(path + filename +'.txt' )
    nt = len(DATOS)
    Y = np.arange(0,nt*dt,dt) 
    XH = DATOS[:,0]
    XV = DATOS[:,1]
    
    plt.plot(Y , XH,)
    plt.grid()
    plt.figure(0)
    
    plt.plot(Y , XV)
    plt.grid()
    plt.figure(1)
    
    return 





































