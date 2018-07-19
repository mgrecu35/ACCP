import cPickle as pickle
import h5py as h5
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
from numpy import *
matplotlib.rcParams.update({'font.size': 14})
x=exp(-3+arange(40)*0.1)
def makeCFAD(pwcL,nodeL):
    cfad=zeros((40,40),float)
    n=pwcL.shape[0]
    for i in range(n):
        if pwcL[i,nodeL[i][1]]>0:
            for k in range(nodeL[i][0],88):
                if pwcL[i,k]>0:
                    ik=int((log(pwcL[i,k])+3)*10)
                    ir=87-k
                    if ik>=0 and ik<40 and ir<40:
                        cfad[ir,ik]+=1
    return cfad

def makeCFADZ(Z,nodeL):
    cfad=zeros((60,45),float)
    n=Z.shape[0]
    for i in range(n):
        if Z[i,nodeL[i][1]]>0:
            for k in range(nodeL[i][0],88):
                if Z[i,k]>0:
                    ik=int(Z[i,k])
                    ir=87-k
                    if ik>=0 and ik<45 and ir<60:
                        cfad[ir,ik]+=1
    return cfad


from numpy import *
from matplotlib.colors import LogNorm


def readdata2(i0,fname):
    [pRateL,paramDSDL,\
     pwcL,nodeL,nwNodeL,log10NwL]=pickle.load(open(fname,'rb'))
    
    
    #i0=0
    pRateL=array(pRateL[i0])
    nwNodeL=array(nwNodeL[i0])
    log10NwL=array(log10NwL[i0])
    pRatem=zeros((176),float)
    for i in range(176):
        a=nonzero(pRateL[:,i]>-0.001)
        pRatem[i]=pRateL[a[0],i].mean()

    #plt.plot(pRatem,(arange(176)*0.125)[::-1])
    #plt.ylim(0,10)

    #plt.figure()
    pwcL=array(pwcL[i0])
    nodeL=array(nodeL[i0])
    pRatem=zeros((88),float)
    dnL=[]
    
    for i in range(pwcL.shape[0]):
        pwcL[i,nodeL[i,4]+1:]=pwcL[i,nodeL[i,4]]
        dn1=zeros((88),float)
        dn1[range(nodeL[i,0],nodeL[i,4]+1)]=interp(arange(nodeL[i,0],nodeL[i,4]+1),nwNodeL[i,:],log10NwL[i,:])
        dn1[nodeL[i,4]:]=dn1[nodeL[i,4]]
        dnL.append(dn1)
    #pwcL[pwcL<0]=0
    #for i in range(88):
    #a=nonzero(pwcL[:,i]>-0.001)
    #pRatem[i]=pwcL[a[0],i].mean()
    print array(dnL).shape, array(pwcL).shape
    
    return pwcL,nodeL,array(dnL)
    plt.plot(pRatem,(arange(88)*0.25)[::-1])
    plt.ylim(0,10)
    
    plt.figure()
    cfad=makeCFAD(pwcL,nodeL)

    plt.pcolormesh(x,arange(40)*0.25,cfad,norm=LogNorm(),cmap='jet')
    plt.xscale('log')


def readdata(i0):
    [pRateL,paramDSDL,\
     pwcL,nodeL]=pickle.load(open('DPRProfsConv21300.pklz','rb'))
    
    
    i0=0
    pRateL=array(pRateL[i0])
    pRatem=zeros((176),float)
    for i in range(176):
        a=nonzero(pRateL[:,i]>-0.001)
        pRatem[i]=pRateL[a[0],i].mean()

    #plt.plot(pRatem,(arange(176)*0.125)[::-1])
    #plt.ylim(0,10)

    #plt.figure()
    pwcL=array(pwcL[i0])
    nodeL=array(nodeL[i0])
    pRatem=zeros((88),float)
    for i in range(pwcL.shape[0]):
        pwcL[i,nodeL[i,4]+1:]=pwcL[i,nodeL[i,4]]
    #pwcL[pwcL<0]=0
    #for i in range(88):
    #a=nonzero(pwcL[:,i]>-0.001)
    #pRatem[i]=pwcL[a[0],i].mean()
    return pwcL,nodeL
    plt.plot(pRatem,(arange(88)*0.25)[::-1])
    plt.ylim(0,10)
    
    plt.figure()
    cfad=makeCFAD(pwcL,nodeL)

    plt.pcolormesh(x,arange(40)*0.25,cfad,norm=LogNorm(),cmap='jet')
    plt.xscale('log')

def plot_zCFAD(cfad,fname,xmax,hmax):
    plt.figure()
    plt.pcolormesh(arange(45),arange(60)*0.25,cfad,norm=LogNorm(),cmap='jet')
    plt.xlim(0,xmax)
    plt.ylim(0,hmax)
    plt.xlabel('dBZ')
    plt.ylabel('Height(km)')
    plt.colorbar()
    plt.savefig(fname)

def saveD(sfcPwcL,z2d_ms,piaWL,pwcL,fname):
    pickle.dump([sfcPwcL,z2d_ms,piaWL,pwcL],open(fname,'wb'))

def loadD(fname):
    [sfcPwcL,z2d_ms,piaWL,pwcL]=pickle.load(open(fname,'rb'))
    return sfcPwcL,z2d_ms,piaWL,pwcL
