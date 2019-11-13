from netCDF4 import Dataset
from numpy import *
import matplotlib.pyplot as plt
icase=2
if icase==1:
    fname_in='wrfJune11_2014.Fields_G4ICE_17:48:00.nc'
if icase==2:
    fname_in='wrfJune11_2014.Fields_G4ICE_19:00:00.nc'
if icase==3:
    fname_in='wrfJune11_2014.Fields_G4ICE_18:24:00.nc'
    
fh=Dataset(fname_in)

z3d_obs=fh['z3d_obs'][:,:,:,:]
z3d=fh['z3d'][:,:,:,:]
v3d=fh['v3d'][:,:,:,:]
v3d_air=fh['wm'][:,:,:]
kext3d=fh['kext3d'][:,:,:,:]
salb3d=fh['salb3d'][:,:,:,:]
asym3d=fh['asym3d'][:,:,:,:]
pia2d=fh['pia2d'][:,:,:]
h=fh['z'][:,:,:]
T=fh['T'][:,:,:]
print(h.shape)
print(T.shape)
print(pia2d.shape)
qr=fh['qr'][:,:,:]
#import julia
#jl=julia.Julia()
#jl.include("retr.jl")
#jl.init()

a=nonzero(pia2d[:,:,1]>10)


from numpy import *
w=zeros((11,11),float)
for i in range(11):
    for j in range(11):
        dx=(i-5)*0.5/2.2
        dy=(j-5)*0.5/2.2
        w[i,j]=exp(-4*log(2)*(dx**2+dy**2))**2
wt=sum(w)


z1dku=[]
z1dka=[]
z1dku_obs=[]
z1dka_obs=[]
z1dka_ms_obs=[]
z1dku_eff=[]
z1dka_eff=[]
z1dkaHR_obs=[]
z1dkuHR_obs=[]
z1dkaHR_eff=[]
z1dkuHR_eff=[]
qrL=[]
hg=arange(120)*0.125+1
piaLKu=[]
piaLKa=[]
zsfcL=[]
kextLKa=[]
kextLKu=[]

theta=0.7/4.
noNorm=0
dr=0.125
alt=400.
freq=35.

import pyHB2 as pyHB2

#for i1,j1 in zip(a[0],a[1]):

z3d_obs[z3d_obs<-60]=-60.
z3d[z3d<-60]=-60.
nx,ny,nz,nf=z3d.shape
z3dms=zeros((nx,ny,120))-99.
doppler_v3d_Ka=zeros((nx,ny,120))-99.
v3d_Ka=zeros((nx,ny,120))-99.
windsL=[]
for j1 in range(ny):
    print(j1)
    for i1 in range(nx):
        h1=h[:,j1,i1]
        #print(h1[30:50])
        hm=0.5*(h1[1:]+h1[:-1])
        zg=interp(hg,hm,z3d[i1,j1,:,1])
        zg_obs=interp(hg,hm,z3d_obs[i1,j1,:,1])[::-1]
        kext1=kext3d[i1,j1,:,1]
        salb1=salb3d[i1,j1,:,1]
        asym1=asym3d[i1,j1,:,1]
        kextg=interp(hg,hm,kext1)
        salbg=interp(hg,hm,salb1)
        asymg=interp(hg,hm,asym1)
        noMS=0
        zms = pyHB2.multiscatterf(kextg[::-1],salbg[::-1],asymg[::-1],\
                                  zg[::-1],dr,noMS,alt,
                                  theta,freq,noNorm)
        a1=nonzero(zg[::-1]<-40)
        if i1>1 and i1<nx-1 and j1>1 and j1<ny-1:
            v1d=(v3d[i1-1:i1+2,j1-1:j1+2,:,1].mean(axis=0)).mean(axis=0)
        else:
            v1d=v3d[i1,j1,:,1]
        if i1>1 and i1<nx-1 and j1>1 and j1<ny-1:
            v1d_air=(v3d_air[i1-1:i1+2,j1-1:j1+2,:].mean(axis=0)).mean(axis=0)
        else:
            v1d_air=v3d_air[i1,j1,:]
        v1dg=interp(hg,hm,v1d)
        v1d_airg=interp(hg,h1,v1d_air)
        doppler_v3d_Ka[i1,j1,:]=v1dg
        v3d_Ka[i1,j1,:]=v1d_airg
            #if zms[a1].max()>0:
        #    stop
        #print(zms-zg[::-1])
        #print(zg[::-1])
        #noMS=1
        #znoms = pyHB2.multiscatterf(kextg[::-1],salbg[::-1],asymg[::-1],\
        #                            zg[::-1],dr,noMS,alt,
        #                            theta,freq,noNorm)
        z3dms[i1,j1,:]=zms.copy()


import pickle
#pickle.dump(z3dms,open('z3dms_19:00.pklz','wb'))
#z3dms=pickle.load(open('z3dms.pklz','rb'))
#stop

for i1,j1 in zip(a[0],a[1]):
    if i1>=5 and i1<nx-5 and j1>=5 and j1<ny-5:
        zku1d=zeros(100)-99
        zku1de=zeros(100)-99
        zka1d=zeros(100)-99
        zka1dms=zeros(120)-99
        zka1de=zeros(100)-99
        kext1d=zeros((100,2))
        for k in range(100):
            if z3d_obs[i1-5:i1+6,j1-5:j1+6,k,0].max()>-90:
                zm=sum(10**(0.1*z3d_obs[i1-5:i1+6,j1-5:j1+6,k,0])*w)/wt
                zku1d[k]=log10(zm)*10.0
                zme=sum(10**(0.1*z3d[i1-5:i1+6,j1-5:j1+6,k,0])*w)/wt
                zku1de[k]=log10(zme)*10.0
                kext1d[k,0]=sum(kext3d[i1-5:i1+6,j1-5:j1+6,k,0]*w)/wt
                kext1d[k,1]=sum(kext3d[i1-5:i1+6,j1-5:j1+6,k,1]*w)/wt
            if z3d_obs[i1-5:i1+6,j1-5:j1+6,k,1].max()>-90:
                zm=sum(10**(0.1*z3d_obs[i1-5:i1+6,j1-5:j1+6,k,1])*w)/wt
                zka1d[k]=log10(zm)*10.0
                zme=sum(10**(0.1*z3d[i1-5:i1+6,j1-5:j1+6,k,1])*w)/wt
                zka1de[k]=log10(zme)*10.0

        for k in range(120):
            if z3dms[i1-5:i1+6,j1-5:j1+6,k].max()>-90:
                zm=sum(10**(0.1*z3dms[i1-5:i1+6,j1-5:j1+6,k])*w)/wt
                zka1dms[k]=log10(zm)*10.0
            else:
                zka1dms[k]=-60.
            #print(zka1dms[k])

        h1=h[:,j1,i1]
        hm=0.5*(h1[1:]+h1[:-1])
        if h[0,j1,i1]<1:
            piaSum=10**(-0.1*pia2d[i1-5:i1+6,j1-5:j1+6,0])*w/wt
            piaW=-log10(sum(piaSum))*10.0
            piaLKu.append(piaW)
            piaSum=10**(-0.1*pia2d[i1-5:i1+6,j1-5:j1+6,1])*w/wt
            piaW=-log10(sum(piaSum))*10.
            piaLKa.append(piaW)
            ind=nonzero((hm-1)*(hm-2.)<0)
            qr1=mean([sum(qr[ik,j1-5:j1+6,i1-5:i1+6]*w)/wt for ik in ind[0]])
            ind=nonzero((hm-2)*(hm-3.)<0)
            qr2=mean([sum(qr[ik,j1-5:j1+6,i1-5:i1+6]*w)/wt for ik in ind[0]])
            ind=nonzero((hm-3)*(hm-4.)<0)
            qr3=mean([sum(qr[ik,j1-5:j1+6,i1-5:i1+6]*w)/wt for ik in ind[0]])
            qrL.append([qr1,qr2,qr3])
            zsfcL.append([zku1d[0],zka1d[0],zku1de[ind].mean(),zka1de[ind].mean()])
            z1dku_obs.append(interp(hg,hm,zku1d))
            z1dka_obs.append(interp(hg,hm,zka1d))
            z1dku_eff.append(interp(hg,hm,zku1de))
            z1dka_eff.append(interp(hg,hm,zka1de))
            z1dkaHR_obs.append(interp(hg,hm,z3d_obs[i1,j1,:,1]))
            z1dkuHR_obs.append(interp(hg,hm,z3d_obs[i1,j1,:,0]))
            z1dkaHR_eff.append(interp(hg,hm,z3d[i1,j1,:,1]))
            z1dkuHR_eff.append(interp(hg,hm,z3d[i1,j1,:,0]))
            kextLKu.append(interp(hg,hm,kext1d[:,0]))
            kextLKa.append(interp(hg,hm,kext1d[:,1]))
            z1dka_ms_obs.append(zka1dms[::-1])
            windsL.append([doppler_v3d_Ka[i1,j1,:],v3d_Ka[i1,j1,:]])
            #plt.plot(z1dka_obs[-1],hg)
            #plt.plot(z1dka_ms_obs[-1],hg)
            #plt.show()
            #stop
        #print(zku1d)
        #print(h[:,j1,i1])
        #print(T[:,i1,j1])
        #stop

        t1=T[:,j1,i1]
        #jl.retrieve(zku1d,h1,t1)
        #stop


z1dku_obs=array(z1dku_obs)
z1dka_obs=array(z1dka_obs)
z1dka_ms_obs=array(z1dka_ms_obs)
z1dku_eff=array(z1dku_eff)
z1dka_eff=array(z1dka_eff)
z1dkaHR_obs=array(z1dkaHR_obs)
z1dkuHR_obs=array(z1dkuHR_obs)
z1dkaHR_eff=array(z1dkaHR_eff)
z1dkuHR_eff=array(z1dkuHR_eff)
piaLKu=array(piaLKu)
piaLKa=array(piaLKa)
import pickle

if icase==1:
    fname_out='zProfsq0125_17:48.pklz'
if icase==2:
    fname_out='zProfsq0125_19:00.pklz'
if icase==3:
    fname_out='zProfsq0125_18:24.pklz'

pickle.dump([z1dku_obs,z1dka_obs,z1dku_eff,z1dka_eff,z1dkuHR_obs,z1dkuHR_eff,\
             z1dkaHR_obs,z1dkaHR_eff,qrL,piaLKu,piaLKa,zsfcL,array(kextLKu),\
             array(kextLKa),z1dka_ms_obs,windsL],open(fname_out,'wb'))
stop
import matplotlib
matplotlib.rcParams.update({'font.size': 13})

plt.plot(10*log10((10**(0.1*z1dkuHR_obs)).mean(axis=0)),hg,label='Ku fine')
plt.plot(10*log10((10**(0.1*z1dkaHR_obs)).mean(axis=0)),hg,label='Ka fine')
plt.plot(10*log10((10**(0.1*z1dku_obs)).mean(axis=0)),hg,label='Ku coarse')
plt.plot(10*log10((10**(0.1*z1dka_obs)).mean(axis=0)),hg,label='Ka coarse')
plt.ylim(0,12)
plt.xlim(0,50)
plt.xlabel("dBZ")
plt.ylabel("Height (km)")
plt.legend()
plt.savefig('nubfZ2_g.png')

from sklearn.cluster import KMeans
z1dku_obs[z1dku_obs<0]=0
nc=10
kmeans = KMeans(n_clusters=nc, random_state=1).fit(z1dku_obs[:,0:95])
for i in range(nc):
    plt.figure()
    a=nonzero(kmeans.labels_==i)
    plt.plot(z1dku_obs[a[0],:].mean(axis=0),hg)
    plt.title('Class %2i'%(i))
    zm=z1dku_obs[a[0],:].mean(axis=0)
    zm2=z1dkuHR_obs[a[0],:].mean(axis=0)
    zs=z1dku_obs[a[0],:].std(axis=0)
    plt.ylabel('Height(km)')
    plt.xlabel('dBZ')
    plt.errorbar(zm, hg, xerr=zs, color='blue')
    plt.errorbar(zm2, hg, xerr=zs, color='green')
    #print (rL[a[0],:].mean(axis=0))
    plt.ylim(0,12)
    plt.xlim(10,50)
    plt.xlabel("dBZ")
    plt.ylabel("Height (km)")
    plt.savefig('class_1748_g_%2.2i.png'%i)
    plt.close('all')

