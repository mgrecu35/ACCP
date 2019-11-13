import pickle
from numpy import *
[z1dku_obs,z1dka_obs,z1dku_eff,z1dka_eff,z1dkuHR_obs,z1dkuHR_eff,\
 z1dkaHR_obs,z1dkaHR_eff,qrL,piaLKu,piaLKa,zsfcL,kext3dKu,kext3dKa,\
 z1dka_ms_obs,windL]=pickle.load(open('zProfsq0125_18:24.pklz','rb'))
[z1dku_obsv,z1dka_obsv,z1dku_effv,z1dka_effv,z1dkuHR_obsv,z1dkuHR_effv,\
 z1dkaHR_obsv,z1dkaHR_effv,qrLv,piaLKuv,piaLKav,zsfcLv,kext3dKuv,kext3dKav,\
 z1dka_ms_obsv,windLv]=pickle.load(open('zProfsq0125_19:00.pklz','rb'))
hg=arange(120)*0.125+1
[z1dku_obsv2,z1dka_obsv2,z1dku_effv2,z1dka_effv2,z1dkuHR_obsv2,z1dkuHR_effv2,\
 z1dkaHR_obsv2,z1dkaHR_effv2,qrLv2,piaLKuv2,piaLKav2,zsfcLv2,kext3dKuv2,\
 kext3dKav2,z1dka_ms_obsv2,windLv2]=pickle.load(open('zProfsq0125_17:48.pklz','rb'))

windL=array(windL)
windLv=array(windLv)
windLv2=array(windLv2)

from sklearn.cluster import KMeans
from sklearn import neighbors

n_neighbors=30

X=[[z1[2],pia+random.randn()*3,z1[30],z1[20],z1[50],z1[60],0.001*sum(10**(0.1*z1[2:30])),0.001*sum(10**(0.1*z1[30:60]))] for z1,pia in zip(z1dka_ms_obs,piaLKa)]
X=[]
nlev=2.0
for z1,pia,wind in zip(z1dka_ms_obs,piaLKa,windL):
    x1=[z1[2],pia+random.randn()*3,z1[30],z1[20],z1[50],z1[60],0.001*sum(10**(0.1*z1[2:30])),0.001*sum(10**(0.1*z1[30:60]))]
    x1.extend(wind[0,30:50][::2]+nlev*random.randn(10))
    X.append(x1)
    
for z1,pia,wind in zip(z1dka_ms_obsv,piaLKav,windLv):
    x1=[z1[2],pia+random.randn()*3,z1[30],z1[20],z1[50],z1[60],0.001*sum(10**(0.1*z1[2:30])),0.001*sum(10**(0.1*z1[30:60]))]
    x1.extend(wind[0,30:50][::2]+nlev*random.randn(10))
    X.append(x1)

Y=[]
for z1,wind,zs in zip(qrL,windL,zsfcL):
    x=list(z1)
    x.append(zs[-1])
    x.extend(wind[1,30:50][::2])
    Y.append(x)
for z1,wind,zs in zip(qrLv,windLv,zsfcLv):
    x=list(z1)
    x.append(zs[-1])
    x.extend(wind[1,30:50][::2])
    Y.append(x)



Xv=[]
    
for z1,pia,wind in zip(z1dka_ms_obsv2,piaLKav2,windLv2):
    x1=[z1[2],pia+random.randn()*3,z1[30],z1[20],z1[50],z1[60],0.001*sum(10**(0.1*z1[2:30])),0.001*sum(10**(0.1*z1[30:60]))]
    x1.extend(wind[0,30:50][::2]+nlev*random.randn(10))
    Xv.append(x1)

Yv=[]

for z1,wind,zs in zip(qrLv2,windLv2,zsfcLv2):
    x=list(z1)
    x.append(zs[-1])
    x.extend(wind[1,30:50][::2])
    Yv.append(x)

    


X=array(X)
Y=array(Y)
Yv=array(Yv)
Xv=array(Xv)

import pickle

n=len(X)
r=random.random(n)
a=nonzero(r<0.58)
b=nonzero(r>0.58)

Xv[:,:8][Xv[:,0:8]<0]=0
X[:,0:8][X[:,0:8]<0]=0
pickle.dump([X,Y,Xv,Yv],open('deepL_Data.pklz','wb'))
stop
xstd=X.std(axis=0)
xm=X.mean(axis=0)
for i in range(8):
    X[:,i]-=xm[i]
    Xv[:,i]-=xm[i]
    X[:,i]/=xstd[i]
    Xv[:,i]/=xstd[i]
knn = neighbors.KNeighborsRegressor(n_neighbors, weights='distance')
knn.fit(X[a[0],:], Y[a])
y=knn.predict(X[b[0],:])
yv=knn.predict(Xv[:,:])

stop
nc=50
zL=[]
qL=[]
for z1,q in zip(z1dka_ms_obs,qrL):
    zL.append(z1)
    qL.append(q)
for z1,q in zip(z1dka_ms_obsv,qrLv):
    zL.append(z1)
    qL.append(q)

zL=array(zL)
qL=array(qL)
zL[zL<0]=0
z1dka_ms_obsv2[z1dka_ms_obsv2<0]=0
kmeans = KMeans(n_clusters=nc, random_state=1).fit(zL)
labels = kmeans.labels_
qc=zeros((nc))
zcL=[]
for i in range(nc):
    a=nonzero(labels==i)
    zcL.append(zL[a[0],:].mean(axis=0))
    qc[i]=qL[a].mean()

classv=kmeans.predict(z1dka_ms_obsv2)

qrp=qc[classv]
print(corrcoef(qrp,qrLv2))

import matplotlib.pyplot as plt
qrp=qc[classv]
plt.plot(z1dka_ms_obsv2[0,:],hg)
plt.plot(zcL[classv[0]],hg)
n_neighbors=20
knn = neighbors.KNeighborsRegressor(n_neighbors, weights='distance')
knn.fit(zL[:,0:90],qL)
yp=knn.predict(z1dka_ms_obsv[:,0:90])
