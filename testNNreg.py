from sklearn.neighbors import NearestNeighbors,KNeighborsRegressor
from readPickles import *
from sklearn.tree import DecisionTreeRegressor
from sklearn import linear_model
from numpy import *
from sklearn.ensemble import GradientBoostingRegressor

est = GradientBoostingRegressor(n_estimators=100, learning_rate=0.1,
                                max_depth=1, random_state=0, loss='ls')
sfcPwcL,z2d_ms,piaWL,pwcL=loadD("dBase_CSAT_0.pklz")
X=z2d_ms[:,57:86][:,::-1]
dx=random.randn(X.shape[0],X.shape[1])*1
X+=dx
X[X<-5]=-5
rL=[]
r=random.random(X.shape[0])
a=nonzero(r<0.75)
b=nonzero(r>0.75)
nbrs = NearestNeighbors(n_neighbors=50, \
                        algorithm='ball_tree').fit(X[a[0],:])


zCSATL,hL=pickle.load(open('CPRCase.pklz','rb'))
pwcRet=[]
distL=[]
zwcL=[]
from hb import *
dr=0.24

piawL=[]
i=80
#zwc,piaW=hb(zCSATL[i,:],hL[i,:],dr)
#stop
i=100
#zwc,piaW,sfcPwc,sfcPwc2=hb(zCSATL[i,:],hL[i,:],dr)
#stop
sfcPwcL=[]
sfcPwcL2=[]
for i in range(zCSATL.shape[0]):
    zwc,piaW,sfcPwc,sfcPwc2=hb(zCSATL[i,:],hL[i,:],dr)
    zwcL.append(zwc)
    piawL.append(piaW)
    sfcPwcL.append(sfcPwc)
    sfcPwcL2.append(sfcPwc2)
zwcL=array(zwcL)
zmm=ma.array(zwcL,mask=zwcL<-20)
plt.figure()
x=array([arange(zwcL.shape[0]) for k in range(125)]).T
ax=plt.subplot(211)
plt.pcolormesh(x,hL,zmm[:,:], \
               vmin=-10,vmax=25,cmap='jet')
plt.ylim(0,10)

zmm2=ma.array(zCSATL,mask=zwcL<-20)

x=array([arange(zwcL.shape[0]) for k in range(125)]).T
ax=plt.subplot(212)
plt.pcolormesh(x,hL,zmm2[:,:], \
               vmin=-10,vmax=25,cmap='jet')
plt.ylim(0,10)
plt.figure()
plt.plot(piawL)
plt.figure()
plt.plot(sfcPwcL2)
plt.plot(sfcPwcL)
stop
for it in range(5):
    
    neigh = KNeighborsRegressor(n_neighbors=100,weights='distance')
    neigh.fit(X[a[0],:],(sfcPwcL[a]))
    yp=neigh.predict(X[b[0],:])
    print corrcoef(yp,(sfcPwcL[b]))
    rL.append(corrcoef(yp,(sfcPwcL[b]))[1,0])
plt.plot(yp,(sfcPwcL[b]),'+')
#plt.xlim(0.01,3)
#plt.ylim(0.01,3)
