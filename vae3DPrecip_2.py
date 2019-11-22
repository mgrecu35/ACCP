from numpy import hstack
from numpy import zeros
from numpy import ones, random
from numpy.random import rand
from numpy.random import randn
import matplotlib.pyplot as plt
import pickle
# Define parameters
batch_size = 128
	
from numpy.random import seed
from numpy import *
import pickle


qrL,qsL,qgL,qhL,zMSL,piaL=pickle.load(open('q3d_zm_NUBF_f2.pklz','rb'))
seed(2771973225)
nt=qrL.shape[0]
r=random.random(nt)
a=nonzero(r<0.2)
b=nonzero(r>0.2)
n1=a[0].shape[0]
n2=b[0].shape[0]
zMSL[zMSL<12]=0

X=zMSL[a[0],30:]+random.randn(n1,90)
w=zeros((7,7),float)
#X[:,0]=piaL[a[0],0]+random.rand(n1)*4
for i in range(7):
    for j in range(7):
        dx=(i-3)*0.5/2.2
        dy=(j-3)*0.5/2.2
        w[i,j]=exp(-4*log(2)*(dx**2+dy**2))**2
wt=sum(w)
Y=[]
for i in a[0]:
    y=[]
    for k in range(0,14):
        y.append(sum(qrL[i,:,:,k]*w/wt))
    y.extend(piaL[i,:]/50.)
    y.append(piaL[i,0]/piaL[i,1])
    Y.append(y)
    
    
Y=array(Y)

Xv=zMSL[b[0],30:]+random.randn(n2,90)
#Xv[:,0]=piaL[b[0],0]+random.rand(n2)*4
Yv=[]
for i in b[0]:
    y=[]
    for k in range(14):
        y.append(sum(qrL[i,:,:,k]*w)/wt)
    y.extend(piaL[i,:]/50.)
    y.append(piaL[i,0]/piaL[i,1])
    Yv.append(y)
    
Yv=array(Yv)


from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel

kernel = 1.0 * RBF(length_scale=100.0, length_scale_bounds=(1e-2, 1e3)) \
    + WhiteKernel(noise_level=1, noise_level_bounds=(1e-10, 1e+1))
gp2 = GaussianProcessRegressor(kernel=kernel,
                              alpha=0.0).fit(X, Y)

y_mean, y_cov = gp2.predict(Xv, return_cov=True)
c=[corrcoef(y_mean[:,i],Yv[:,i])[0,1] for i in range(Yv.shape[1])]

from sklearn.neighbors import NearestNeighbors
nbrs = NearestNeighbors(n_neighbors=30, algorithm='ball_tree').fit(Y[:])
distances, indices = nbrs.kneighbors(y_mean[:])
nubf_factor=[]
for i,pia2 in zip(indices,piaL[b]):
    nubf_factor.append([mean(piaL[a[0][i][0:10],0]/piaL[a[0][i][0:10],1]),pia2[0]/pia2[1]])

nubf_factor=array(nubf_factor)
print(corrcoef(nubf_factor.T))

cfad_nubf=zeros((40,40),float)
for i in range(nubf_factor.shape[0]):
    i0=int(nubf_factor[i,0]*40)
    j0=int(nubf_factor[i,1]*40)
    cfad_nubf[i0,j0]+=1

import matplotlib

matplotlib.rcParams.update({'font.size': 13})

fig=plt.figure()
ax = fig.add_subplot(111,aspect='equal')
ax.set_aspect('equal')
cfad_nubfm=ma.array(cfad_nubf,mask=cfad_nubf<=0)
plt.pcolormesh(arange(40)/40.,arange(40)/40.,cfad_nubfm, cmap='jet')
plt.xlabel("True NUBF factor")
plt.ylabel("Estimated NUBF factor")
plt.savefig('nubfFactor_wopia.png')


NMRSE_pia=array([0.68251261, 0.65411281, 0.61761185, 0.58164088, 0.55024405,
                 0.52088946, 0.49204292, 0.46341001, 0.43444081, 0.4062908 ,
                 0.38319957, 0.3718346 , 0.37435627, 0.39045523, 0.14322331,
                 0.44259635, 0.58870795])
corrcoef_pia=[0.7323324870751698, 0.7576055063091602, 0.7871794443239651, 0.8137619389321328, 0.8351466268179478, 0.8536699727393021, 0.8705730539333423, 0.8861547564823182, 0.9007553413521351, 0.9138684516691825, 0.9239048805983566, 0.9287144412728238, 0.9279239471758639, 0.9219749491903266, 0.9897235070316489, 0.8967220281105653, 0.808444840489621]

corrcoef_wopia=[0.6394076926439117, 0.6764901352518913, 0.7198423538485385, 0.7592631533127842, 0.7916758623554808, 0.8200319842061442, 0.846293539385844, 0.868743173534532, 0.8873612540385395, 0.9026006166612514, 0.9133485754426953, 0.9179072029447671, 0.9164107754730408, 0.9049775003537401, 0.7705582873390772, 0.8570626038958451, 0.7243256138246357]


NMRSE_wopia=array([0.77046169, 0.73748841, 0.69460327, 0.65092259, 0.61097096,
                   0.57231782, 0.53275638, 0.49540132, 0.46134818, 0.43091097,
                   0.40780718, 0.39766384, 0.40139683, 0.42754218, 0.6377759 ,
                   0.51521253, 0.6896317 ])


NMRSE_wpia=(y_mean-Yv).std(axis=0)/Yv.std(axis=0)


NMRSE_wopia=(y_mean-Yv).std(axis=0)/Yv.std(axis=0)
