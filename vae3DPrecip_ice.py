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

qtL=qsL+qgL+qhL

X=zMSL[a[0],30:]+random.randn(n1,90)
w=zeros((7,7),float)
X[:,0]=piaL[a[0],0]+random.rand(n1)*4
for i in range(7):
    for j in range(7):
        dx=(i-3)*0.5/2.2
        dy=(j-3)*0.5/2.2
        w[i,j]=exp(-4*log(2)*(dx**2+dy**2))**2
wt=sum(w)
Y=[]
for i in a[0]:
    y=[]
    for k in range(4,38):
        y.append(sum(qtL[i,:,:,k]*w/wt))
    y.extend(piaL[i,:]/50.)
    y.append(piaL[i,0]/piaL[i,1])
    Y.append(y)
    
    
Y=array(Y)

Xv=zMSL[b[0],30:]+random.randn(n2,90)
Xv[:,0]=piaL[b[0],0]+random.rand(n2)*4
Yv=[]
for i in b[0]:
    y=[]
    for k in range(4,38):
        y.append(sum(qtL[i,:,:,k]*w)/wt)
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




corrcoef_wopia=[0.4954362575744359, 0.5820357853306992, 0.6621209764991621, 0.7202952442574447, 0.7634616029318976, 0.7984862644958252, 0.8288389216553678, 0.8556002107060687, 0.8902612168114912, 0.9200853647366091, 0.9415281272027706, 0.9514540371056689, 0.9565373482073456, 0.959924252789471, 0.9633677513552807, 0.9666515657526683, 0.9698769780798527, 0.9733211366285381, 0.9771250761360798, 0.9797216258732434, 0.9817687116832428, 0.9834866515620272, 0.9840682867061673, 0.9847955650247638, 0.9850576078881734, 0.9854038339715032, 0.9862325971826897, 0.9863145306886922, 0.9876980230859771, 0.9880367823596887, 0.9895180501752703, 0.9894907187182402, 0.9902625624907277, 0.9900259052877853, 0.7634427417191134, 0.8524960955572056, 0.7169000020553452]

nmrse=(y_mean-Yv).std(axis=0)/Yv.std(axis=0)

NMRSE_wopia=array([0.87482667, 0.815555  , 0.7501521 , 0.69388668, 0.64597523,
                    0.60209961, 0.55953394, 0.51764373, 0.45546771, 0.39180997,
                    0.33723746, 0.30826944, 0.29227369, 0.28116723, 0.26945907,
                    0.2578005 , 0.24567978, 0.23190607, 0.21497745, 0.20216109,
                    0.19121168, 0.18192014, 0.17862835, 0.1746713 , 0.17328227,
                    0.17123178, 0.16611944, 0.16527355, 0.15660546, 0.1543387 ,
                    0.14463105, 0.14496509, 0.13969121, 0.1413795 , 0.64707413,
                    0.52286956, 0.69774422])

corrcoef_pia=[0.5484390415393415, 0.6256724483936378, 0.6963610164159777, 0.7454227016390631, 0.7831946620292328, 0.8158563222102962, 0.8451193005917939, 0.8704817089029699, 0.8990517457120405, 0.9239737763134759, 0.9437284551210164, 0.9532119180620061, 0.9581615314888468, 0.961340214414287, 0.9644515790165227, 0.9672384861508156, 0.9699775193564112, 0.97305660707398, 0.9768452684001268, 0.9796168200755738, 0.9819183431902033, 0.9837806494807064, 0.9844525290746865, 0.9851593448815086, 0.9853495636623842, 0.985712867466324, 0.986492026374648, 0.986524540656558, 0.987945518188427, 0.988269651024845, 0.9896663565385903, 0.9895598056975701, 0.9902632029141086, 0.990011755872055, 0.991710535182077, 0.8936670550728032, 0.802911980557977]
NMRSE_pia=array([0.84204642, 0.78346677, 0.71900931, 0.66701462, 0.62197553,
                 0.57834827, 0.5345992 , 0.49220078, 0.43785413, 0.38251707,
                 0.33099229, 0.30276413, 0.28685379, 0.27617936, 0.2653922 ,
                 0.25537243, 0.2450137 , 0.23270785, 0.21598361, 0.20248745,
                 0.19031227, 0.18022028, 0.1764006 , 0.17249904, 0.17149601,
                 0.16938603, 0.16460603, 0.16412966, 0.15513665, 0.15289894,
                 0.14362807, 0.14442856, 0.13954148, 0.14126409, 0.1285386 ,
                 0.44887666, 0.59609836])
