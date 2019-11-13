from numpy import hstack
from numpy import zeros
from numpy import ones, random
from numpy.random import rand
from numpy.random import randn
from keras.models import Sequential
from keras.layers import Dense, Dropout, Input, Concatenate
from matplotlib import pyplot
from keras import optimizers
from keras.models import Model


import pickle



def define_pmodel(cvar_dim, n_outputs=3, mc=True):
    input1= Input(shape=(cvar_dim,))
    layer2=Dense(10, activation='sigmoid',
                 kernel_initializer='he_uniform', input_dim=cvar_dim)(input1)
    layer2=Dropout(0.1)(layer2, training=mc)
    layer3=Dense(10, activation='sigmoid')(layer2)
    layer3=Dropout(0.1)(layer3, training=mc)
    #layer4=Dense(10, activation='sigmoid')(layer3)
    #layer4=Dropout(0.2)(layer4, training=True)
    output=Dense(n_outputs, activation='linear')(layer3)
    model=Model(inputs=input1, outputs=output)
    model.compile(loss='mean_squared_error',optimizer='adam')
    return model

def generate_real_samples(X,Y,n):
    # generate inputs in [-0.5, 0.5]
    nx,nc=X.shape
    #print(nx,nc)
    ind=random.choice(nx,n)
    class_real=ones((n,1))
    return X[ind,:]+randn(n,X.shape[1])/8.,Y[ind,:]+randn(n,Y.shape[1])/8.

[X,Y,Xv,Yv]=pickle.load(open('deepL_Data.pklz','rb'))
xm=X.mean(axis=0)
xstd=X.std(axis=0)
for i in range(18):
    X[:,i]-=xm[i]
    Xv[:,i]-=xm[i]
    X[:,i]/=xstd[i]
    Xv[:,i]/=xstd[i]
X=X[:,:10]
Xv=Xv[:,:10]
cvar_dim=10
Y=Y[:,4:6]
Yv=Yv[:,4:6]
p_model=define_pmodel(cvar_dim, n_outputs=2, mc=True)

n_epochs=30000
n_batch=128
n_eval=1000

half_batch=int(n_batch/2)
from numpy import *
for i in range(n_epochs):
    # prepare real samples
    x_real, y_real = generate_real_samples(X,Y,n_batch)
    # prepare fake examples1
    p_model.train_on_batch(x_real, y_real)
    if i%n_eval==0:
        x_real, y_real = generate_real_samples(Xv,Yv,10*n_batch)
        yp=p_model.predict(x_real)
        print(corrcoef(yp[:,-1],y_real[:,-1]))
