import tensorflow as tf
from tensorflow import keras
from numpy import hstack
from numpy import zeros
from numpy import ones, random
from numpy.random import rand
from numpy.random import randn
from keras.models import Sequential
from keras.layers import Dense, Dropout, Input, Concatenate, Conv2D, Flatten, Activation, Lambda, \
    MaxPooling2D, Conv2DTranspose
from matplotlib import pyplot
from keras import optimizers
from keras.models import Model
from keras import backend as K

# Define parameters
batch_size = 128


def encoder( nx,  ny, nz, latent_dim):
    activation = tf.nn.relu
    X=Input(shape=(nx,ny,nz,))
    x = Conv2D(filters=16, kernel_size=2, strides=1, padding='same', activation=activation)(X)
    x = MaxPooling2D(pool_size=(2, 2))(x)
    #x = Conv2D(filters=3, kernel_size=2, strides=1, padding='same', activation=activation)(x)
    #x = MaxPooling2D(pool_size=(2, 2))(x)
    x = Flatten()(x)
    
    # Local latent variables
    mean_ = Dense(latent_dim)(x)
    std_dev = Dense(latent_dim, activation='softplus')(x)  # softplus to force >0

    def sampling(args):
        z_mean, z_log_sigma = args
        epsilon = K.random_normal(shape=(latent_dim,),
                                  mean=0., stddev=1.0)
        return z_mean + z_log_sigma * epsilon

    z = Lambda(sampling, output_shape=(latent_dim,))([mean_,std_dev])

#    epsilon = tf.random.normal(tf.stack([tf.shape(x)[0], latent_dim]), name='epsilon')
#    z = mean_ + tf.multiply(epsilon, std_dev)
    enc=Model(inputs=X,outputs=[z,mean_,std_dev])
    return enc
nl=4
enc=encoder(6,6,1,nl)
#stop
def decoder(latent_dim,nx,ny,nz):
    activation = tf.nn.relu
    inputLayer=Input(shape=(latent_dim,))
    nx2=int(nx/2)
    ny2=int(ny/2)
    def reshapeL(layer):
         layer = tf.reshape(layer, [-1, nx-1, ny-1, 1])
         return layer
     
    x = Dense(nx2*ny2,activation='relu')(inputLayer)
    x = Dense((nx-1)*(ny-1),activation='relu')(inputLayer)
    x = Lambda(reshapeL)(x)
    x=Conv2DTranspose(filters=1,kernel_size=2, strides=(1, 1), padding='valid', activation='relu',dilation_rate=(1, 1))(x)
    m=Model(inputs=inputLayer,outputs=x)
    return m

dec=decoder(nl,6,6,1)

def t_model(encoder,decoder,nx,ny,nz):
    input_batch= Input(shape=(nx,ny,nz,))
    z, mean_, std_dev = encoder(input_batch)
    output = decoder(z)
    
    flat_output = tf.reshape(output, [-1, nx * ny *nz])
    flat_input = tf.reshape(input_batch, [-1, nx * ny * nz])

    model=Model(inputs=input_batch,outputs=output)

    kl_loss = 0.5 * tf.reduce_sum(tf.square(mean_) + tf.square(std_dev) -\
                                          tf.compat.v1.log(tf.square(std_dev)) - 1, 1)
    
    #img_loss = tf.math.reduce_euclidean_norm(flat_input  - flat_output, 1)

    #vae_loss = K.mean(latent_loss + img_loss)
    #xent_loss = K.sum(K.binary_crossentropy(x_in, x_out), axis=[1, 2, 3])
    #p = K.clip(K.softmax(logits, -1), K.epsilon(), 1 - K.epsilon())
    # 假设先验分布为均匀分布，那么kl项简化为负熵
    #kl_loss = K.sum(p * K.log(p), axis=[1, 2])
    #vae_loss = K.mean(xent_loss + kl_loss)

    # add_loss是新增的方法，用于更灵活地添加各种loss
    #model.add_loss(vae_loss)
    return model, kl_loss

nubf_m,vae_kl_loss=t_model(enc,dec,6,6,1)

def reconstruction_loss(x, t_decoded):
    #return K.sum(K.binary_crossentropy(
    #    K.batch_flatten(x), 
    #    K.batch_flatten(t_decoded)), axis=-1)
    return K.sum(K.square(
        K.batch_flatten(x)-
        K.batch_flatten(t_decoded)), axis=-1)

def vae_loss(x, t_decoded):
    '''Total loss for the plain VAE'''
    return K.mean(reconstruction_loss(x, t_decoded) + vae_kl_loss)

nubf_m.compile(optimizer='adam', loss=vae_loss)


import pickle
#[X,Y,Xv,Yv]=pickle.load(open('deepL_Data_June11vMay23.pklz','rb'))
[z1dku_obs,z1dka_obs,z1dku_eff,z1dka_eff,z1dkuHR_obs,z1dkuHR_eff,\
 z1dkaHR_obs,z1dkaHR_eff,qrL,piaLKu,piaLKa,zsfcL,kext3dKu,kext3dKa,\
 z1dka_ms_obs,windL,qnubf]=pickle.load(open('zProfsq0125_11June2014_17:48.pklz','rb'))
[z1dku_obsv,z1dka_obsv,z1dku_effv,z1dka_effv,z1dkuHR_obsv,z1dkuHR_effv,\
 z1dkaHR_obsv,z1dkaHR_effv,qrLv,piaLKuv,piaLKav,zsfcLv,kext3dKuv,kext3dKav,\
 z1dka_ms_obsv,windLv,qnubfv]=pickle.load(open('zProfsq0125_11June2014_19:00.pklz','rb'))

x_input=[]
from numpy import *
for x in qnubf:
    x_input.append(x[0:1,:,:].T)
for x in qnubfv:
    x_input.append(x[0:1,:,:].T)
x_input=array(x_input)
x_input/=8.
def generate_real_samples(X,n):
    # generate inputs in [-0.5, 0.5]
    nt,nx,ny,nz=X.shape
    #print(nx,nc)
    ind=random.choice(nt,n)
    return X[ind,:,:]

n_epochs=10000
n_batch=128
n_eval=1000

half_batch=int(n_batch/2)
from numpy import *
for i in range(n_epochs):
    if i%100==0:
        print(i)
    # prepare real samples
    x_real = generate_real_samples(x_input,n_batch)
    # prepare fake examples1
    nubf_m.train_on_batch(x_real,x_real)
    if i%100==0:
        xp=nubf_m.predict(x_real)
        print(corrcoef(xp[:,3,3,0],x_real[:,3,3,0]))
        print(corrcoef(xp[:,0,0,0],x_real[:,0,0,0]))

    #if i%n_eval==0:
        
