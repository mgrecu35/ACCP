import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import pyHB2 as pyHB2
pyHB2.initt()
icase=1
if icase==1:
    fh=Dataset('wrfout_d04_2014-06-11_17:48:00')
if icase==2:
    fh=Dataset('wrfout_d04_2014-06-11_19_00_00')
if icase==3:
    fh=Dataset('wrfout_d04_2014-06-11_18:24:00')
if icase==4:
    fh=Dataset('wrfout_d04_2014-06-12_18:00:00')
if icase==5:
    fh=Dataset('wrfout_d04_2014-05-23_21:36:00')

qr=fh['QRAIN'][0,0,500:700,50:350]
from numpy import log
qrm=np.ma.array(qr,mask=qr<1e-6)
#plt.pcolormesh(np.log10(qrm*1e3),vmin=-3,vmax=0,cmap='jet')

#stop
lonlat=np.loadtxt('lonlats.txt')


if icase==1:
    ny1=0
    ny2=300
    nx1=870
    nx2=930

#----------------19:00:##########
if icase==2:
    ny1=350
    ny2=750
    nx1=35
    nx2=296
################################
#----------------18:24:##########

if icase==3:
    ny1=350
    ny2=900
    nx1=35
    nx2=296
    
if icase==4:
    ny1=500
    ny2=900
    nx1=400
    nx2=600

if icase==5:
    ny1=500
    ny2=700
    nx1=50
    nx2=350
#ny1=0
#ny2=300
#nx1=850
#nx2=950
#ny1=550
#ny2=850
#nx1=135
#nx2=136
import matplotlib
from readWRF import *
qv,qr,qs,qc,qg,qh,ncr,ncs,ncg,\
    nch,th,prs,T,t2c,h,xlong,xlat,rho=readWRF(fh,ny1,ny2,nx1,nx2)

from numpy import log10

#qs=qs+qh
#qh*=0.
nyX=0
hm=0.5*(h[:-1,50,nyX]+h[1:,50,nyX])
#plt.subplot(411)
#plt.pcolormesh(np.arange(300),hm,qr[:,0:300,nyX],norm=matplotlib.colors.LogNorm(),vmin=0.01e-3,vmax=6*1e-3,cmap='jet')
#plt.subplot(412)
#plt.pcolormesh(np.arange(300),hm,qh[:,0:300,nyX],norm=matplotlib.colors.LogNorm(),vmin=0.01e-3,vmax=6*1e-3,cmap='jet')
#plt.subplot(413)
#plt.pcolormesh(np.arange(300),hm,qg[:,0:300,nyX],norm=matplotlib.colors.LogNorm(),vmin=0.01e-3,vmax=6*1e-3,cmap='jet')
#plt.subplot(414)
#plt.pcolormesh(np.arange(300),hm,qs[:,0:300,nyX],norm=matplotlib.colors.LogNorm(),vmin=0.01e-3,vmax=6*1e-3,cmap='jet')

w=fh['W'][0,:,ny1:ny2,nx1:nx2]
wm=w[:,:,:]
nz,ny,nx=qv.shape

n0w=0.
nfreq=8
zKu_L=[]
zKa_L=[]
zW_L=[]
vKu_L=[]
vKa_L=[]
vW_L=[]
j=200



x1L=[]
from getNw import *
dr=0.25

fd=Dataset('nw_dm_dsd/dm_z.nc')
dm_z=fd['dm_z'][:,:]

fd=Dataset('nw_dm_dsd/doppler_2.nc')
vt_R=fd['rain'][:,:]
vt_S=fd['snow'][:,:]

import matplotlib.pyplot as plt
import matplotlib.colors as col

freqs=[13.8,35.5,94.]
ny2-=ny1
ny1-=ny1
nx2-=nx1
nx1-=nx1
from scipy.ndimage.filters import gaussian_filter
fshape=qh.shape
from numpy import random
dnh=random.randn(fshape[0],fshape[1],fshape[2])
dns=random.randn(fshape[0],fshape[1],fshape[2])
dng=random.randn(fshape[0],fshape[1],fshape[2])
dnr=random.randn(fshape[0],fshape[1],fshape[2])

dnri=gaussian_filter(dnr,sigma=2)*5
dngi=gaussian_filter(dns,sigma=2)*5
dnsi=gaussian_filter(dng,sigma=2)*5
dnhi=gaussian_filter(dnh,sigma=2)*5

#stop


from fields3D_WRF4ICE import *
z3d,z3d_obs,kext3d,salb3d,asym3d,\
    v3d,dnr,dns,dng,\
    rainDBL,snowDBL,pia2d=radarFields_3d(nx1,nx2,ny1,ny2,qs,qg,\
                                         qh,qr,qc,qv,T,prs,ncs,ncg,nch,ncr,\
                                         rho,wm,h,dm_z,vt_R,vt_S,pyHB2,nz,freqs,\
                                         dnri,dnsi,dngi,dnhi)

wmT=wm[:,ny1:ny2,nx1:nx2].T
import xarray as xr
z3dx=xr.DataArray(z3d,dims=['nx','ny','nz','nfreq'])
wmTx=xr.DataArray(wmT,dims=['nx','ny','nz1'])
kext3dx=xr.DataArray(kext3d,dims=['nx','ny','nz','nfreq'])
salb3dx=xr.DataArray(salb3d,dims=['nx','ny','nz','nfreq'])
asym3dx=xr.DataArray(asym3d,dims=['nx','ny','nz','nfreq'])
v3dx=xr.DataArray(v3d,dims=['nx','ny','nz','nfreq'])
z1=h[:,ny1:ny2,nx1:nx2]
temp=T[:,ny1:ny2,nx1:nx2]
zx=xr.DataArray(z1,dims=['nz1','ny','nx'])
tempx=xr.DataArray(temp,dims=['nz','ny','nx'])
pia2dx=xr.DataArray(pia2d,dims=['nx','ny','nz3'])
qrx=xr.DataArray(qr[:,ny1:ny2,nx1:nx2]*rho[:,ny1:ny2,nx1:nx2]*1e3,dims=['nz','ny','nx'])
qsx=xr.DataArray(qs[:,ny1:ny2,nx1:nx2]*rho[:,ny1:ny2,nx1:nx2]*1e3,dims=['nz','ny','nx'])
qgx=xr.DataArray(qg[:,ny1:ny2,nx1:nx2]*rho[:,ny1:ny2,nx1:nx2]*1e3,dims=['nz','ny','nx'])
qcx=xr.DataArray(qc[:,ny1:ny2,nx1:nx2]*rho[:,ny1:ny2,nx1:nx2]*1e3,dims=['nz','ny','nx'])

dnrx=xr.DataArray(dnr,dims=['nx','ny','nz'])
dnsx=xr.DataArray(dns,dims=['nx','ny','nz'])
dngx=xr.DataArray(dng,dims=['nx','ny','nz'])

z3d_obsx=xr.DataArray(z3d_obs,dims=['nx','ny','nz','nfreq'])
d=xr.Dataset({'z3d':z3dx,'z3d_obs':z3d_obsx,'kext3d':kext3dx,'salb3d':salb3dx,'asym3d':asym3dx,\
              'v3d':v3dx,'z':zx, 'T': tempx,\
              'qr':qrx, 'qs':qsx, 'qg':qgx, 'qc':qcx, 'dnr':dnrx, \
              'dns':dnsx, 'dng':dngx, 'pia2d':pia2dx, 'wm':wmTx})

if icase==1:
    fname_out='wrfJune11_2014.Fields_G4ICE_17:48:00.nc'
if icase==2:
    fname_out='wrfJune11_2014.Fields_G4ICE_19:00:00.nc'
if icase==3:
    fname_out='wrfJune11_2014.Fields_G4ICE_18:24:00.nc'
if icase==4:
    fname_out='wrfJune12_2014.Fields_G4ICE_18:00:00.nc'
if icase==5:
    fname_out='wrfMay23_2014.Fields_G4ICE_21:36:00.nc'
d.to_netcdf(fname_out)
stop
from numpy import *
w=zeros((11,11),float)
for i in range(11):
    for j in range(11):
        dx=(i-5)*0.5/4.5
        dy=(j-5)*0.5/4.5
        w[i,j]=exp(-4*log(2)*(dx**2+dy**2))**2
w=sum(w)
a=nonzero(pia2d[:,:,1]>30)
z1dku=[]
z1dka=[]
z1dku_obs=[]
z1dka_obs=[]
for i1,j1 in zip(a[0],a[1]):
    if i1>=5 and i1<294 and j1>=5 and j1<54:
        zku1d=zeros(100)-99
        zka1d=zeros(100)-99
    z1dku_obs.append(zku1d)
    z1dka_obs.append(zka1d)

z1dku_obs=array(z1dku_obs)
z1dka_obs=array(z1dka_obs)
plt.plot(z1dku_obs.mean(axis=0),h[1:,0,0])
plt.plot(z1dka_obs.mean(axis=0),h[1:,0,0])

 

#plt.figure()
#plt.contourf(arange(60),z[1:,60,60],z3d[:,34,:,0].T,vmin=0,levels=arange(13)*5,cmap='jet')
#plt.subplot(211)
z3dm=ma.array(z3d_obs,mask=z3d_obs<0)
plt.subplot(211)
plt.pcolormesh(arange(ny1,ny2),h[1:,100,0],z3dm[0,:,:,0].T,vmin=0,vmax=45,cmap='jet')
plt.ylim(0,15)
plt.colorbar()
plt.subplot(212)
plt.pcolormesh(arange(ny1,ny2),h[1:,100,0],z3dm[0,:,:,1].T,vmin=0,vmax=40,cmap='jet')
plt.ylim(0,15)
plt.colorbar()
#plt.show()

a=nonzero(pia2d[:,:,1]>30)

stop
import bokeh.io
import bokeh.plotting
import bokeh.transform
import bokeh.util.hex

bins = bokeh.util.hex.hexbin( z3dm[:,:,:,1].T.flatten(),\
                              h[1:,:,nx1:nx2].flatten(),\
                              0.5)

from bokeh.models import Range1d, LogColorMapper, LogTicker, ColorBar
p = bokeh.plotting.figure(plot_width=700,plot_height=700,\
                          tools="wheel_zoom,reset,save", match_aspect=True, \
                          x_range=Range1d(10, 50), y_range=Range1d(0,15),\
                          background_fill_color="white", aspect_ratio=1)

p.xaxis.axis_label_text_font_size = "20pt"
p.yaxis.axis_label_text_font_size = "20pt"
p.yaxis.major_label_text_font_size = "20pt"
p.xaxis.major_label_text_font_size = "20pt"
p.grid.visible = False
p.title.text = " "
p.title.align = "right"
p.title.text_color = "orange"
p.title.text_font_size = "25px"
p.xaxis.axis_label = 'dBZ'
p.yaxis.axis_label = 'Height(km)'
color_mapper = LogColorMapper(palette="Spectral11", low=1, \
                              high=max(bins.counts))
from bokeh.transform import linear_cmap, log_cmap
linearc=log_cmap('counts', 'Spectral11', 1, max(bins.counts))

p.hex_tile(q="q", r="r", size=0.5, line_color=None, source=bins,
           fill_color=linearc)

color_mapper = LogColorMapper(palette="Spectral11", low=1, high=max(bins.counts))
color_bar = ColorBar(color_mapper=color_mapper, ticker=LogTicker(),
                     label_standoff=12, border_line_color=None, location=(0,0))
color_bar.major_label_text_font_size="20pt"
color_bar.title="Counts"
color_bar.title_text_font_size="20pt"
p.add_layout(color_bar, 'right')

bokeh.io.export_png(p, filename="cfadKa_att.png")
bokeh.io.output_file("hex_tile.html")
bokeh.io.show(p)



bins2 = bokeh.util.hex.hexbin( z3dm[:,:,:,0].T.flatten(),\
                              h[1:,:,nx1:nx2].flatten(),\
                              0.5)

from bokeh.models import Range1d, LogColorMapper, LogTicker, ColorBar

p2 = bokeh.plotting.figure(plot_width=700,plot_height=700,\
                          tools="wheel_zoom,reset,save", match_aspect=True, \
                          x_range=Range1d(10, 50), y_range=Range1d(0,15),\
                          background_fill_color="white", aspect_ratio=1)

p2.xaxis.axis_label_text_font_size = "20pt"
p2.yaxis.axis_label_text_font_size = "20pt"
p2.yaxis.major_label_text_font_size = "20pt"
p2.xaxis.major_label_text_font_size = "20pt"
p2.grid.visible = False
p2.title.text = " "
p2.title.align = "right"
p2.title.text_color = "orange"
p2.title.text_font_size = "25pt"
p2.xaxis.axis_label = 'dBZ'
p2.yaxis.axis_label = 'Height(km)'
color_mapper = LogColorMapper(palette="Spectral11", low=1, \
                              high=max(bins.counts))
from bokeh.transform import linear_cmap, log_cmap
linearc=log_cmap('counts', 'Spectral11', 1, max(bins.counts))

p2.hex_tile(q="q", r="r", size=0.5, line_color=None, source=bins2,
           fill_color=linearc)

color_mapper = LogColorMapper(palette="Spectral11", low=1, high=max(bins.counts))
color_bar = ColorBar(color_mapper=color_mapper, ticker=LogTicker(),
                     label_standoff=12, border_line_color=None, location=(0,0))
color_bar.major_label_text_font_size="20pt"
color_bar.title="Counts"
color_bar.title_text_font_size="20pt"
p2.add_layout(color_bar, 'right')

bokeh.io.export_png(p2, filename="cfadKu_att.png")
bokeh.io.output_file("hex_tile.html")
bokeh.io.show(p2)
