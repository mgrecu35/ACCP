from scipy.special import gamma
from numpy import log10
pi=3.1415
rhosn=200.
rhog=400.
rhow=1000.
cs=rhosn*pi/6.
ds=3.
cg=rhog*pi/6
dg=3.
cons1=gamma(1+ds)*cs
cons2=gamma(1+dg)*cg

def getZ_W(dm,zKu_S,nws,dm_z,zWS):
    i0dm=int((dm-0.2)/0.02)
    if i0dm>0 and i0dm<190:
        dZ=zKu_S-10*log10(nws/0.08)-dm_z[i0dm,4]
        zWS=dm_z[i0dm,6]+dZ+10*log10(nws/0.08)
    return zWS

def getn0s(qs,ncs):
    lams1=(cons1*ncs/qs)**(1./ds)
    if lams1<500:
        lams1=500
    n0s=ncs*lams1
    return lams1, n0s

def getn0g(qg,ncg):
    lams2=(cons2*ncg/qg)**(1./dg)
    if lams2<500:
        lams2=500
    n0g=ncg*lams2
    return lams2, n0g

def getn0w(qr,ncr):
    lamr = (pi*rhow*ncr/qr)**(1./3.)
    if lamr<500:
        lamr=500
    n0w = ncr*lamr

    return lamr, n0w
