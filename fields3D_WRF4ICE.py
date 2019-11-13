from getNw import *
from numpy import *
def radarFields_3d(nx1,nx2,ny1,ny2,qs,qg,qh,qr,qc,qv,T,prs,ncs,ncg,nch,ncr,rho,wm,z,dm_z,vt_R,vt_S,pyHB2,nz,freqs):
    nfreq=8
    z3d=zeros((nx2-nx1,ny2-ny1,nz,3),float)
    z3dm=zeros((nx2-nx1,ny2-ny1,nz,3),float)
    v3d=zeros((nx2-nx1,ny2-ny1,nz,3),float)
    kext3d=zeros((nx2-nx1,ny2-ny1,nz,3),float)
    salb3d=zeros((nx2-nx1,ny2-ny1,nz,3),float)
    pia3d=zeros((nx2-nx1,ny2-ny1,3),float)
    asym3d=zeros((nx2-nx1,ny2-ny1,nz,3),float)
    dns=zeros((nx2-nx1,ny2-ny1,nz),float)
    dnr=zeros((nx2-nx1,ny2-ny1,nz),float)
    dng=zeros((nx2-nx1,ny2-ny1,nz),float)
    qag=qh.copy()*0.
    ncag=nch.copy()*0.
    freqsT=[10.000000,19.000000,22.000000,37.000000,85.000000,165, 186.31,190.31]
    snowDBL=[]
    rainDBL=[]
    snowDBL=[]
    for j in range(ny1,ny2):
        print(j)
        for i in range(nx1,nx2):
            zKu_1d=zeros((nz),float)-99
            zKa_1d=zeros((nz),float)-99
            zW_1d=zeros((nz),float)-99
            zKu_1dm=zeros((nz),float)-99
            zKa_1dm=zeros((nz),float)-99
            zW_1dm=zeros((nz),float)-99
            vW_1d=zeros((nz),float)-99
            vKa_1d=zeros((nz),float)-99
            vKu_1d=zeros((nz),float)-99
            piaKu=0.
            piaKa=0.
            piaW=0.
            for k in range(nz-1,-1,-1):
                iwc=(qs[k,j,i]+qg[k,j,i]+qh[k,j,i]+qag[k,j,i])*rho[k,j,i]*1e3
                pwc=qr[k,j,i]*rho[k,j,i]*1e3
                dpwc=0.
                if k==nz-1:
                    dr=(z[nz-1,j,i]-z[nz-2,j,i])/1.
                else:
                    dr=(z[k+1,j,i]-z[k,j,i])/1.
                if iwc>1e-4:
                    t2c=T[k,j,i]-273.15
                    fract=1.0
                    if t2c>0:
                        f=1-t2c/3.
                        if fract<0.1:
                            fract=0.1
                        fract=0.01
                    if(t2c<-20):
                        t2c=-20
                    #fract=1.
                    nws=0.08*exp(-0.12*t2c)
                    nwg=0.08*exp(-0.12*t2c)
                    nwh=0.08#*exp(-0.12*t2c)
                    dns[i-nx1,j-ny1,k]=log10(nws/0.08)
                    dng[i-nx1,j-ny1,k]=log10(nwg/0.08)
                    iwc1=fract*(qs[k,j,i])*rho[k,j,i]*1e3
                    iwc2=fract*(qg[k,j,i])*rho[k,j,i]*1e3
                    iwc3=(qh[k,j,i])*rho[k,j,i]*1e3
                    #if t2c>0:
                    #    fract=0.5
                    #fract=1.
                    dpwc=(1-fract)*(qs[k,j,i]+qg[k,j,i])*rho[k,j,i]*1e3
                    #if t2c<-10:
                        #nws=0.08
                    kextS1,salbS1,asymS1,dpiaKu_S1,dpiaKa_S1,zKu_S1,zKa_S1,zW_S1,dmS1 = pyHB2.getsnowp2(log10(nws/0.08),(iwc1),nfreq)
                    if t2c>-10:
                        kextS2,salbS2,asymS2,dpiaKu_S2,dpiaKa_S2,zKu_S2,zKa_S2,zW_S2,dmS2 = pyHB2.getsnowp2(log10(nwg/0.08),(iwc2),nfreq)
                    else:
                        kextS2,salbS2,asymS2,dpiaKu_S2,dpiaKa_S2,zKu_S2,zKa_S2,zW_S2,dmS2 = pyHB2.getsnowp2(log10(nwg/0.08),(iwc2),nfreq)

                    kextS3,salbS3,asymS3,dpiaKu_S3,dpiaKa_S3,zKu_S3,zKa_S3,zW_S3,dmS3 = pyHB2.gethail(log10(nwh/0.08),(iwc3),nfreq)
                    zW_S1=getZ_W(dmS1,zKu_S1,nws,dm_z,zW_S1)
                    zW_S2=getZ_W(dmS2,zKu_S2,nwg,dm_z,zW_S2)
                    zW_S3=getZ_W(dmS3,zKu_S3,nwh,dm_z,zW_S3)
                    zKu_S=log10(10.**(0.1*zKu_S1)+10.**(0.1*zKu_S2)+10.**(0.1*zKu_S3))*10
                    zKa_S=log10(10.**(0.1*zKa_S1)+10.**(0.1*zKa_S2)+10.**(0.1*zKa_S3))*10
                    zWS=log10(10.**(0.1*zW_S1)+10.**(0.1*zW_S2)+10.**(0.1*zW_S3))*10
                    if iwc1>0.01:
                        snowDBL.append([dpiaKu_S1+dpiaKu_S2,dpiaKa_S1+dpiaKa_S2,
                                        zKu_S,zKa_S,iwc1+iwc2])
                    i0s1=int((zKu_S1-10*log10(nws/0.08)+12.)/0.5)
                    i0s1=max(i0s1,0)
                    i0s1=min(114,i0s1)
                    i0s2=int((zKu_S2-10*log10(nwg/0.08)+12.)/0.5)
                    i0s2=max(i0s2,0)
                    i0s2=min(114,i0s2)
                    i0s3=int((zKu_S3-10*log10(nwh/0.08)+12.)/0.5)
                    i0s3=max(i0s3,0)
                    i0s3=min(114,i0s3)
                    vs1=vt_S[i0s1,1:]
                    vs2=vt_S[i0s2,1:]
                    vs3=vt_S[i0s3,1:]
                    z_S1=array([zKu_S1,zKa_S1,zW_S1])
                    z_S2=array([zKu_S2,zKa_S2,zW_S2])
                    z_S3=array([zKu_S3,zKa_S3,zW_S3])
                    vS=(10**(0.1*z_S1)*vs1+10**(0.1*z_S2)*vs2+\
                        10**(0.1*z_S3)*vs3)/\
                        (10**(0.1*z_S1)+10**(0.1*z_S2)+10**(0.1*z_S3))
                
                    dpiaW_S1=kextS1[4]*4.343
                    dpiaW_S2=kextS2[4]*4.343
                    dpiaW_S3=kextS3[4]*4.343

                else:
                    zKu_S=-99
                    zKa_S=-99
                    zWS=-99
                    vS=array([0.,0.,0.])
                    dpiaKu_S1=0
                    dpiaKu_S2=0
                    dpiaKu_S3=0
                    dpiaKa_S1=0
                    dpiaKa_S2=0
                    dpiaKa_S3=0
                    dpiaW_S1=0
                    dpiaW_S2=0
                    dpiaW_S3=0
                    kextS1=zeros((nfreq),float)
                    kextS2=zeros((nfreq),float)
                    kextS3=zeros((nfreq),float)
                    kextS4=zeros((nfreq),float)
                    asymS1=zeros((nfreq),float)
                    asymS2=zeros((nfreq),float)
                    asymS3=zeros((nfreq),float)
                    asymS4=zeros((nfreq),float)
                    salbS1=zeros((nfreq),float)
                    salbS2=zeros((nfreq),float)
                    salbS3=zeros((nfreq),float)
                    salbS4=zeros((nfreq),float)
                n0w=0.
                if pwc+dpwc>1e-9:
                    t2c=T[k,j,i]-273.15
                    if t2c<-5:
                        t2c=-5
                    if t2c<10:
                        nwr=0.04#*exp((10-t2c)/6.)
                    else:
                        nwr=0.04
                    dnr[i-nx1,j-ny1,k]=log10(nwr/0.08)
                    
                    kextR,salbR,asymR,dpiaKu_R,dpiaKa_R,zKu_R,zKa_R,zWR, dmR = pyHB2.getrainp2(log10(nwr/0.08),(pwc+dpwc),nfreq)
                    x=[dpiaKu_R,dpiaKa_R,zKu_R,zKa_R,pwc]
                    rainDBL.append(x)
                    i0dm=int((dmR-0.2)/0.02)
                    if i0dm>0 and i0dm<190:
                        dZ=zKu_R-10*log10(nwr/0.08)-dm_z[i0dm,1]
                        zWR=dm_z[i0dm,3]+dZ+10*log10(nwr/0.08)
                    i0r=int((zKu_R-10*log10(nwr/0.08)+12)/0.5)
                    i0r=max(i0r,0)
                    i0r=min(124,i0r)
                    vR=vt_R[i0r,1:]
                    if zKu_R>300:
                        stop
                    dpiaW_R=kextR[4]*4.343
                else:
                    zKu_R=-99
                    zKa_R=-99
                    zWR=-99
                    dpiaKu_R=0
                    dpiaKa_R=0.
                    dpiaW_R=0.
                    dmR=0
                    vR=array([0,0.,0])
                    kextR=zeros((nfreq),float)
                    salbR=zeros((nfreq),float)
                    asymR=zeros((nfreq),float)
                a=2.65
                b=0.098
                dpiaL=[]
                ireturn=0
                asymg=(asymS1*salbS1*kextS1+asymS2*salbS2*kextS2+asymS3*salbS3*kextS3+asymS4*salbS4*kextS4+asymR*salbR*kextR)/\
                       (salbS1*kextS1+salbS2*kextS2+salbS3*kextS3+salbS4*kextS4+salbR*kextR+1e-5)
                salbg=(salbS1*kextS1+salbS2*kextS2+salbS3*kextS3+salbS4*kextS4+salbR*kextR)/\
                       (kextS1+kextS2+kextS3+kextS4+kextR+1e-5)
                #if zKu_R<0 and zKu_S>20:
                    #print(asymg)
                    #print(salbg)
                    #print(salbS1)
                    #print(zKu_S1,zKa_S1,zW_S1)
                    #print(salbS2)
                    #print(zKu_S2,zKa_S2,zW_S2)
                    #print(salbS3)
                    #print(zKu_S3,zKa_S3,zW_S3)
                    #print(salbS4)
                    #stop
                asymRadFreq=interp(freqs,freqsT,asymg)
                salbRadFreq=interp(freqs,freqsT,salbg)
                kextL=[]
                for freq in freqs:
                    cldw=qc[k,j,i]*1e3*rho[k,j,i]
                    z_clw = pyHB2.gcloud(freq,T[k,j,i],cldw)
                    absair,abswv = pyHB2.gasabsr98(freq,T[k,j,i],qv[k,j,i],prs[k,j,i],ireturn)
                    dpiaL.append((z_clw+absair+abswv)*4.343)
                    kextL.append((z_clw+absair+abswv))
                zKu_1d[k]=log10(10.**(0.1*zKu_S)+10.**(0.1*zKu_R))*10-piaKu-(dpiaKu_S1+dpiaKu_S2+dpiaKu_S3+dpiaKu_R+dpiaL[0])*dr
                zKa_1d[k]=log10(10.**(0.1*zKa_S)+10.**(0.1*zKa_R))*10-piaKa-(dpiaKa_S1+dpiaKa_S2+dpiaKa_S3+dpiaKa_R+dpiaL[1])*dr
                zW_1d[k]=log10(10.**(0.1*zWS)+10.**(0.1*zWR))*10.-piaW-(dpiaW_S1+dpiaW_S2+dpiaW_S3+dpiaW_R+dpiaL[2])*dr
                piaKu+=2*(dpiaKu_S1+dpiaKu_S2+dpiaKu_R+dpiaL[0])*dr  # attenuation due to cloud and  water vapor
                piaKa+=2*(dpiaKa_S1+dpiaKa_S2+dpiaKa_R+dpiaL[1])*dr  # not included yet
                piaW+=2*(dpiaW_S1+dpiaW_S2+dpiaW_R+dpiaL[2])*dr
                z3d[i-nx1,j-ny1,k,0]=log10(10.**(0.1*zKu_S)+10.**(0.1*zKu_R))*10
                z3d[i-nx1,j-ny1,k,1]=log10(10.**(0.1*zKa_S)+10.**(0.1*zKa_R))*10
                z3d[i-nx1,j-ny1,k,2]=log10(10.**(0.1*zWS)+10.**(0.1*zWR))*10
               
                kext3d[i-nx1,j-ny1,k,0]=(dpiaKu_S1+dpiaKu_S2+dpiaKu_S3+\
                                         dpiaKu_R+dpiaL[0])/4.343
                kext3d[i-nx1,j-ny1,k,1]=(dpiaKa_S1+dpiaKa_S2+dpiaKa_S3+\
                                         dpiaKa_R+dpiaL[1])/4.343
                kext3d[i-nx1,j-ny1,k,2]=(dpiaW_S1+dpiaW_S2+dpiaW_S3+\
                                         dpiaW_R+dpiaL[2])/4.343
                if piaKu<0:
                    stop
                vT=(10**(0.1*zKu_S)*vS+10**(0.1*zKu_R)*vR)/\
                    (10**(0.1*zKu_S)+10**(0.1*zKu_R))
                vW_1d[k]=vT[2]-wm[k,j,i]
                vKu_1d[k]=vT[0]-wm[k,j,i]
                vKa_1d[k]=vT[1]-wm[k,j,i]
                v3d[i-nx1,j-ny1,k,:]=vT-wm[k,j,i]
                asym3d[i-nx1,j-ny1,k,:]=asymRadFreq
                salb3d[i-nx1,j-ny1,k,:]=salbRadFreq
                #print(piaKu,piaKa,piaW)
            z3dm[i-nx1,j-ny1,:,0]=zKu_1d
            z3dm[i-nx1,j-ny1,:,1]=zKa_1d
            z3dm[i-nx1,j-ny1,:,2]=zW_1d
            pia3d[i-nx1,j-ny1,:]=[piaKu,piaKa,piaW]
    return z3d,z3dm,kext3d,salb3d,asym3d,v3d, dnr,dns,dng, rainDBL,snowDBL,pia3d

