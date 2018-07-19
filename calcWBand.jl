using PyCall
PyDict(pyimport("matplotlib")["rcParams"])["font.size"] = [12]
using PyPlot
unshift!(PyVector(pyimport("sys")["path"]), "")
@pyimport readPickles as r
@pyimport numpy as np
(pwcL,nodeL,dnL)=r.readdata2(0,"newDPRConv/DPRProfsConv21300.pklz")
include("SDSU_Wrap.jl")
nx,nz=size(pwcL);
imelt=0
whail=0.
wgrpl=0.
icldw=0.
fmelt=0.
temp=273.
z2d=Array{Float32}(nx,88)
z2d_noms=Array{Float32}(nx,88)
z2d_ms=Array{Float32}(nx,88)
dr=0.25
i4=4
# multiscatterF(nrange,kext,salb,g,zTrue,zMS,dr,noMS,&
#      alt,dr,theta,freq,noNorm,alt)
kextA=getKextWintGPM()
sfcPwcL=Array{Float32}(nx)
piaWL=Array{Float32}(nx)
println(sizeof(dnL))
dn1=-0.5
dn2=-0.3
zsL=Array{Float32}(0)
zmL=Array{Float32}(0)
zrL=Array{Float32}(0)
piasL=Array{Float32}(0)
piarL=Array{Float32}(0)
pwcR=Array{Float32}(0)
piamL=Array{Float32}(0)
zwL=Array{Float32}(0)
for i=1:nx
    z=Array{Float32}(88)*0-99
    piaW=0.
    zm=Array{Float32}(88)*0-99
    kext1d=Array{Float32}(88)*0
    for k=1:88
        kext1d[k]=kext1d[k]+Float32(kextA[89-k])
    end
    salb1d=Array{Float32}(88)*0
    g1d=Array{Float32}(88)*0
    if nodeL[i,1]-20<nodeL[i,2]
        println(i)
        sfcPwcL[i]=pwcL[i,nodeL[i,5]]
        for k=nodeL[i,1]+1:nodeL[i,2]
            #println(k)
            dn=dn1+(dnL[i,k]-6.9)
            wsnow=pwcL[i,k]/10^dn
            wrain=0.
            dBZ,kexttot,dBZm,kexttotm,
            salbtot,asymtot=getScatt(imelt,wrain,wsnow,
                                     wgrpl,whail,icldw,temp,fmelt);
            z[k]=dBZ[i4]+10*dn
            piaW+=kexttot[i4]*4.343*dr*10^dn
            zm[k]=z[k]-piaW
            piaW+=kexttot[i4]*4.343*dr*10^dn
            kext1d[k]=kext1d[k]+kexttot[i4]*10^dn
            salb1d[k]=salbtot[i4]
            g1d[k]=asymtot[i4]
            if(wsnow>0.0001)
                push!(zsL,z[k])
                push!(piasL,kexttot[i4]*10^dn*4.343)
            end
        end
        for k = nodeL[i,2]+1: nodeL[i,4]
            f=(nodeL[i,4]+1-k)/(nodeL[i,4]-nodeL[i,2]+0.)
            dn=(1-f)*dn2+f*dn1+(dnL[i,k]-6.9)
            wsnow=pwcL[i,k]*f/10^dn
            wrain=(1-f)*pwcL[i,k]/10^dn
            #println(f,"wsnow=",wsnow, "wrain=", wrain)
            dBZ,kexttot,dBZm,kexttotm,
            salbtot,asymtot=getScatt(imelt,wrain,wsnow,
                                     wgrpl,whail,icldw,temp,fmelt);
            z[k]=dBZ[i4]+10*dn
            piaW+=kexttot[i4]*4.343*dr*10^dn
            zm[k]=z[k]-piaW
            piaW+=kexttot[i4]*4.343*dr*10^dn
            kext1d[k]=kext1d[k]+kexttot[i4]*10^dn
            salb1d[k]=salbtot[i4]
            g1d[k]=asymtot[i4]
            if(pwcL[i,k]>0.0001)
                push!(zmL,z[k])
                push!(piamL,kexttot[i4]*10^dn*4.343)
            end

#println(dBZ)
        end
        for k = nodeL[i,4]+1:88#nodeL[i,5]+1
            dn=dn2+(dnL[i,k]-6.9)
            wrain=pwcL[i,k]/10^dn
            wsnow=0.
            dBZ,kexttot,dBZm,kexttotm,
            salbtot,asymtot=getScatt(imelt,wrain,wsnow,
                                     wgrpl,whail,icldw,temp,fmelt);
            z[k]=dBZ[i4]+10*dn
            piaW+=kexttot[i4]*4.343*dr*10^dn
            zm[k]=z[k]-piaW
            piaW+=kexttot[i4]*4.343*dr*10^dn
            kext1d[k]=kext1d[k]+kexttot[i4]*10^dn
            salb1d[k]=salbtot[i4]
            g1d[k]=asymtot[i4]
            if(pwcL[i,k]>0.0001)
                push!(zrL,z[k])
                push!(piarL,kexttot[i4]*10^dn*4.343)
                push!(pwcR,wrain)
                push!(zwL,z[k])
            end
        end
        noMS=Int32(0)
        freq=Float32(94.6)
        dr=Float32(0.25)
        noNorm=Int32(0)
        theta=Float32(0.108)
        zW=z
        #print(piaW,"  ")
        #println(typeof(kext1d))
        zMS=msCalc(kext1d,salb1d,g1d,zW,noMS,freq,theta,dr,noNorm)
        #msCalc(kext,salb,g,zw,noMS,freq,theta,dr,noNorm)
    end
    z2d[i,:]=z
    z2d_noms[i,:]=zm
    z2d_ms[i,:]=zMS
    piaWL[i]=piaW
end
zcfad=r.makeCFADZ(z2d,nodeL)
zmax=30
fname="zcfad.png"
hmax=15
r.plot_zCFAD(zcfad,fname,zmax+5,hmax)
zcfadm=r.makeCFADZ(z2d_noms,nodeL)
fname="zcfad_noMS.png"
r.plot_zCFAD(zcfadm,fname,zmax,hmax)
zcfadm=r.makeCFADZ(z2d_ms,nodeL)
fname="zcfadMS.png"
r.plot_zCFAD(zcfadm,fname,zmax,hmax)
r.saveD(sfcPwcL,z2d_ms,piaWL,pwcL,"dBase_CSAT_0.pklz")
@pyimport numpy as np
a=np.polyfit(0.1*zsL,log10.(piasL),1)
b=np.polyfit(0.1*zmL,log10.(piamL),1)
c=np.polyfit(0.1*zrL,log10.(piarL),1)
d=np.polyfit(log10.(piarL),log10.(pwcR),1)
e=np.polyfit(0.1*zwL,log10.(pwcR),1)
#dBZ,kexttot,dBZm,kexttotm=getScatt(imelt,wrain,wsnow,wgrpl,whail,icldw,temp,fmelt);
#push!(wrainL,wrain)
#    push!(dbzL,dBZ[4])
#    push!(kextL,kexttot[4])


