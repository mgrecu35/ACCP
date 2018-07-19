using PyCall
PyDict(pyimport("matplotlib")["rcParams"])["font.size"] = [12]
using PyPlot
@pyimport matplotlib.colors as col

unshift!(PyVector(pyimport("sys")["path"]), "")

#@pyimport testSDSU as sdsu
#@pyimport testSDSU as sdsu
#sdsu.init()
imd=2
ccall((:setparams_, "SDSU"), Void, () )    
ccall((:simradarlut_, "SDSU"), Void, (Ptr{Int32},), &imd)

function getKExtWint()
    kext=[0.08216629, 0.07987797, 0.07557384, 0.07255906, 0.06967789,
          0.06681726, 0.06410733, 0.0579389 , 0.0535338 , 0.05497366,
          0.04758603, 0.04699508, 0.04475723, 0.04301621, 0.0412752 ,
          0.03960225, 0.03795199, 0.03521832, 0.03287926, 0.03057081,
          0.02841735, 0.02642137, 0.02526997, 0.0246093 , 0.02303219,
          0.02217608, 0.0213438 , 0.02051153, 0.01967925, 0.01884698,
          0.01811531, 0.01741973, 0.01671269, 0.01577828, 0.01484388,
          0.01403029, 0.01323571, 0.0125283 , 0.011831  , 0.01113371,
          0.01050953, 0.01017607, 0.00990444, 0.00963282, 0.00926262,
          0.00873582, 0.00813506, 0.00680915, 0.00657739, 0.00638145,
          0.00615319, 0.00587261, 0.00559203, 0.00531145, 0.00503086,
          0.00458927, 0.00440378, 0.00420007, 0.00394702, 0.00380493,
          0.00365312, 0.00348672, 0.00336368, 0.00330825, 0.00325281,
          0.00319737, 0.00306326, 0.00292377, 0.00278428, 0.00264478,
          0.0025104 , 0.00241081, 0.00231121, 0.00223822, 0.00216396,
          0.00209095, 0.00201794, 0.00193494, 0.00183698, 0.00174985]
    return kext
end
function getKextWintGPM()
    kext=[8.21662946e-02, 7.55738447e-02, 6.96778881e-02, 6.41073307e-02,
          5.35338041e-02, 4.75860271e-02, 4.47572256e-02, 4.12751959e-02,
          3.79519929e-02, 3.28792570e-02, 2.84173511e-02, 2.52699735e-02,
          2.30321881e-02, 2.13438024e-02, 1.96792514e-02, 1.81153095e-02,
          1.67126890e-02, 1.48438794e-02, 1.32357142e-02, 1.18310039e-02,
          1.05095288e-02, 9.90444264e-03, 9.26262388e-03, 8.13505938e-03,
          6.57739045e-03, 6.15319080e-03, 5.59202721e-03, 5.03086362e-03,
          4.40378107e-03, 3.94701658e-03, 3.65311502e-03, 3.36368262e-03,
          3.25280970e-03, 3.06326061e-03, 2.78427622e-03, 2.51039903e-03,
          2.31121460e-03, 2.16396105e-03, 2.01793567e-03, 1.83697749e-03,
          1.66902493e-03, 1.54053277e-03, 1.43668544e-03, 1.30044650e-03,
          1.21366131e-03, 1.12104987e-03, 1.02332281e-03, 9.47407738e-04,
          8.86880292e-04, 8.26022914e-04, 7.72601122e-04, 7.19179330e-04,
          6.60842315e-04, 6.20694319e-04, 5.72564688e-04, 5.43907313e-04,
          5.14138487e-04, 4.76742716e-04, 4.38730887e-04, 3.84415069e-04,
          3.47759608e-04, 3.18882607e-04, 2.91012331e-04, 2.72114102e-04,
          2.58439039e-04, 2.43519336e-04, 2.27899523e-04, 2.12279710e-04,
          1.99132661e-04, 1.86132302e-04, 1.73918219e-04, 1.62477930e-04,
          1.48804159e-04, 1.35429434e-04, 1.25819754e-04, 1.18628387e-04,
          1.11888704e-04, 1.04919122e-04, 9.65958559e-05, 8.86443754e-05,
          8.07605225e-05, 7.34793581e-05, 6.79371083e-05, 6.66605336e-05,
          6.02330162e-05, 5.38578605e-05, 4.96013161e-05, 4.56343579e-05]
    return kext
end
function getcloud(freq,cldw)
    t=Float32(273.)
    z_clw=Array{Float32}(1)
    cldwi=Float32(cldw)
    freqi=Float32(freq)
    #print(z_clw)
    ccall((:gcloud_,"SDSU"),Void,(Ptr{Float32},Ptr{Float32},
                                  Ptr{Float32},Ptr{Float32}),
          &freqi,&t,&cldwi,z_clw)
    z_clw1=z_clw[1]
    #print(z_clw1)
    #print(x)
    #SUBROUTINE GCLOUD(FREQY,T,CLW,Z_CLW)
    return z_clw1
end

function msCalc(kext,salb,g,zw,noMS,freq,theta,dr,noNorm)
    pi4=pi^4
    lambd=300./freq/1e3
    lamb4=lambd^4
    Z=10.^(0.1*zw)*pi4/1e18/lamb4/4*0.93
    dZ=log10(pi4/1e18/lamb4/4*0.93)*10;
    #println(typeof(kext))
    kext=kext*1e-3
    #println(typeof(kext))
    if noNorm==1
        salb=salb*1e-3
        salb=salb/kext
    end
    ext2bscatt=Array{Float32}(88)*0.+1000.
    for k=1:88
        if zw[k]>-10
            ext2bscatt[k]=kext[k]/Z[k]
        else
            ext2bscatt[k]=100000.
        end
        #println(kext[k]," ",salb[k]," ", g[k], " ", zw[k])
    end
        
    alt=705.
    dr=dr*1000
    bscatt=zeros(Float32,88)
    #println(typeof(bscatt))
    #println(noMS,freq," ",theta)
    #println(typeof(kext))
    #println(typeof(ext2bs))
    nr=Int32(88)
    kext=Array{Float32}(kext)
    ext2bscatt=Array{Float32}(ext2bscatt)
    salb=Array{Float32}(salb)
    g=Array{Float32}(g)
    
    ccall((:multiscatterf_,"SDSU"),
          Void,(Ptr{Int32},Ptr{Float32},Ptr{Float32},Ptr{Float32},
                Ptr{Float32},Ptr{Float32},Ptr{Float32},
                Ptr{Int32},Ptr{Float32},Ptr{Float32},Ptr{Float32},
                Ptr{Float32}),
          &nr,kext,ext2bscatt,salb,g,bscatt,&lambd,
          &noMS,&alt,&dr,&theta,&freq)
    #bscatt = sdsu.multiscatterf(kext,ext2bscatt,salb,g,bscatt,lambd,noMS,\
    #                       alt,dr,theta,freq)
    zMS=Array{Float32}(88)*0-99
    #println(bscatt)
    for k=1:88
        if bscatt[k]>0
            zMS[k]=log10(bscatt[k])*10-dZ
        end
    end
    return zMS
end

function getScatt(imelt,wrain,wsnow,wgrpl,whail,icldw,temp,fmelt)
    wcldw=0.
    nrdfreq=6
    kexttot=Array{Float32}(nrdfreq)
    asymtot=Array{Float32}(nrdfreq)
    salbtot=Array{Float32}(nrdfreq)
    kexttotm=Array{Float32}(nrdfreq)
    asymtotm=Array{Float32}(nrdfreq)
    salbtotm=Array{Float32}(nrdfreq)
    asymtotm=Array{Float32}(nrdfreq)
    salbtotm=Array{Float32}(nrdfreq)
    dBZ=Array{Float32}(nrdfreq)
    dBZm=Array{Float32}(nrdfreq)
  
    ccall((:fromlktables2_,"SDSU"),Void, (Ptr{Int32},Ptr{Float32},Ptr{Float32},
                                          Ptr{Float32},Ptr{Float32},Ptr{Float32},
                                          Ptr{Float32},Ptr{Float32},Ptr{Float32},
                                          Ptr{Float32},Ptr{Float32},Ptr{Float32},
                                          Ptr{Float32},Ptr{Float32},Ptr{Float32},
                                       Ptr{Float32},Ptr{Float32},Ptr{Int32}),
           &imelt,&wrain,&wsnow,&wgrpl,&whail,&wcldw,&icldw,
           &temp, kexttot,salbtot,asymtot,kexttotm,
          salbtotm, asymtotm, dBZ, dBZm,&fmelt,&nrdfreq)
    if(wrain+wsnow+wgrpl+whail<0.0000001)
        dBZ=dBZ*0
        kexttot=kexttot*0
    end
    return dBZ,kexttot,dBZm,kexttotm, salbtot, asymtot
          
end


function gauss_newton(zobs,ifreq)
    imelt=0
    wrain=0.1
    wsnow=0.
    wgrpl=0.
    whail=0.
    icldw=0.
    temp=274.
    fmelt=0.
    dn=0.4
    eps=1e3
    for it=1:5
        dBZ,kexttot,dBZm,kexttotm=getScatt(imelt,wrain/dn,wsnow,wgrpl,whail,
                                           icldw,temp,fmelt)
        dBZ=dBZ+log10(dn)*10.
        kexttot=kexttot*dn
        wrain1=1.1*wrain
        dBZ1,kexttot1,dBZm1,kexttotm1=getScatt(imelt,wrain1/dn,wsnow,wgrpl,
                                               whail,icldw,temp,fmelt)
        dBZ1=dBZ1+log10(dn)*10.
        kexttot1=kexttot1*dn
        
        eps=dBZ[ifreq]-zobs
        gradZ=(dBZ1[ifreq]-dBZ[ifreq])/(0.1*wrain)
        dwrain=eps*gradZ/(gradZ*gradZ+0.001)
        wrain=wrain-dwrain
        if wrain<0.01
            wrain=0.01
        end
    end
  
    return wrain, eps
end
