from numpy import *

def hb(zw,h1,dr):
    betaS=1.13
    betaM=1.13
    betaR=1.12
    zetaS=0
    alphaS=10**(-2.005)
    alphaM=10**(-1.87)
    alphaR=10**(-1.77)
    zwc=zw.copy()
    zeta=zw.copy()*0.
    betav=h1.copy()*0.+betaS
    zmax=-30
    beta=betaM
    zs=-99
    b=1.0269  
    a=10**(-0.890047)
    a=10**(-0.60585)
    b=1.09

    for k in range(h1.shape[0]):
        if(zw[k]>-15):
            if h1[k]>2 and zw[k]>zmax:
                zmax=zw[k]
            if h1[k]>2.5:
                q=0.2*log(10.)*betaS
                zeta[k]=zetaS+q*alphaS*10.**(0.1*zw[k]*betaS)*dr
                zetaS=zeta[k]
                betav[k]=betaS
            else:
                if h1[k]>2:
                    q=0.2*log(10.)*betaS
                    zeta[k]=zetaS+q*alphaS*10.**(0.1*zw[k]*betaM)*dr
                    zetaS=zeta[k]
                    betav[k]=betaM
                else:
                    if h1[k]>0.48:
                        q=0.2*log(10.)*betaR
                        zeta[k]=zetaS+q*alphaR*10.**(0.1*zw[k]*betaR)*dr
                        zetaS=zeta[k]
                        betav[k]=betaR
            beta=betav[k]
    if(zetaS>0.98):
        eps=0.98/zetaS
    else:
        eps=1.
    piaW=-log10(1-eps*zetaS)*10/beta
    #print -log10(1-eps*zetaS)*10/beta, beta, piaW
    for it in range(2):
        iset=0
        for k in range(h1.shape[0]):
            if h1[k]>2 and zwc[k]>zmax:
                zmax=zwc[k]
            if(zwc[k]>-15):
                zwc[k]=zw[k]-log10(1-eps*zeta[k])*10/betav[k]
            if h1[k]>0.48+0.24:
                zsfc=zwc[k]
                zsfc_m=zw[k]
                iset=1
        print zsfc, zmax
        if piaW>4 and iset==1 and zsfc<zmax-1:
            eps=(1-10**(0.1*(zsfc_m-(zmax-1))*beta))/zetaS

        piaAv=0
        ic=0
        dpiaS=-99
        sfcPwc=0
        for k in range(h1.shape[0]):
            if(zwc[k]>-15):
                if(h1[k]>0.48+0.24):
                    zwc[k]=zw[k]-log10(1-eps*zeta[k])*10/betav[k]
                    if h1[k]<2.0:
                        sfcPwc=10.**(-2.5271+0.1*1.03789*zwc[k])
                        sfcPwc=10.**(-2.4502+0.1*1.17625*zwc[k])
                        if(dpiaS<0):
                            dpiaS=zwc[k]-zw[k]
                        else:
                            piaAv=zwc[k]-zw[k]-dpiaS
                            ic+=1
                else:
                    zwc[k]=-99
    
        if ic>1:
            piaAv=piaAv/(dr*2)/ic
            sfcPwc2=a*piaAv**b
        else:
            sfcPwc2=0.
        
    piaW=-log10(1-eps*zetaS)*10/beta
    print piaAv,piaW,dpiaS
    return zwc,piaW,sfcPwc,sfcPwc2
            
