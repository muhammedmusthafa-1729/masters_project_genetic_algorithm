#Optimization for minimum weight for given inlet conditions and outlet
#Optimization parameters are Tube material, tube size, ST/do, SL/do, NT1, dNT, NL-  7 parameters 
#Enthalpy method from Barron, using cross flow local equations - Property variation in both streams with temperature is considered
#Constraints- pressure drop percentage, min bend diameter, pitch angle, min tube thickness etc
from numpy import *
import random
import CoolProp.CoolProp as CP
from math import *
from time import perf_counter
t1=perf_counter()
def fitness(arr):
    #Enthalpy method from Barron, using cross flow local equations
    Fluid="Nitrogen"
    mh=500e-3
    mc=475e-3
    Thi=300
    Tci=80
    Tco=295
    Ph_i=2e6
    Phi=2e6
    Pci=100e3
    I_hi=CP.PropsSI('H','P',Phi,'T',Thi,Fluid)
    I_ci=CP.PropsSI('H','P',Pci,'T',Tci,Fluid)
    I_co=CP.PropsSI('H','P',Pci,'T',Tco,Fluid)
    Q=mc*(I_co-I_ci)
    I_ho=I_hi-Q/mh
    Tho=CP.PropsSI('T','H',I_ho,'P',Phi,Fluid)
    ch=(I_ho-I_hi)/(Tho-Thi)
    cc=(I_co-I_ci)/(Tco-Tci)
    Cc=cc*mc
    Ch=ch*mh
    #Cc mixed, Ch unmixed
    if Cc>Ch:
        Cmin=Ch
        Cmax=Cc
    else:
        Cmin=Cc
        Cmax=Ch
    if Ch==Cmax:
        #Cc is minimum capacity
        I_co_min=CP.PropsSI('H','P',Pci,'T',Thi,Fluid)
        Effectiveness=(I_co-I_ci)/(I_co_min-I_ci)
        #print(Effectiveness)
    CR=Cmin/Cmax
    # Geometry
    x1=arr[0]
    xm1=["Copper","Aluminum"]
    y2=0
    for i in range(1,4):
        y2a=arr[i]*2**(3-i)
        y2=y2+y2a
    x2=y2
    xdo2=[3.175,4.7625,6.35,7.9375,9.525,12.7,15.875,19.05]
    xdi2=[[1.651,3.2512,4.826,6.2992,7.8994,11.0744,14.097,16.9164],[1.397,2.9845,4.572,6.1595,7.747,10.922,13.3858,16.5608]]
    y3=0
    for i in range(4,8):
        y3a=arr[i]
        y3=y3+y3a
    x3=y3
    xv3=[1.125,1.25,1.5,2.0,3.0]
    y4=0
    for i in range(8,12):
        y4a=arr[i]
        y4=y4+y4a
    x4=y4
    xv4=[1.125,1.25,1.5,2.0,3.0]
    y5=0
    for i in range(12,15):
        y5a=arr[i]*2**(14-i)
        y5=y5+y5a
    x5=y5+1
    y6=0
    for i in range(15,18):
        y6a=arr[i]*2**(17-i)
        y6=y6+y6a
    x6=y6+1
    y7=0
    for i in range(18,22):
        y7a=arr[i]*2**(21-i)
        y7=y7+y7a
    x7=y7+1
    Material=xm1[x1]   #1
    if Material=="Copper":
        Kwall=400
        AllowStress=98.6e6/1.5
        RhoTube=8940
        Ewall=117e9 #Young's Modulus
        nuwall=0.33
    elif Material=="Aluminum":
        Kwall=237
        AllowStress=47.6e6/1.5
        RhoTube=2770
        Ewall=73e9 
        nuwall=0.33
    do=xdo2[x2]*1e-3  #2
    di=xdi2[x1][x2]*1e-3
    t=(do-di)/2
    if do/t<15:
        Db=4*do
    else:
        Db=5*do
    minthickness1=Phi*do/(2*(AllowStress+0.4*Phi))
    minthickness2=do*((1-nuwall**2)*Pci/(2*Ewall))**(1/3)
    ST_do=xv3[x3]   #3
    ST=ST_do*do
    aT=ST-do
    SL_do=xv4[x4] #4
    SL=SL_do*do
    NT1=x5          #5
    dNT=x6           #6
    NTn=NT1
    NL=x7              #7
    NT=[]
    for i in range(NL):
        NT.append(NTn)
        NTn=NTn+dNT
    Dpc1=NT1*2*ST/(dNT)
    Dpcn=Dpc1
    Dpc=[]
    for i in range(NL):
        Dpc.append(Dpcn)
        Dpcn=Dpcn+2*ST
    Do=Dpc[NL-1]+ST
    Di=Dpc[0]-ST
    if Di<=0 or Dpc1<Db:
        return 0
    Dm=average(Dpc)
    theta=atan(NT1*SL/(pi*Dpc1))
    if theta*180/pi>20:
        return 0
    T_cmat=[]
    I_cmat=[]
    N=50
    for i in range(N+1):
        Tc=Tco-(i/(N))*(Tco-Tci)
        T_cmat.append(Tc)
    for i in range(len(T_cmat)):
        I_c=CP.PropsSI('H','P',Pci,'T',T_cmat[i],Fluid)
        I_cmat.append(I_c)
    I_hmat=[I_hi]
    for i in range(len(I_cmat)-1):
        I_h=I_hmat[len(I_hmat)-1]-(mc/mh)*(I_cmat[len(I_hmat)-1]-I_cmat[len(I_hmat)])
        I_hmat.append(I_h)
    T_hmat=[]
    for i in range(len(I_hmat)):
        Th=CP.PropsSI('T','P',Phi,'H',I_hmat[i],Fluid)
        T_hmat.append(Th)
    eff=[]
    Area=[]
    Npass=[]
    P_h=[]
    pdropH=[]
    Hcmat=[]

    # Calculating for each of the N heat exchanger segments
    for i in range(N):
        Th_i=T_hmat[i]
        Th_o=T_hmat[i+1]
        Tc_i=T_cmat[i+1]
        Tc_o=T_cmat[i]
        if Phi<1.9e6:
            return 0
        Ic_o=CP.PropsSI('H','P',Pci,'T',Tc_o,Fluid)
        Ic_i=CP.PropsSI('H','P',Pci,'T',Tc_i,Fluid)
        Ih_o=CP.PropsSI('H','P',Phi,'T',Th_o,Fluid)
        Ih_i=CP.PropsSI('H','P',Phi,'T',Th_i,Fluid)
        c_c=(Ic_i-Ic_o)/(Tc_i-Tc_o)
        c_h=(Ih_i-Ih_o)/(Th_i-Th_o)
        C_c=c_c*mc
        C_h=c_h*mh
        #Cc mixed, Ch unmixed
        if C_c>C_h:
            C_min=C_h
            C_max=C_c
        else:
            C_min=C_c
            C_max=C_h
        if C_h==C_max:
            #Cc is minimum capacity
            I_co_min=CP.PropsSI('H','P',Pci,'T',Th_i,Fluid)
            Effectiveness=(Ic_o-Ic_i)/(I_co_min-Ic_i)
        C_R=C_min/C_max
        eff.append(Effectiveness)
        NTU=(1/C_R)*log(abs(1/(1-C_R*log(1/(1-Effectiveness)))))
        UA=NTU*C_min
        #Overall Heat Transfer calculations
        #Calculation of ho
        Dm=Di+NL*ST
        AffC=pi*Dm*NL*aT
        GmaxC=mc/AffC
        Tcavg=(Tc_i+Tc_o)/2
        Thavg=(Th_i+Th_o)/2
        Cpc=CP.PropsSI('C','T',Tcavg,'P',Pci,Fluid)
        MuC=CP.PropsSI('V','T',Tcavg,'P',Pci,Fluid)
        DC=CP.PropsSI('D','T',Tcavg,'P',Pci,Fluid)
        KC=CP.PropsSI('L','T',Tcavg,'P',Pci,Fluid)
        PrC=CP.PropsSI('Prandtl','T',Tcavg,'P',Pci,Fluid)
        ReC=GmaxC*do/MuC
        c=0.386
        n=0.408
        if 1<=SL/do<1.19:
            if 1<=ST/do<1.19:
                c=0.401
                n=0.403
            elif 1.19<=ST/do:
                c=0.377
                n=0.596
        elif 1.19<=SL/do<1.375:
            if 1<=ST/do<1.19:
                c=0.426
                n=0.426
            elif 1.19<=ST/do<1.375:
                c=0.386
                n=0.408
            elif 1.375<=ST/do<1.75:
                c=0.305
                n=0.392
            elif 1.75<=ST/do<2.5:
                c=0.111
                n=0.296
            elif 2.5<=ST/do:
                c=0.0703
                n=0.248
        elif 1.375<=SL/do<1.75:
            if 1<=ST/do<1.19:
                c=0.470
                n=0.470
            elif 1.19<=ST/do<1.375:
                c=0.407
                n=0.414
            elif 1.375<=ST/do<1.75:
                c=0.278
                n=0.380
            elif 1.75<=ST/do<2.5:
                c=0.112
                n=0.298
            elif 2.5<=ST/do:
                c=0.0753
                n=0.256
        elif 1.75<=SL/do<2.5:
            if ST/do<1.375:
                c=0.464
                n=0.430
            elif 1.375<=ST/do<1.75:
                c=0.332
                n=0.398
            elif 1.75<=ST/do<2.5:
                c=0.254
                n=0.368
            elif 2.5<=ST/do:
                c=0.220
                n=0.352
        elif 2.5<=SL/do:
            if ST/do<1.375:
                c=0.332
                n=0.399
            elif 1.375<=ST/do<1.75:
                c=0.396
                n=0.416
            elif 1.75<=ST/do<2.5:
                c=0.415
                n=0.419
            elif 2.5<=ST/do:
                c=0.317
                n=0.392
        JHC=c*ReC**(-n)
        ho=JHC*Cpc*GmaxC/PrC**(2/3)
        #Calculate hi
        NTsum=sum(NT)
        AffH=NTsum*(pi/4)*di**2
        GmaxH=mh/AffH
        Cph=CP.PropsSI('C','T',Thavg,'P',Phi,Fluid)
        MuH=CP.PropsSI('V','T',Thavg,'P',Phi,Fluid)
        DH=CP.PropsSI('D','T',Thavg,'P',Phi,Fluid)
        KH=CP.PropsSI('L','T',Thavg,'P',Phi,Fluid)
        PrH=CP.PropsSI('Prandtl','T',Thavg,'P',Phi,Fluid)
        ReH=GmaxH*di/MuH
        c=1+3.6*(1-di/Dm)*(di/Dm)**0.8
        JHH=c*0.023*ReH**(-0.2)
        hi=JHH*Cph*GmaxH/PrH**(2/3)
        #Calculate Uo and NP
        Uo=(1/ho+do/(di*hi)+do*log(do/di)/(2*Kwall))**(-1)
        Ao=UA/Uo
        Area.append(Ao)
        L=Ao/(pi*do)
        Lt=L/NTsum
        V_tube=L*pi*(do**2-di**2)/4
        M_tube=V_tube*RhoTube
        Nturns1=Lt/(pi*Dpc1/cos(theta))
        N_P=Nturns1*NT1
        Npass.append(N_P)
        L_straight=SL*N_P
        Hcmat.append(L_straight)
        #Pressure drop calculation
        if ReH<2100:
                    fH=16/ReH
        elif ReH<4000:
                fH=0.0054+2.3e-8*ReH**(2/3)
        else:
                fH=0.00128+0.1143*ReH**(-0.311)
        pdrop1=fH*L*GmaxH**2/(2*DH*di)
        p_dropH=pdrop1*(1+0.823*(1+(di/Dm))*(di/Dm)**0.53*(ReH)**0.25)
        P_h.append(Phi)
        if Phi<1.9e6:
            return 0
        Phi=Phi-p_dropH
        pdropH.append(p_dropH)
    Ao=sum(Area) #total heat transfer area
    L=Ao/(pi*do)
    Lt=L/NTsum
    V_tube=L*pi*(do**2-di**2)/4
    M_tube=V_tube*RhoTube
    Nturns1=Lt/(pi*Dpc1/cos(theta))
    NP=Nturns1*NT1
    L_straight=SL*NP
    Hc=sum(Hcmat)
    p_dropH=sum(pdropH)
    #Calculation of cold side pressure drop, average
    Dm=Di+NL*ST
    AffC=pi*Dm*NL*aT
    GmaxC=mc/AffC
    Tcavg=(Tci+Tco)/2
    Thavg=(Thi+Tho)/2
    Tw=(Tcavg+Thavg)/2
    Tcfavg=(Tcavg+Tw)/2
    Cpc=CP.PropsSI('C','T',Tcfavg,'P',Pci,Fluid)
    MuC=CP.PropsSI('V','T',Tcfavg,'P',Pci,Fluid)
    DC=CP.PropsSI('D','T',Tcfavg,'P',Pci,Fluid)
    KC=CP.PropsSI('L','T',Tcfavg,'P',Pci,Fluid)
    PrC=CP.PropsSI('Prandtl','T',Tcfavg,'P',Pci,Fluid)
    ReC=GmaxC*do/MuC
    if ReC<100:
        fC=0.435*ReC**(-0.133)*(sin(theta))**(-0.36)
        p_dropC=fC*L_straight*GmaxC**2/(2*DC*do)
    else:
        n=0.43+1.13*do/SL
        fC=(0.176+0.32*(ST/do)*((SL/do)-1)**(-n))*ReC**(-0.15)
        p_dropC=fC*NL*GmaxC**2/(2*9.81*DC)
    percentpdropC=p_dropC/Pci
    percentpdropH=p_dropH/Ph_i
    if percentpdropC>0.05 or percentpdropH>0.05:
        return 0
    fit=1/M_tube
    return fit
def pitchfit(arr):
    x1=arr[0]
    xm1=["Copper","Aluminum"]
    y2=0
    for i in range(1,4):
        y2a=arr[i]*2**(3-i)
        y2=y2+y2a
    x2=y2
    xdo2=[3.175,4.7625,6.35,7.9375,9.525,12.7,15.875,19.05]
    xdi2=[[1.651,3.2512,4.826,6.2992,7.8994,11.0744,14.097,16.9164],[1.397,2.9845,4.572,6.1595,7.747,10.922,13.3858,16.5608]]
    y3=0
    for i in range(4,8):
        y3a=arr[i]
        y3=y3+y3a
    x3=y3
    xv3=[1.125,1.25,1.5,2.0,3.0]
    y4=0
    for i in range(8,12):
        y4a=arr[i]
        y4=y4+y4a
    x4=y4
    xv4=[1.125,1.25,1.5,2.0,3.0]
    y5=0
    for i in range(12,15):
        y5a=arr[i]*2**(14-i)
        y5=y5+y5a
    x5=y5+1
    y6=0
    for i in range(15,18):
        y6a=arr[i]*2**(17-i)
        y6=y6+y6a
    x6=y6+1
    y7=0
    for i in range(18,22):
        y7a=arr[i]*2**(21-i)
        y7=y7+y7a
    x7=y7+1
    do=xdo2[x2]*1e-3  #2
    di=xdi2[x1][x2]*1e-3
    t=(do-di)/2
    if do/t<15:
        Db=4*do
    else:
        Db=5*do
    ST_do=xv3[x3]   #3
    ST=ST_do*do
    aT=ST-do
    SL_do=xv4[x4] #4
    SL=SL_do*do
    NT1=x5          #5
    dNT=x6           #6
    NTn=NT1
    NL=x7              #7
    Dpc1=NT1*2*ST/(dNT)
    theta=atan(NT1*SL/(pi*Dpc1))
    if Dpc1<Db or theta*180/pi>20:
        return 0
    else:
        return 1    

popsize=30
matingpoolsize=30
numelites=2
mutationrate=0.10
generations=1000

elite=[]        #elite
pop=[]
write=open("VOID.txt",'w')
for h in range(popsize):
    c=[]
    for i in range(22):
        a=random.randint(0,1)
        c.append(a)
    pop.append(c)
print("initial population",pop)
for x in range(generations):
    weightage=[]
    fitnessmax=0
    for i in range(len(pop)):
        d=fitness(pop[i])
        if d>=fitnessmax:    #elite
            hu=pop[i]       #elite
            fitnessmax=d    #elite
        weightage.append(round(d,5))
    #extra
    if fitnessmax==0:
        for h in range(popsize):
            c=[]
            for i in range(22):
                a=random.randint(0,1)
                c.append(a)
            pop.append(c)
        print("initial population",pop)
        continue
    #extra
    mating_pool=random.choices(pop,k=matingpoolsize,weights=weightage)
    popweight=dict(zip(weightage,pop))
    weightage.sort()
    weightage.reverse()
    sortedpop=[]
    for i in weightage:
        sortedpop.append(popweight[i]) 
    pop=[]
    #Adding elite chromosomes
    for i in range(numelites):
        pop.append(sortedpop[i]) 
    #crossover
    while len(pop)<=popsize:
        a1=random.randint(0,popsize-1)
        a2=random.randint(0,popsize-1)
        str1=mating_pool[a1]
        str2=mating_pool[a2]
        pos=random.randint(1,len(str1))
        child1=str1[pos:]+str2[:pos]
        child2=str2[pos:]+str1[:pos]
        if pitchfit(child1)!=0:
            pop.append(child1)
        if pitchfit(child2)!=0:
            pop.append(child2)
    for i in range(numelites,popsize):
        for j in range(22):
            n=random.random()
            if n<mutationrate:
                if pop[i][j]==0:
                    pop[i][j]=1
                else:
                    pop[i][j]=0
    wr1=round(1/max(weightage),5)
    if average(weightage)!=0:
        print(max(weightage),"Min weight = ",wr1,x+1,pop[0])
    wr2=str(wr1)
    write.write(wr2)
    write.write("\n")
    print(weightage)
    print(x+1)
def maxfitness(arr):
    #Enthalpy method from Barron, using cross flow local equations
    Fluid="Nitrogen"
    mh=500e-3
    mc=475e-3
    Thi=300
    Tci=80
    Tco=295
    Phi=2e6
    Ph_i=2e6
    Pci=100e3
    I_hi=CP.PropsSI('H','P',Phi,'T',Thi,Fluid)
    I_ci=CP.PropsSI('H','P',Pci,'T',Tci,Fluid)
    I_co=CP.PropsSI('H','P',Pci,'T',Tco,Fluid)
    Q=mc*(I_co-I_ci)
    print("Q = ",Q)
    I_ho=I_hi-Q/mh
    Tho=CP.PropsSI('T','H',I_ho,'P',Phi,Fluid)
    ch=(I_ho-I_hi)/(Tho-Thi)
    cc=(I_co-I_ci)/(Tco-Tci)
    Cc=cc*mc
    Ch=ch*mh
    #Cc mixed, Ch unmixed
    if Cc>Ch:
        Cmin=Ch
        Cmax=Cc
    else:
        Cmin=Cc
        Cmax=Ch
    if Ch==Cmax:
        #Cc is minimum capacity
        I_co_min=CP.PropsSI('H','P',Pci,'T',Thi,Fluid)
        Effectiveness=(I_co-I_ci)/(I_co_min-I_ci)
        #print(Effectiveness)
        Overalleffectiveness=Effectiveness
    # Geometry
    x1=arr[0]
    xm1=["Copper","Aluminum"]
    y2=0
    for i in range(1,4):
        y2a=arr[i]*2**(3-i)
        y2=y2+y2a
    x2=y2
    xdo2=[3.175,4.7625,6.35,7.9375,9.525,12.7,15.875,19.05]
    xdi2=[[1.651,3.2512,4.826,6.2992,7.8994,11.0744,14.097,16.9164],[1.397,2.9845,4.572,6.1595,7.747,10.922,13.3858,16.5608]]
    y3=0
    for i in range(4,8):
        y3a=arr[i]
        y3=y3+y3a
    x3=y3
    xv3=[1.125,1.25,1.5,2.0,3.0]
    y4=0
    for i in range(8,12):
        y4a=arr[i]
        y4=y4+y4a
    x4=y4
    xv4=[1.125,1.25,1.5,2.0,3.0]
    y5=0
    for i in range(12,15):
        y5a=arr[i]*2**(14-i)
        y5=y5+y5a
    x5=y5+1
    y6=0
    for i in range(15,18):
        y6a=arr[i]*2**(17-i)
        y6=y6+y6a
    x6=y6+1
    y7=0
    for i in range(18,22):
        y7a=arr[i]*2**(21-i)
        y7=y7+y7a
    x7=y7+1
    Material=xm1[x1]   #1
    if Material=="Copper":
        Kwall=400
        AllowStress=98.6e6/1.5
        RhoTube=8940
        Ewall=117e9 #Young's Modulus
        nuwall=0.33
    elif Material=="Aluminum":
        Kwall=237
        AllowStress=47.6e6/1.5
        RhoTube=2770
        Ewall=73e9 
        nuwall=0.33
    do=xdo2[x2]*1e-3  #2
    di=xdi2[x1][x2]*1e-3
    t=(do-di)/2
    if do/t<15:
        Db=4*do
    else:
        Db=5*do
    minthickness1=Phi*do/(2*(AllowStress+0.4*Phi))
    minthickness2=do*((1-nuwall**2)*Pci/(2*Ewall))**(1/3)
    ST_do=xv3[x3]   #3
    ST=ST_do*do
    aT=ST-do
    SL_do=xv4[x4] #4
    SL=SL_do*do
    NT1=x5          #5
    dNT=x6           #6
    NTn=NT1
    NL=x7              #7
    NT=[]
    for i in range(NL):
        NT.append(NTn)
        NTn=NTn+dNT
    Dpc1=NT1*2*ST/(dNT)
    Dpcn=Dpc1
    Dpc=[]
    for i in range(NL):
        Dpc.append(Dpcn)
        Dpcn=Dpcn+2*ST
    Do=Dpc[NL-1]+ST
    Di=Dpc[0]-ST
    if Di<=0 or Dpc1<Db:
        return 0
    Dm=average(Dpc)
    theta=atan(NT1*SL/(pi*Dpc1))
    if theta*180/pi>20:
        return 0
    T_cmat=[]
    I_cmat=[]
    N=50
    for i in range(N+1):
        Tc=Tco-(i/(N))*(Tco-Tci)
        T_cmat.append(Tc)
    for i in range(len(T_cmat)):
        I_c=CP.PropsSI('H','P',Pci,'T',T_cmat[i],Fluid)
        I_cmat.append(I_c)
    I_hmat=[I_hi]
    for i in range(len(I_cmat)-1):
        I_h=I_hmat[len(I_hmat)-1]-(mc/mh)*(I_cmat[len(I_hmat)-1]-I_cmat[len(I_hmat)])
        I_hmat.append(I_h)
    T_hmat=[]
    for i in range(len(I_hmat)):
        Th=CP.PropsSI('T','P',Phi,'H',I_hmat[i],Fluid)
        T_hmat.append(Th)
    eff=[]
    Area=[]
    Npass=[]
    P_h=[]
    pdropH=[]
    NTUmat=[]
    Hcmat=[]
    for i in range(N):
        Th_i=T_hmat[i]
        Th_o=T_hmat[i+1]
        Tc_i=T_cmat[i+1]
        Tc_o=T_cmat[i]
        if Phi<1.9e6:
            return 0
        Ic_o=CP.PropsSI('H','P',Pci,'T',Tc_o,Fluid)
        Ic_i=CP.PropsSI('H','P',Pci,'T',Tc_i,Fluid)
        Ih_o=CP.PropsSI('H','P',Phi,'T',Th_o,Fluid)
        Ih_i=CP.PropsSI('H','P',Phi,'T',Th_i,Fluid)
        c_c=(Ic_i-Ic_o)/(Tc_i-Tc_o)
        c_h=(Ih_i-Ih_o)/(Th_i-Th_o)
        C_c=c_c*mc
        C_h=c_h*mh
        #Cc mixed, Ch unmixed
        if C_c>C_h:
            C_min=C_h
            C_max=C_c
        else:
            C_min=C_c
            C_max=C_h
        if C_h==C_max:
            #Cc is minimum capacity
            I_co_min=CP.PropsSI('H','P',Pci,'T',Th_i,Fluid)
            Effectiveness=(Ic_o-Ic_i)/(I_co_min-Ic_i)
        C_R=C_min/C_max
        eff.append(Effectiveness)
        NTU=(1/C_R)*log(abs(1/(1-C_R*log(1/(1-Effectiveness)))))
        NTUmat.append(NTU)
        UA=NTU*C_min
        #Overall Heat Transfer calculations
        #Calculation of ho
        Dm=Di+NL*ST
        AffC=pi*Dm*NL*aT
        GmaxC=mc/AffC
        Tcavg=(Tc_i+Tc_o)/2
        Thavg=(Th_i+Th_o)/2
        Cpc=CP.PropsSI('C','T',Tcavg,'P',Pci,Fluid)
        MuC=CP.PropsSI('V','T',Tcavg,'P',Pci,Fluid)
        DC=CP.PropsSI('D','T',Tcavg,'P',Pci,Fluid)
        KC=CP.PropsSI('L','T',Tcavg,'P',Pci,Fluid)
        PrC=CP.PropsSI('Prandtl','T',Tcavg,'P',Pci,Fluid)
        ReC=GmaxC*do/MuC
        c=0.386
        n=0.408
        if 1<=SL/do<1.19:
            if 1<=ST/do<1.19:
                c=0.401
                n=0.403
            elif 1.19<=ST/do:
                c=0.377
                n=0.596
        elif 1.19<=SL/do<1.375:
            if 1<=ST/do<1.19:
                c=0.426
                n=0.426
            elif 1.19<=ST/do<1.375:
                c=0.386
                n=0.408
            elif 1.375<=ST/do<1.75:
                c=0.305
                n=0.392
            elif 1.75<=ST/do<2.5:
                c=0.111
                n=0.296
            elif 2.5<=ST/do:
                c=0.0703
                n=0.248
        elif 1.375<=SL/do<1.75:
            if 1<=ST/do<1.19:
                c=0.470
                n=0.470
            elif 1.19<=ST/do<1.375:
                c=0.407
                n=0.414
            elif 1.375<=ST/do<1.75:
                c=0.278
                n=0.380
            elif 1.75<=ST/do<2.5:
                c=0.112
                n=0.298
            elif 2.5<=ST/do:
                c=0.0753
                n=0.256
        elif 1.75<=SL/do<2.5:
            if ST/do<1.375:
                c=0.464
                n=0.430
            elif 1.375<=ST/do<1.75:
                c=0.332
                n=0.398
            elif 1.75<=ST/do<2.5:
                c=0.254
                n=0.368
            elif 2.5<=ST/do:
                c=0.220
                n=0.352
        elif 2.5<=SL/do:
            if ST/do<1.375:
                c=0.332
                n=0.399
            elif 1.375<=ST/do<1.75:
                c=0.396
                n=0.416
            elif 1.75<=ST/do<2.5:
                c=0.415
                n=0.419
            elif 2.5<=ST/do:
                c=0.317
                n=0.392
        JHC=c*ReC**(-n)
        ho=JHC*Cpc*GmaxC/PrC**(2/3)
        #Calculate hi
        NTsum=sum(NT)
        AffH=NTsum*(pi/4)*di**2
        GmaxH=mh/AffH
        Cph=CP.PropsSI('C','T',Thavg,'P',Phi,Fluid)
        MuH=CP.PropsSI('V','T',Thavg,'P',Phi,Fluid)
        DH=CP.PropsSI('D','T',Thavg,'P',Phi,Fluid)
        KH=CP.PropsSI('L','T',Thavg,'P',Phi,Fluid)
        PrH=CP.PropsSI('Prandtl','T',Thavg,'P',Phi,Fluid)
        ReH=GmaxH*di/MuH
        c=1+3.6*(1-di/Dm)*(di/Dm)**0.8
        JHH=c*0.023*ReH**(-0.2)
        hi=JHH*Cph*GmaxH/PrH**(2/3)
        #Calculate Uo and NP
        Uo=(1/ho+do/(di*hi)+do*log(do/di)/(2*Kwall))**(-1)
        Ao=UA/Uo
        Area.append(Ao)
        L=Ao/(pi*do)
        Lt=L/NTsum
        V_tube=L*pi*(do**2-di**2)/4
        M_tube=V_tube*RhoTube
        Nturns1=Lt/(pi*Dpc1/cos(theta))
        N_P=Nturns1*NT1
        Npass.append(N_P)
        L_straight=SL*N_P
        Hcmat.append(L_straight)
        #Pressure drop calculation
        if ReH<2100:
                    fH=16/ReH
        elif ReH<4000:
                fH=0.0054+2.3e-8*ReH**(2/3)
        else:
                fH=0.00128+0.1143*ReH**(-0.311)
        pdrop1=fH*L*GmaxH**2/(2*DH*di)
        p_dropH=pdrop1*(1+0.823*(1+(di/Dm))*(di/Dm)**0.53*(ReH)**0.25)
        P_h.append(Phi)
        if Phi<1.9e6:
            return 0
        Phi=Phi-p_dropH
        pdropH.append(p_dropH)
    Ao=sum(Area)
    L=Ao/(pi*do)
    Lt=L/NTsum
    V_tube=L*pi*(do**2-di**2)/4
    M_tube=V_tube*RhoTube
    Nturns1=Lt/(pi*Dpc1/cos(theta))
    NP=Nturns1*NT1
    L_straight=SL*NP
    Hc=sum(Hcmat)
    p_dropH=sum(pdropH)
    #Calculation of cold side pressure drop, average
    Dm=Di+NL*ST
    AffC=pi*Dm*NL*aT
    GmaxC=mc/AffC
    Tcavg=(Tci+Tco)/2
    Thavg=(Thi+Tho)/2
    Tw=(Tcavg+Thavg)/2
    Tcfavg=(Tcavg+Tw)/2
    Cpc=CP.PropsSI('C','T',Tcfavg,'P',Pci,Fluid)
    MuC=CP.PropsSI('V','T',Tcfavg,'P',Pci,Fluid)
    DC=CP.PropsSI('D','T',Tcfavg,'P',Pci,Fluid)
    KC=CP.PropsSI('L','T',Tcfavg,'P',Pci,Fluid)
    PrC=CP.PropsSI('Prandtl','T',Tcfavg,'P',Pci,Fluid)
    ReC=GmaxC*do/MuC
    if ReC<100:
        fC=0.435*ReC**(-0.133)*(sin(theta))**(-0.36)
        p_dropC=fC*L_straight*GmaxC**2/(2*DC*do)
    else:
        n=0.43+1.13*do/SL
        fC=(0.176+0.32*(ST/do)*((SL/do)-1)**(-n))*ReC**(-0.15)
        p_dropC=fC*NL*GmaxC**2/(2*9.81*DC)
    percentpdropC=p_dropC/Pci
    percentpdropH=p_dropH/Ph_i
    if percentpdropC>0.05 or percentpdropH>0.05:
        return 0
    V0=pi*Do**2*L_straight/4
    Compactness=Ao/V0
    fit=1/M_tube
    print("Material, Tube size, ST/do, SL/do, NT1 =", Material, do,ST_do,SL_do,NT1)
    print("Increment in number of tubes with each layer =", dNT)
    print("Total number of Tubes NTsum = ", NTsum)
    print("Number of Layers", NL)
    print("Number of Passes", NP)
    print("Heat exchanger length = ", L_straight)
    print("Total and single tube length =",L, Lt)
    print("Mandrel diameter, Shell diameter =", Di, Do)
    print("Effectiveness = ", Overalleffectiveness)
    print("Compactness = ", Compactness)
    print("theta =", theta*180/pi)
    print("Mass of heat exchanger = ", M_tube)
    print("Pressure drop in shell side, tube side = ",percentpdropC*100 ,"(",p_dropC,")", percentpdropH*100,"(",p_dropH,")")
    print("Volume of tubing = ", V_tube)
    print(Thi,Tho,Tci,Tco)
    print("NTU =",sum(NTUmat))
    print("mc, mh ", mc,mh)
    return(fit)
maxfitness(pop[0])
print(pop[0])
t2=perf_counter()
print("time taken = ", t2-t1)
write.write("time")
write.write(str(t2-t1))
write.close()
#print(weightagewithtime)