# -*- coding: utf-8 -*-
##################
def epsilon(u):
    return 0.5*(nabla_grad(u) + nabla_grad(u).T)
def sigmaL(u):
    return lambda_*div(u)*Identity(dim) + Fmuequiv*epsilon(u)
ac =0.7;MGrande=50*mu0;MGrande2=(MGrande/(1e-1*mu0))**2;
def eta(alpha) :
    return eta1 + eta2 / (ac - alpha) 
def deta(alpha) :
    return eta2 / (ac - alpha)**2
def etaM(alpha) :
    return MGrande*((alpha-ac)   +pow((alpha-ac)**2+1/MGrande2,0.5))
def detaM(alpha) :
    return MGrande*(1+ (alpha-ac)/pow((alpha-ac)**2+1/MGrande2,0.5))
def I1(u) :
    return div(u)
def I2(u) :
    return inner(epsilon (u),epsilon(u))
def Xi(u) :
    return I1(u) /( sqrt(I2(u)) + 1e-9) 
def Max(a, b) : 
    return (a+b+abs(a-b))/2
def Min(a, b) : 
    return (a+b-abs(a-b))/2
    
    
aparabola=1e2;dac=0.1;    
def amas(alpha) :
    return Max(0,alpha-ac+dac)     
def eta(alpha) :
    return aparabola*amas(alpha)**2 
def deta(alpha) :
    return 2*aparabola*amas(alpha)
##################    
    
    
    
def RHSalpha(u_nl) :
    RHSalpha=  mur * I2(u_nl) 
    if (VersionDivEta==1):
        RHSalpha = (RHSalpha+ gammar*    div(u_nl)   * sqrt(I2(u_nl)))*Max(Xi(u_nl)-xi0,0)/(abs(Xi(u_nl) - xi0)+1e-10 )
    elif (VersionDivEta==2) :
        RHSalpha = RHSalpha+ gammar *abs(div(u_nl))  * sqrt(I2(u_nl));
    elif (VersionDivEta==3) :
        RHSalpha = RHSalpha+ gammar *Max(div(u_nl),0)* sqrt(I2(u_nl));
    
    
    tmp=project(RHSalpha,Va)
    print "RHSalpha_min=",tmp.vector().min(),"max=",tmp.vector().max()
    return RHSalpha
    

def eta_deta(alpha_nl):
    global etaalpha,detaalpha,VersionEta
    if (VersionEta==1) :
        etaalpha=eta(alpha_nl);
        detaalpha=deta(alpha_nl);
    elif (VersionEta==2) :
        etaalpha=etaM(alpha_nl);
        detaalpha=detaM(alpha_nl);
    return
    

def cxi1(u_nl,alpha_nl):
    global gammar
    lcxi1 =   gammar  * sqrt( I2(u_nl) ) * alpha_nl; 
    return lcxi1
    

def muequiv(u_nl,alpha_nl):
    if (UsoEpsilonUpunto==1) :
        lmuequiv= 2.0*mu0  - alpha_nl*(2*mur + gammar*Xi(u_nl) ) + 1/dt * etaalpha
    else :
        lmuequiv= 2.0*mu0  - alpha_nl*(2*mur + gammar*Xi(u_nl) ) +  etaalpha             
    return lmuequiv

