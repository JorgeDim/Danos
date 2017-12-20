# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 20:48:30 2017

@author: user1
"""

kappaD_exp=Expression('kappaD',degree=1,kappaD=kappaD)
kappaD2_exp=Expression('factor',degree=1,factor=100)


vd,w = TestFunctions(Vua)
ud, alpha=TrialFunctions(Vua)
dim = ud.geometric_dimension()  # space dimension

def um(ud) :
    return  up+(dt/4)*(udp+ud)

######## Lado IZQUIERDO #####
FV        =     (rho/dt)*inner(ud,vd)*dx 
FV        +=    inner(sigmaL(um(ud)), epsilon(vd))*dx  
if (UsoEpsilonUpunto==1) :
    FV    +=    Fetaalpha*inner(epsilon(   ud   ), epsilon(vd))*dx
else :
    FV    +=    Fetaalpha*inner(epsilon(um(ud)), epsilon(vd))*dx
FV        +=   (1/C/dt)*alpha*w *dx 
FV        +=    inner(kappaD_exp*grad(alpha),grad(w))*dx 

###### Lado DERECHO ###############
FV        += - dot(FuerzaDot, vd)*dx \
             - Fcxi1 * div(vd) *dx
if UsoDalphaDt!=0:
    FV    += - FRHSalpha *w*dx 
###### Separacion ################    
Bilineal,Lineal = lhs(FV), rhs(FV)

#####################Modelo Estacionario#########################
FactorGrande = Function(Va)
FactorGrande.assign(interpolate(Constant(2000),Va))
v,w = TestFunctions(Vua)
u, alpha=TrialFunctions(Vua)
FVE =           inner(sigmaL(u), epsilon(v))*dx 
FVE        +=   (FactorGrande+FnEtaPrima)*alpha*w*dx 
 
if (UsoEpsilonUpunto==0) :
    FVE    +=   Fetaalpha*inner(epsilon(u), epsilon(v))*dx
FVE        +=   inner(kappaD_exp*grad(alpha),grad(w))*dx 
###### Lado DERECHO ###############
FVE        += - dot(Fuerza, v)*dx \
              - Fcxi1 * div(v) *dx
if UsoDalphaDt!=0:
    FVE    += - FRHSalpha *w*dx 
###### Separacion ################    
BilinealE,LinealE = lhs(FVE), rhs(FVE)
##############################################
        
####################Problema Estacionario Alpha###################
w = TestFunction(Va)
alpha=TrialFunction(Va)

FVEA        =   inner(kappaD_exp*grad(alpha),grad(w))*dx 
###### Lado DERECHO ###############
FVEA       += - FRHSalpha *w*dx 
###### Separacion ################    
BilinealEA,LinealEA = lhs(FVEA), rhs(FVEA)
##############################################

#####################Modelo Velocidad (Alpha==dato) #########################
FactorGrande = Function(Va)
FactorGrande.assign(interpolate(Constant(2000),Va))
v = TestFunction(Vu)
u =TrialFunction(Vu)
FV =           inner(sigmaL(u), epsilon(v))*dx 
 
if (UsoEpsilonUpunto==0) :
    FV    +=   Fetaalpha*inner(epsilon(u), epsilon(v))*dx
###### Lado DERECHO ###############
FV        += - dot(Fuerza, v)*dx \
              - Fcxi1 * div(v) *dx
###### Separacion ################    
BilinealEV,LinealEV = lhs(FV), rhs(FV)
##############################################