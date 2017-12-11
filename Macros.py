# -*- coding: utf-8 -*-
##################
def epsilon(u):
    return 0.5*(nabla_grad(u) + nabla_grad(u).T)
def sigmaL(u):
    return lambda_*div(u)*Identity(dim) + Fmuequiv*epsilon(u)
ac =0.71;MGrande=500*mu0;MGrande2=(MGrande/(1e-1*mu0))**2;
def eta(alpha) :
    return eta1 + eta2 / (ac - alpha) 
def deta(alpha) :
    return eta2 / (ac - alpha)**2
def etaM(alpha) :
    return MGrande*((alpha-ac)   +sqrt((alpha-ac)**2+1/MGrande2))
def detaM(alpha) :
    return MGrande*(1+ (alpha-ac)/sqrt((alpha-ac)**2+1/MGrande2))
def I1(u) :
    return div(u)
def I2(u) :
    return inner(epsilon (u),epsilon(u))
def Xi(u) :
    return I1(u) /( sqrt(I2(u)) + 1e-9) 
def Max(a, b): return (a+b+abs(a-b))/Constant(2)
##################    
