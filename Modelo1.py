# -*- coding: utf-8 -*-
# JSM: Tipos de modelo a testear:

from __future__ import print_function
from fenics import *
from mshr import *
import numpy as np
import matplotlib.pyplot as plt


VersionEta=2              #0: no hay eta, 
                        #1: Se usa la función eta(alpha) del articulo. 
                        #2: Se usa la función eta(alpha) propuesta por JSM
UsoEpsilonUpunto=0 #1: Se usa el modelo con el término \eta(\alpha)*E(u):E(udot). 
                        #0: Se usa E(u) en lugar de E(udot) en dicha parte del modelo "propuesto por JSM"
VersionDivEta=2    #1: Se usa termino: RHSalpha+= gammar *    div(u1nl,u2nl,u3nl)             * sqrt( I2(u1nl,u2nl,u3nl) );
                        #2: Se usa término: RHSalpha+= gammar *    abs(div(u1nl,u2nl,u3nl))        * sqrt( I2(u1nl,u2nl,u3nl) );
                        #3: Se usa término: RHSalpha+= gammar *    max(div(u1nl,u2nl,u3nl),0*u1nl) * sqrt( I2(u1nl,u2nl,u3nl) )
UsoDalphaDt=0;                         
#runfile("Defaults.py")
prefix=	""+`VersionEta`+`UsoEpsilonUpunto`+`VersionDivEta`+`UsoDalphaDt`+"_";			

print('prefix=',prefix)

#-----------------------------------------
# Model parameters.
#-----------------------------------------
MP = 1e6 ;
mu0    = 30000 * MP ;
lambda_ = 30000 * MP;
rho    = 2700;
E  = mu0 *( 3 * lambda_ + 2 * mu0) / (lambda_ + mu0); 
nu = lambda_ /  (2 * ( lambda_ + mu0)); 
tEnd = 10;
Cgravedad = 9.8e0;
xi0 = -0.008;     


gammar = 35000 * MP; 
mur =  1.8 * gammar;
C = 1e-6; 
kappa = 1e-9;
eta1 = 1e-2;
eta2 = 1e-1;


dt = 1e-3 ;
x0 = -100;
y0 = -100;
z0 = -100;
x1 =  100;
y1 =  100;
z1 =  100;
sigmaP=10*MP/0.05;

mesh = BoxMesh(Point(x0, y0, z0), Point(x1, y1, z1), 10, 10, 10)
P1v = VectorElement('P', tetrahedron, 1)
P1e = FiniteElement('P', tetrahedron, 1)
element =MixedElement([P1v, P1e])
Vu = FunctionSpace(mesh, P1v)
Vua = FunctionSpace(mesh, element)
Va= FunctionSpace(mesh, P1e)

alphanl  = Function(Va)
Fmuequiv  = Function(Va)
Fetaalpha = Function(Va)
Fcxi1     = Function(Va)
FRHSalpha = Function(Va)
un       = Function(Vu)
unl      = Function(Vu)
Fuerza      = Function(Vu)


##################
def epsilon(u):
    return 0.5*(nabla_grad(u) + nabla_grad(u).T)
def sigma0(u):
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
tol = 1E-1
def lados(x, on_boundary):
    return on_boundary and (near(x[0],x0,tol) or near(x[0],x1,tol))
def fondo(x, on_boundary):
    return on_boundary and (near(x[2],z0,tol) )
def frentefondo(x, on_boundary):
    return on_boundary and (near(x[1],y0,tol) or near(x[1],y1,tol))

bc1 = DirichletBC(Vua.sub(0).sub(0), Constant((0)), lados)
bc2 = DirichletBC(Vua.sub(0).sub(2), Constant((0)), fondo)
bc3 = DirichletBC(Vua.sub(0).sub(1), Constant((0)), frentefondo)
bc=[bc1,bc2,bc3]


ZeroVector = Expression(('0','0','0'),degree=1)
ZeroEscalar= Expression('0',degree=1)


v,w = TestFunctions(Vua)
u, alpha=TrialFunctions(Vua)
dim = u.geometric_dimension()  # space dimension

T = Constant((0, 0, 0))
Bilineal = (rho/dt/dt)*inner(u,v)*dx + inner(sigma0(u), epsilon(v))*dx +   \
    +(1/C/dt)*alpha*w *dx + inner(kappa*grad(alpha),grad(w))*dx 
L1 = dot(Fuerza, v)*dx + dot(T, v)*ds 
L2 = Fcxi1 * nabla_div(v) *dx
L3=0
if UsoEpsilonUpunto != 0:
    L3 = (1/dt*Fetaalpha)*inner(epsilon(un),epsilon(v))*dx 
L4=0
if UsoDalphaDt!=0:
    L4 = FRHSalpha *w*dx 
Lineal=L1+L2+L3+L4

        

# Initial condition


uo = interpolate(ZeroVector, Vu)
un.assign(project(uo + dt * Constant((0, 0, 0)), Vu))
alphan = project(ZeroEscalar, Va)

time = 0
k = 1
etaalpha=0;detaalpha=0
nt = int(tEnd/dt)

zmaxdet=[0];
alphamax=[0];
alphamin=[0];
tiempos=[0];
#%%  Iteraciones Temporales
for i in range(0,nt):
    print("\n=================\nIteración ",i,"\n=================")
    
    time = time + dt; 
    print(" time = ",time," dt = ",dt)
    
    # Non-linear terms are implemented by fixed-point iterations (u1nl...alphanl)
    
    fv = Constant((0, 0, -rho*Cgravedad)) + (2*un-uo)*rho /dt/dt
    fv=project(fv,Vu)
    Fuerza.assign(fv)  ## Esta es la funcion usada en la FV

    u_1_, u_2_, u_3_ = Fuerza.split(deepcopy=True)
    print("f3_min=",u_3_.vector().min(),"max=",u_3_.vector().max())
    
    unl.assign(un)
    alphanl.assign(alphan)
    errnl=1;errnl2=10;relaj=1;itnl=0
    while ((errnl>1e-8 or itnl<2) and itnl<50 ) :
        itnl+=1
        print("itnl = ",itnl) 
        if (VersionEta==1) :
            etaalpha=eta(alphanl);
            detaalpha=deta(alphanl);
        elif (VersionEta==2) :
            etaalpha=etaM(alphanl);
            detaalpha=detaM(alphanl);


        cxi1 =   gammar  * sqrt( I2(unl) ) * alphanl; 
        
        #tmp=project(I2(unl) ,Va)
        #print("mur=",mur," I2(unl) _min=",tmp.vector().min(),"max=",tmp.vector().max())
        #tmp=project(RHSalpha,Va)
        #print("RHSalpha1_min=",tmp.vector().min(),"max=",tmp.vector().max())
        #print('RHSalpha)=',RHSalpha)
        
    
        RHSalpha= 1/C/dt*alphan + mur * I2(unl) 
        if (VersionDivEta==1):
            RHSalpha = (RHSalpha+ gammar*    div(unl)   * sqrt(I2(unl)))*Max(Xi(unl)-xi0,0)/(abs(Xi(unl) - xi0)+1e-10 )
        elif (VersionDivEta==2) :
            RHSalpha = RHSalpha+ gammar *abs(div(unl))  * sqrt(I2(unl));
        elif (VersionDivEta==3) :
            RHSalpha = RHSalpha+ gammar *Max(div(unl),0)* sqrt(I2(unl));
        if (UsoEpsilonUpunto==1) :
            muequiv= 2.0*mu0  - alphanl*(2*mur + gammar*Xi(unl) ) + 1/dt * etaalpha
            RHSalpha = RHSalpha - detaalpha* inner(epsilon(unl), epsilon(unl-un)) /dt
        else :
            muequiv= 2.0*mu0  - alphanl*(2*mur + gammar*Xi(unl) ) +  etaalpha
            RHSalpha = RHSalpha - detaalpha* inner(epsilon(unl),epsilon(unl))             


        tmp=project(muequiv,Va)
        print("muequiv_min=",tmp.vector().min(),"max=",tmp.vector().max())
        tmp=project(RHSalpha,Va)
        print("RHSalpha_min=",tmp.vector().min(),"max=",tmp.vector().max())
        # Compute solution
        
        FRHSalpha.assign(project(RHSalpha,Va))
        Fetaalpha.assign(project(etaalpha,Va))
        Fcxi1.assign(project(cxi1,Va))
        Fmuequiv.assign(project(muequiv,Va))
        
        ua = Function(Vua)
        solve(Bilineal == Lineal, ua, bc)
        u,alpha=ua.split(deepcopy=True)
        u_1_, u_2_, u_3_ = u.split(deepcopy=True)
        
        print("\tuxmin=",u_1_.vector().min(),"max=",u_1_.vector().max())
        print("\tuymin=",u_2_.vector().min(),"max=",u_2_.vector().max())
        print("\tuzmin=",u_3_.vector().min(),"max=",u_3_.vector().max())
        print("\talmin=",alpha.vector().min(),"max=",alpha.vector().max())
        
        
        error_L2= norm(alpha.vector()-alphanl.vector(), 'linf')     \
                +norm(u.vector()-unl.vector(), 'linf')
        print('error_L2 = ', error_L2)
        errnl=error_L2
        unl.assign(u)
        alphanl.assign(alpha)
        
        
    # Update solution
    uo.assign(un) 
    un.assign(u)        
    alphan.assign(alpha)
    
        
    zmaxdet.append(u_3_.vector().min());
    alphamax .append(alpha.vector().max());
    alphamin .append(alpha.vector().min());
    tiempos.append(time);
    
    #print("zmaxdet=",zmaxdet)
    plt.figure(1,figsize=(20, 10))
    plt.clf()
    
    plt.subplot(211)
    plt.plot(tiempos,zmaxdet, 'bo-')
    plt.ylabel('Zmin')
    plt.grid()
    plt.subplot(212)
    plt.plot(tiempos,alphamin, 'r.-')
    plt.plot(tiempos,alphamax, 'b.-')
    plt.ylabel('Alpha Min/Max')
    plt.xlabel('Tiempo')
    plt.grid()
    plt.pause(0.005)
    #plt.show()



print(epsilon(u))

z=Expression('x[2]',degree=1)
V0 = FunctionSpace(mesh, 'P', 1)
z=project(z,V0)
plot(z)
plot(u)

interactive()
    