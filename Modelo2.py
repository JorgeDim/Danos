# -*- coding: utf-8 -*-

runfile("InicioDatosbasicos.py")
runfile("Malla100x100x100.py")

#%% Espacios Elementos Finitos
P1v = VectorElement('P', tetrahedron, 1)
P1e = FiniteElement('P', tetrahedron, 1)
element =MixedElement([P1v, P1v, P1e])
Vu = FunctionSpace(mesh, P1v)
Vvua = FunctionSpace(mesh, element)
Va= FunctionSpace(mesh, P1e)

alphan  = Function(Va)
alphanl  = Function(Va)
Fmuequiv  = Function(Va)
Fetaalpha = Function(Va)
Fcxi1     = Function(Va)
FRHSalpha = Function(Va)
un       = Function(Vu)
unl      = Function(Vu)
udn       = Function(Vu)
udnl      = Function(Vu)
Fuerza      = Function(Vu)
FuerzaDot      = Function(Vu)

runfile("Macros.py")
runfile("Malla100x100x100_BC_vua.py")



ZeroVector = Constant((0,0,0))
ZeroEscalar= Constant(0)


vd,v,w = TestFunctions(Vvua)
ud,u, alpha=TrialFunctions(Vvua)
dim = u.geometric_dimension()  # space dimension

Bilineal = (rho/dt)*inner(ud,vd)*dx \
    + inner(sigma0(u), epsilon(vd))*dx +   \
    + (rho/dt)*inner(u,v)*dx \
    +(1/C/dt)*alpha*w *dx \
    + inner(kappa*grad(alpha),grad(w))*dx 

L1 = dot(FuerzaDot, vd)*dx + dot(Fuerza, v)*dx  
L2 = Fcxi1 * div(v) *dx
L3=0
if UsoEpsilonUpunto != 0:
    L3 = (1/dt*Fetaalpha)*inner(epsilon(un),epsilon(v))*dx 
L4=0
if UsoDalphaDt!=0:
    L4 = FRHSalpha *w*dx 
Lineal=L1+L2+L3+L4

        

# Initial condition


un    .assign(project(ZeroVector , Vu))
udn   .assign(project(ZeroVector , Vu))
alphan.assign(project(ZeroEscalar, Va))

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
    
    fv = Constant((0, 0, -rho*Cgravedad)) + udn*rho/dt
    FuerzaDot.assign(project(fv,Vu)) 
    fv = udn*rho+un*rho/dt
    Fuerza.assign(project(fv,Vu)) 

    #u_1_, u_2_, u_3_ = Fuerza.split(deepcopy=True)
    #print("f3_min=",u_3_.vector().min(),"max=",u_3_.vector().max())
    
    udnl.assign(udn)
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
        
    
        RHSalpha= 1/C/dt*alphan + mur * I2(unl) 
        if (VersionDivEta==1):
            RHSalpha = (RHSalpha+ gammar*    div(unl)   * sqrt(I2(unl)))*Max(Xi(unl)-xi0,0)/(abs(Xi(unl) - xi0)+1e-10 )
        elif (VersionDivEta==2) :
            RHSalpha = RHSalpha+ gammar *abs(div(unl))  * sqrt(I2(unl));
        elif (VersionDivEta==3) :
            RHSalpha = RHSalpha+ gammar *Max(div(unl),0)* sqrt(I2(unl));
        if (UsoEpsilonUpunto==1) :
            muequiv= 2.0*mu0  - alphanl*(2*mur + gammar*Xi(unl) ) + 1/dt * etaalpha
            RHSalpha = RHSalpha - detaalpha* inner(epsilon(unl), epsilon(udnl)) 
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
        
        ua = Function(Vvua)
        solve(Bilineal == Lineal, ua, bc)
        ## 
        ud,u,alpha=ua.split(deepcopy=True)
        u_1_, u_2_, u_3_ = u.split(deepcopy=True)
        
        print("\tuxmin=",u_1_.vector().min(),"max=",u_1_.vector().max())
        print("\tuymin=",u_2_.vector().min(),"max=",u_2_.vector().max())
        print("\tuzmin=",u_3_.vector().min(),"max=",u_3_.vector().max())
        print("\talmin=",alpha.vector().min(),"max=",alpha.vector().max())
        
        
        error_L2= norm(alpha.vector()-alphanl.vector(), 'linf')     \
                +norm(u.vector()-unl.vector(), 'linf')
        print('error_L2 = ', error_L2)
        errnl=error_L2
        udnl.assign(ud)
        unl.assign(u)
        alphanl.assign(alpha)
        
        
    # Update solution
    udn.assign(ud) 
    un.assign(u)        
    alphan.assign(alpha)
    
        
    zmaxdet .append(u_3_.vector().min());
    alphamax.append(alpha.vector().max());
    alphamin.append(alpha.vector().min());
    tiempos .append(time);
    
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
    