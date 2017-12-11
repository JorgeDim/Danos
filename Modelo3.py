# -*- coding: utf-8 -*-

runfile("InicioDatosbasicos.py")
runfile("Malla100x100x100.py")

#%% Espacios Elementos Finitos
P1v = VectorElement('P', tetrahedron, 1)
P1e = FiniteElement('P', tetrahedron, 1)
element =MixedElement([P1v, P1e])
Vua = FunctionSpace(mesh, element)
Vu = FunctionSpace(mesh, P1v)
Va= FunctionSpace(mesh, P1e)

alphap  = Function(Va)
alpha_nl  = Function(Va)
Fmuequiv  = Function(Va)
Fetaalpha = Function(Va)
Fcxi1     = Function(Va)
FRHSalpha = Function(Va)
up       = Function(Vu)
u_nl      = Function(Vu)
udp       = Function(Vu)
ud_nl      = Function(Vu)
Fuerza      = Function(Vu)
FuerzaDot      = Function(Vu)

runfile("Macros.py")
runfile("Malla100x100x100_BC.py")



ZeroVector = Constant((0,0,0))
ZeroEscalar= Constant(0)


vd,w = TestFunctions(Vua)
ud, alpha=TrialFunctions(Vua)
dim = ud.geometric_dimension()  # space dimension

Bilineal = (rho/dt)*inner(ud,vd)*dx \
    + inner(sigmaL((u_nl+up)/2), epsilon(vd))*dx   \
    +(1/C/dt)*alpha*w *dx \
    + inner(kappaD*grad(alpha),grad(w))*dx 
    
    
if (UsoEpsilonUpunto==1) :
    Bilineal= Bilineal + Fetaalpha*inner(epsilon(ud), epsilon(vd))*dx
else :
    Bilineal= Bilineal + Fetaalpha*inner(epsilon((u_nl+up)/2), epsilon(vd))*dx
  

L1= dot(FuerzaDot, vd)*dx 
L2 = Fcxi1 * div(vd) *dx
L3=0
if UsoEpsilonUpunto != 0:
    L3 = (1/dt*Fetaalpha)*inner(epsilon(up),epsilon(v))*dx 
L4=0
if UsoDalphaDt!=0:
    L4 = FRHSalpha *w*dx 
Lineal=L1+L2+L3+L4
F=Bilineal-Lineal
Bilineal,Lineal = lhs(F), rhs(F)

        

# Initial condition


up    .assign(project(ZeroVector , Vu))
udp   .assign(project(ZeroVector , Vu))
alphap.assign(project(ZeroEscalar, Va))

tiempo = 0
kiter = 0
etaalpha=0;detaalpha=0
nt = int(tEnd/dt)

zmaxdet=[0];
alphamax=[0];
alphamin=[0];
tiempos=[0];
#%%  Iteraciones Temporales
for i in range(0,nt):
    print("\n=================\nIteración ",i,"\n=================")
    
    tiempo = tiempo + dt; 
    kiter+=1
    print(" time = ",tiempo," dt = ",dt)
    
    # Non-linear terms are implemented by fixed-point iterations (u_nl...alpha_nl)
    
    fv = Constant((0, 0, -rho*Cgravedad)) + udp*rho/dt
    FuerzaDot.assign(project(fv,Vu)) 
    fv = udp*rho+up*rho/dt
    Fuerza.assign(project(fv,Vu)) 

    #u_1_, u_2_, u_3_ = Fuerza.split(deepcopy=True)
    #print("f3_min=",u_3_.vector().min(),"max=",u_3_.vector().max())
    
    ud_nl.assign(udp)
    u_nl.assign(up)
    alpha_nl.assign(alphap)
    errnl=1;errnl2=10;relaj=1;itnl=0
    while ((errnl>1e-8 or itnl<2) and itnl<50 ) :
        itnl+=1
        u_nl.assign(project(up+(dt/2)*(udp+ud_nl),Vu))
        print("itnl = ",itnl) 
        if (VersionEta==1) :
            etaalpha=eta(alpha_nl);
            detaalpha=deta(alpha_nl);
        elif (VersionEta==2) :
            etaalpha=etaM(alpha_nl);
            detaalpha=detaM(alpha_nl);


        cxi1 =   gammar  * sqrt( I2(u_nl) ) * alpha_nl; 
        
    
        RHSalpha= 1/C/dt*alphap + mur * I2(u_nl) 
        if (VersionDivEta==1):
            RHSalpha = (RHSalpha+ gammar*    div(u_nl)   * sqrt(I2(u_nl)))*Max(Xi(u_nl)-xi0,0)/(abs(Xi(u_nl) - xi0)+1e-10 )
        elif (VersionDivEta==2) :
            RHSalpha = RHSalpha+ gammar *abs(div(u_nl))  * sqrt(I2(u_nl));
        elif (VersionDivEta==3) :
            RHSalpha = RHSalpha+ gammar *Max(div(u_nl),0)* sqrt(I2(u_nl));
        if (UsoEpsilonUpunto==1) :
            muequiv= 2.0*mu0  - alpha_nl*(2*mur + gammar*Xi(u_nl) ) + 1/dt * etaalpha
            RHSalpha = RHSalpha - detaalpha* inner(epsilon(u_nl), epsilon(ud_nl)) 
        else :
            muequiv= 2.0*mu0  - alpha_nl*(2*mur + gammar*Xi(u_nl) ) +  etaalpha
            RHSalpha = RHSalpha - detaalpha* inner(epsilon(u_nl),epsilon(u_nl))             


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
        ud,alpha=ua.split(deepcopy=True)
        u_1_, u_2_, u_3_ = ud.split(deepcopy=True)
        
        print("\tuxmin=",u_1_.vector().min(),"max=",u_1_.vector().max())
        print("\tuymin=",u_2_.vector().min(),"max=",u_2_.vector().max())
        print("\tuzmin=",u_3_.vector().min(),"max=",u_3_.vector().max())
        print("\talmin=",alpha.vector().min(),"max=",alpha.vector().max())
        
        
        error_L2= norm(alpha.vector()-alpha_nl.vector(), 'linf')     \
                +norm(ud.vector()-ud_nl.vector(), 'linf')
        print('error_L2 = ', error_L2)
        errnl=error_L2
        ud_nl.assign(ud)
        u_nl.assign(project(up+(dt/2)*(udp+ud_nl),Vu))
        alpha_nl.assign(alpha)
        
        
    # Update solution
    udp.assign(ud) 
    up.assign(u_nl)        
    alphap.assign(alpha)
    
        
    u_1_, u_2_, u_3_ = u_nl.split(deepcopy=True)
    zmaxdet.append(u_3_.vector().min());
    alphamax .append(alpha.vector().max());
    alphamin .append(alpha.vector().min());
    tiempos.append(tiempo);
    
    #print("zmaxdet=",zmaxdet)
    plt.figure(1,figsize=(20, 10))
    plt.clf()
    
    plt.subplot(211)
    plt.title("Cputime=%.2fs (%.2f/TimeStep)" % ((time.time()-time0),(time.time()-time0)/kiter))
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




z=Expression('x[2]',degree=1)
V0 = FunctionSpace(mesh, 'P', 1)
z=project(z,V0)
plot(z)
plot(u_nl)

interactive()
    