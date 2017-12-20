# -*- coding: utf-8 -*-

execfile("InicioDatosbasicos.py")
execfile("Malla100x100x100.py")

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
FnEtaPrima= Function(Va)
up       = Function(Vu)
u       = Function(Vu)			## No es parte de la FV
u_nl      = Function(Vu)
udp       = Function(Vu)
ud_nl      = Function(Vu)
Fuerza      = Function(Vu)
FuerzaDot      = Function(Vu)

execfile("Macros.py")
execfile("Malla100x100x100_BC.py")



ZeroVector = Constant((0,0,0))
ZeroEscalar= Constant(0)



execfile("FV.py")

        
####################Problema Estacionario ###################
    
Fuerza.assign(project(Constant((0, 0, -rho*Cgravedad)),Vu)) 

u_nl   .assign(project(ZeroVector , Vu))
alpha_nl.assign(project(ZeroEscalar, Va))

#%%

errnl=1;errnl2=10;relaj=1;itnl=0
while ((errnl>1e-8 or itnl<2) and itnl<5000 ) :
    itnl+=1
    print "itnl = ",itnl
    if (VersionEta==1) :
        print 'VersionEta==1'
        etaalpha=eta(alpha_nl);
        detaalpha=deta(alpha_nl);
    elif (VersionEta==2) :
        print 'VersionEta==2'
        etaalpha=etaM(alpha_nl);
        detaalpha=detaM(alpha_nl);

### FnEtaPrima 
    ArrAlpha=alpha_nl.vector().array()
    if (VersionEta==1) :
        Arrdetaalpha=deta(ArrAlpha);
    elif (VersionEta==2) :
        Arrdetaalpha=detaM(ArrAlpha);
    ArrEtaPrima=0*Arrdetaalpha;
    if (UsoEpsilonUpunto==0) :
        tmp=project(inner(epsilon(u_nl),epsilon(u_nl)),Va)
        ArrEtaPrima = Arrdetaalpha* tmp.vector().array()/(abs(ArrAlpha)+1e-3)  
        ArrEtaPrima = Arrdetaalpha/(abs(ArrAlpha)+1e-3)  
    FnEtaPrima.vector()[:]=ArrEtaPrima
    print "ArrAlpha    :   min,max=\033[6;30;42m",ArrAlpha.min(),",",ArrAlpha.max(), '\033[0m'
    print "Arrdetaalpha:   min,max=\033[6;30;42m",Arrdetaalpha.min(),",",Arrdetaalpha.max(), '\033[0m'
    print "ArrEtaPrima :   min,max=\033[6;30;42m",ArrEtaPrima.min(),",",ArrEtaPrima.max(), '\033[0m'
    print "FnEtaPrima  : min,max=",FnEtaPrima.vector().min(),",",FnEtaPrima.vector().max()
##################

    u_1_, u_2_, u_3_ = Fuerza.split(deepcopy=True)
    print("f3 : min=",u_3_.vector().min(),"max=",u_3_.vector().max())
    
    cxi1 =   gammar  * sqrt( I2(u_nl) ) * alpha_nl; 
    

    RHSalpha=  mur * I2(u_nl) 
    if (VersionDivEta==1):
        RHSalpha = (RHSalpha+ gammar*    div(u_nl)   * sqrt(I2(u_nl)))*Max(Xi(u_nl)-xi0,0)/(abs(Xi(u_nl) - xi0)+1e-10 )
    elif (VersionDivEta==2) :
        RHSalpha = RHSalpha+ gammar *abs(div(u_nl))  * sqrt(I2(u_nl));
    elif (VersionDivEta==3) :
        RHSalpha = RHSalpha+ gammar *Max(div(u_nl),0)* sqrt(I2(u_nl));
    if (UsoEpsilonUpunto==1) :
        muequiv= 2.0*mu0  - alpha_nl*(2*mur + gammar*Xi(u_nl) ) 
    else :
        muequiv= 2.0*mu0  - alpha_nl*(2*mur + gammar*Xi(u_nl) ) +  etaalpha        


    tmp=project(muequiv-2.0*mu0 ,Va)
    print "muequiv_min=",tmp.vector().min(),"max=",tmp.vector().max()
    tmp=project(RHSalpha,Va)
    print "RHSalpha_min=",tmp.vector().min(),"max=",tmp.vector().max()
    # Compute solution
    
    FRHSalpha.assign(project(RHSalpha,Va))
    Fetaalpha.assign(project(etaalpha,Va))
    Fcxi1.assign(project(cxi1,Va))
    Fmuequiv.assign(project(muequiv,Va))



    ################### resuelve solo alpha #########################

    if (itnl>1) :
        alphasolo=Function(Va)
        solve(BilinealEA == LinealEA, alphasolo, bcalpha)
        print "\t\033[6;30;42malplaSolo: min,max=",alphasolo.vector().min(),",",alphasolo.vector().max(), '\033[0m'
        if alphasolo.vector().max()>1 :
            kappaD_exp.kappaD *= alphasolo.vector().max()
            print "\t\033[6;30;42mNuevo kappaD=",kappaD_exp.kappaD, '\033[0m'
            
        solve(BilinealEA == LinealEA, alphasolo, bcalpha)
        print "\t\033[6;30;42malplaSolo: min,max=",alphasolo.vector().min(),",",alphasolo.vector().max(), '\033[0m'

    ################### resuelve solo alpha #########################

    
    A, b = assemble_system(BilinealE, LinealE, bc)
    print 'Norma(A)=',A.norm('linf'),' Norma(b)=',b.norm('linf')
    
    
    
    
    ua = Function(Vua)
    solve(BilinealE == LinealE, ua, bc)
    u,alpha=ua.split(deepcopy=True)
    u_1_, u_2_, u_3_ = u.split(deepcopy=True)
    
    print "\tuxmin=",u_1_.vector().min(),"max=",u_1_.vector().max()
    print "\tuymin=",u_2_.vector().min(),"max=",u_2_.vector().max()
    print "\tuzmin=",u_3_.vector().min(),"max=",u_3_.vector().max()
    print "\talmin=",alpha.vector().min(),"max=",alpha.vector().max()
    
    
    errnl= norm(alpha.vector()-alpha_nl.vector(), 'linf')     \
            +norm(u.vector()-u_nl.vector(), 'linf')
    print 'errnl = \033[6;30;42m', errnl  , '\033[0m, FactorGrande =',FactorGrande.vector().max()
    if (errnl<1e-8 and FactorGrande.vector().max()>1e-12) :
        if FactorGrande.vector().max()>3000 :
            FactorGrande.assign(project(FactorGrande*0.6,Va))
        else :
            FactorGrande.assign(project(FactorGrande*0.95,Va))
        errnl=1;
    if (alpha.vector().max()>1 or alpha.vector().min<-0.1) :
        print '\033[6;30;43m Fuera de Rango ==> Corrijo \033[0m' 
        FactorGrande.assign(project(FactorGrande*1.8,Va))
        u.assign(u_nl)
        alpha.assign(alpha_nl)
        u_1_, u_2_, u_3_ = u.split(deepcopy=True)
        print "\tuxmin=",u_1_.vector().min(),"max=",u_1_.vector().max()
        print "\tuymin=",u_2_.vector().min(),"max=",u_2_.vector().max()
        print "\tuzmin=",u_3_.vector().min(),"max=",u_3_.vector().max()
        print "\talmin=",alpha.vector().min(),"max=",alpha.vector().max()
        errnl=1;
    if FactorGrande.vector().max()<3000 :
        u_nl.assign(project(u_nl+(u-u_nl)*0.2,Vu))
        alpha_nl.assign(project(alpha_nl+(alpha-alpha_nl)*0.2,Va))
    else :
        u_nl.assign(project(u_nl+(u-u_nl)*1,Vu))
        alpha_nl.assign(project(alpha_nl+(alpha-alpha_nl)*1,Va))
    
    



##################### Problema de Evolucion #################
# Initial condition


up    .assign(project(ZeroVector , Vu))
udp   .assign(project(ZeroVector , Vu))
alphap.assign(project(ZeroEscalar, Va))
up.assign(u)
alphap.assign(alpha)

tiempo = 0
kiter = 0
etaalpha=0;detaalpha=0
nt = int(tEnd/dt)

zmaxdet=[up.vector().min()];alphamax=[alphap.vector().max()];alphamin=[alphap.vector().min()];
tiempos=[0];Niter_nl=[0];Residuos=[0];
#%%  Iteraciones Temporales
for i in range(0,nt*0+1):
    print '\n=================\nIteración ',i,'\n================='
    
    tiempo = tiempo + dt; 
    kiter+=1
    print "time = ",tiempo," dt = ",dt
    
    # Non-linear terms are implemented by fixed-point iterations (u_nl...alpha_nl)
    
    fv = Constant((0, 0, -rho*Cgravedad)) + (rho/dt)*udp
    FuerzaDot.assign(project(fv,Vu)) 
    fv = udp*rho+up*rho/dt
    Fuerza.assign(project(fv,Vu)) 

    u_1_, u_2_, u_3_ = Fuerza.split(deepcopy=True)
    print("f3 : min=",u_3_.vector().min(),"max=",u_3_.vector().max())
    u_1_, u_2_, u_3_ = FuerzaDot.split(deepcopy=True)
    print("fd3: min=",u_3_.vector().min(),"max=",u_3_.vector().max())
    
    ud_nl.assign(udp)
    u_nl.assign(up)
    alpha_nl.assign(alphap)
    errnl=1;errnl2=10;relaj=1;itnl=0
    while ((errnl>1e-8 or itnl<2) and itnl<5 ) :
        itnl+=1
        u.assign(project(up+(dt/2)*(udp+ud_nl),Vu))
        u_nl.assign(u)
        print "itnl = ",itnl
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


        tmp=project(muequiv-2.0*mu0,Va)
        print "muequiv*_min=",tmp.vector().min(),"max=",tmp.vector().max()
        tmp=project(RHSalpha-1/C/dt*alphap ,Va)
        print "RHSalpha*_min=",tmp.vector().min(),"max=",tmp.vector().max()
        # Compute solution
        
        FRHSalpha.assign(project(RHSalpha,Va))
        Fetaalpha.assign(project(etaalpha,Va))
        Fcxi1.assign(project(cxi1,Va))
        Fmuequiv.assign(project(muequiv,Va))
        
        
        A, b = assemble_system(Bilineal, Lineal, bc)
        print 'Norma(A)=',A.norm('linf'),' Norma(b)=',b.norm('linf')
        ua = Function(Vua)
        solve(Bilineal == Lineal, ua, bc)
        ud,alpha=ua.split(deepcopy=True)
        u_1_, u_2_, u_3_ = ud.split(deepcopy=True)
        
        print "\tudx : min=",u_1_.vector().min(),"max=",u_1_.vector().max()
        print "\tudy : min=",u_2_.vector().min(),"max=",u_2_.vector().max()
        print "\tudz : min=",u_3_.vector().min(),"max=",u_3_.vector().max()
        print "\talph: min=",alpha.vector().min(),"max=",alpha.vector().max()
        
        u_1_, u_2_, u_3_ = u.split(deepcopy=True)
        
        print "\tux : min=",u_1_.vector().min(),"max=",u_1_.vector().max()
        print "\tuy : min=",u_2_.vector().min(),"max=",u_2_.vector().max()
        print "\tuz : min=",u_3_.vector().min(),"max=",u_3_.vector().max()
        
        
        errnla= norm(alpha.vector()-alpha_nl.vector(), 'linf')     
        errnlu= norm(ud.vector()-ud_nl.vector(), 'linf')
        errnl= errnlu+errnla
        print 'errnl = \033[6;30;42m ', errnl  , '=',errnlu,'+',errnla,'\033[0m'
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
    Niter_nl.append(itnl);
    Residuos.append(errnl);
    
    #print("zmaxdet=",zmaxdet)
    fig = plt.figure(1,figsize=(15, 10))
    plt.clf()
    
    plt.subplot2grid((6,1), (0,0),  rowspan=2)
    plt.title("Modelo3: Cputime=%.2fs (%.2f s/TimeStep), [dt=%g]" % ((time.time()-time0),(time.time()-time0)/kiter,dt))
    plt.plot(tiempos,zmaxdet, 'bo-')
    plt.ylabel('Zmin')
    plt.grid()
    plt.subplot2grid((6,1), (2,0),  rowspan=2)
    plt.plot(tiempos,alphamin, 'r.-')
    plt.plot(tiempos,alphamax, 'b.-')
    plt.ylabel('Alpha Min/Max')
    plt.grid()
    plt.subplot(615)
    plt.plot(tiempos,Niter_nl, 'b.-')
    plt.ylabel('Niter_NL')
    #plt.xlabel('Tiempo')
    plt.grid()
    plt.subplot(616)
    plt.plot(tiempos,Residuos, 'b.-')
    plt.ylabel('ErrNL')
    plt.xlabel('Tiempo')
    plt.grid()
    fig.tight_layout()
    #fig.set_size_inches(w=11,h=7)
    plt.pause(0.00005)
    fig_name = prefix+('Modelo3_[dt=%g]'%dt)+'.png'
    fig.savefig(fig_name)
    #plt.show()




z=Expression('x[2]',degree=1)
V0 = FunctionSpace(mesh, 'P', 1)
z=project(z,V0)
plot(z)
plot(u_nl)
plot(alpha_nl)

interactive()
    