# -*- coding: utf-8 -*-
"""

"""
relaj=0.05

alpha_nl2=alpha_nl
u_nl2=u_nl
errnl=1;errnl2=10;itnl=0
while ((errnl>1e-8 or itnl<2) and itnl<5000000 ) :
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
        if alphasolo.vector().max()>ac+0.1 or alphasolo.vector().max()<ac-0.1:
            kappaD_exp.kappaD *= alphasolo.vector().max()/ac
            print "\t\033[6;30;42mNuevo kappaD=",kappaD_exp.kappaD, '\033[0m'
            
            solve(BilinealEA == LinealEA, alphasolo, bcalpha)
            print "\t\033[6;30;42malplaSolo: min,max=",alphasolo.vector().min(),",",alphasolo.vector().max(), '\033[0m'
        plot(alphasolo, key='alphasolo',title='alphasolo')
        if primeravez :
            interactive()
            primeravez=0


    ################### FIN>resuelve solo alpha #########################

    
    A, b = assemble_system(BilinealE, LinealE, bc)
    print 'Norma(A)=',A.norm('linf'),' Norma(b)=',b.norm('linf')
    
    
    
    
    ua = Function(Vua)
    solve(BilinealE == LinealE, ua, bc)
    uu,aalpha=ua.split(deepcopy=True)
    u.assign(uu)
    alpha.assign(aalpha)
    alpha.vector()[:]=Min(alpha.vector().array(),ac)
    plot(alpha, key='alpha',title='alpha')
    u_1_, u_2_, u_3_ = u.split(deepcopy=True)
    
    print "\tuxmin=",u_1_.vector().min(),"max=",u_1_.vector().max()
    print "\tuymin=",u_2_.vector().min(),"max=",u_2_.vector().max()
    print "\tuzmin=",u_3_.vector().min(),"max=",u_3_.vector().max()
    print "\talmin=",alpha.vector().min(),"max=",alpha.vector().max()
    
    
    errnl= norm(alpha.vector()-alpha_nl.vector(), 'linf')     \
            +norm(u.vector()-u_nl.vector(), 'linf')
    print 'errnl = \033[6;30;42m', errnl  , '\033[0m, errnl2 = \033[6;30;42m', errnl2  , '\033[0m, FactorGrande =',FactorGrande.vector().max()
    
    errnl_arr.append(errnl)
    fig = plt.figure(1,figsize=(15, 10))
    plt.clf()
    plt.title("Errores")
    plt.plot(np.log10(errnl_arr), 'bo-')
    plt.ylabel('ErrNL')
    plt.grid()
    fig.tight_layout()
    #fig.set_size_inches(w=11,h=7)
    plt.pause(0.00005)
#    fig_name = prefix+('Modelo3_[dt=%g]'%dt)+'.png'
#    fig.savefig(fig_name)
    #plt.show()
    
    if (errnl>errnl2*1.01) :
        relaj = max(relaj/1.2,.01)
        u_nl.assign(u_nl2)
        alpha_nl.assign(alpha_nl2)
        print "\033[6;30;42mrelajUP=",relaj, '\033[0m' 
    else:
        if (errnl<errnl2*0.99) :
            relaj=min(relaj*1.1,1)    

        if (errnl<1e-8 and FactorGrande.vector().max()>1e-12) :
            if FactorGrande.vector().max()>3000 :
                FactorGrande.assign(project(FactorGrande*0.6,Va))
            else :
                FactorGrande.assign(project(FactorGrande*0.8,Va))
            errnl=1;
            relaj=0.5
        if (alpha.vector().max()>1 or alpha.vector().min<-0.1) :
            print '\033[6;30;43m Fuera de Rango ==> Corrijo \033[0m' 
            #FactorGrande.assign(project(FactorGrande*1.8,Va))
            u.assign(u_nl)
            alpha.assign(alpha_nl)
            u_1_, u_2_, u_3_ = u.split(deepcopy=True)
            print "\tuxmin=",u_1_.vector().min(),"max=",u_1_.vector().max()
            print "\tuymin=",u_2_.vector().min(),"max=",u_2_.vector().max()
            print "\tuzmin=",u_3_.vector().min(),"max=",u_3_.vector().max()
            print "\talmin=",alpha.vector().min(),"max=",alpha.vector().max()
            errnl=1;
            
        u_nl2.assign(u_nl)
        alpha_nl2.assign(alpha_nl)
        print "\033[6;30;42mrelaj=",relaj, '\033[0m' 
        u_nl.assign(project(u_nl+(u-u_nl)*relaj,Vu))
        alpha_nl.assign(project(alpha_nl+(alpha-alpha_nl)*relaj,Va))

    errnl2=errnl    
