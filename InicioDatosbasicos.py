# -*- coding: utf-8 -*-

from __future__ import print_function
from fenics import *
from mshr import *
import numpy as np
import matplotlib.pyplot as plt
import time

time0=time.time()


VersionEta=1              #0: no hay eta, 
                        #1: Se usa la función eta(alpha) del articulo. 
                        #2: Se usa la función eta(alpha) propuesta por JSM
UsoEpsilonUpunto=0 #1: Se usa el modelo con el término \eta(\alpha)*E(u):E(udot). 
                        #0: Se usa E(u) en lugar de E(udot) en dicha parte del modelo "propuesto por JSM"
VersionDivEta=2    #1: Se usa termino: RHSalpha+= gammar *    div(u1nl,u2nl,u3nl)             * sqrt( I2(u1nl,u2nl,u3nl) );
                        #2: Se usa término: RHSalpha+= gammar *    abs(div(u1nl,u2nl,u3nl))        * sqrt( I2(u1nl,u2nl,u3nl) );
                        #3: Se usa término: RHSalpha+= gammar *    max(div(u1nl,u2nl,u3nl),0*u1nl) * sqrt( I2(u1nl,u2nl,u3nl) )
UsoDalphaDt=1;                         
#runfile("Defaults.py")
tm=time.localtime(time.time())
prefix=	""+('D%d-%d.%d_' % (tm.tm_yday,tm.tm_hour,tm.tm_min))+`VersionEta`+`UsoEpsilonUpunto`+`VersionDivEta`+`UsoDalphaDt`+"_";			

print ('prefix=',prefix)
SVersionEta=['No hay Eta' ,\
             'Se usa Eta del Articulo' , \
             'Se usa EtaM' ]
SUsoEpsilonUpunto = ['Se usa E(u) en lugar de E(udot)', \
                     'Se usa E(udot)' ]
SVersionDivEta=['0 no se usa', 
               'RHSalpha: se Usa div(unl)',
               'RHSalpha: se Usa |div(unl)|',
               'RHSalpha: se Usa max{0,div(unl)}',]
SUsoDalphaDt=['Alpha=cte=0', 
               'Alpha evoluciona']
print (VersionEta,':',SVersionEta[VersionEta])
print (UsoEpsilonUpunto,':',SUsoEpsilonUpunto[UsoEpsilonUpunto])
print (VersionDivEta,':',SVersionDivEta[VersionDivEta])
print (UsoDalphaDt,':',SUsoDalphaDt[UsoDalphaDt])

#-----------------------------------------
# Model parameters.
#-----------------------------------------
MP = 1e6 ;
mu0    = 30000 * MP ;
lambda_ = 30000 * MP;
rho    = 2700;
E  = mu0 *( 3 * lambda_ + 2 * mu0) / (lambda_ + mu0); 
nu = lambda_ /  (2 * ( lambda_ + mu0)); 
Cgravedad = 9.8e0;
xi0 = -0.008;  
sigmaP=10*MP/0.05;
   


gammar = 35000 * MP; 
mur =  0.8 * gammar;

gammar /= 100
kappaD = 1e6;
C = 1e-6; 
eta1 = 1e-2;
eta2 = 1e-1;


tEnd = 50;
dt = 0.001 ;