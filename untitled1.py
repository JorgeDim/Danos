# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 16:35:23 2017

@author: user1
"""

#%% Test 


valEtaPrima = detaalpha* inner(epsilon(u_nl),epsilon(u_nl))/(abs(alpha_nl)+1e-6)   

FnEtaPrima.assign(project(valEtaPrima,Va))

print "FnEtaPrima: min=",FnEtaPrima.vector().min(),"max=",FnEtaPrima.vector().max()