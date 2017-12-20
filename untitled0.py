# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 15:57:56 2017

@author: user1
"""

#%%

ds = Measure("ds")[boundaries]
for caso in range(2,3):
    if caso == 1:
        un = 1.0
        
    elif caso == 2:
        un = DeformaExpression(degree=2)
        un.set_funcdeforma(DeformaCaso2)
        un = interpolate(un, Q)
    
    elif caso == 3:
        un = DeformaExpression(degree=2)
        un.set_funcdeforma(DeformaCaso3)
        un = interpolate(un, Q)
    
    #Dn(f) := dot(grad(f), n)
    
    F = ( lambda_*div(up)*div(vp) \
        + mu2*inner(epsilon(up),epsilon(vp)) )*dx \
        - ( -un*inner(dot(elem_op(Dn, sigma(u)), n), vp)
            + inner(dot(sigma(u),grad(un)-inner(grad(un),n)*n), vp)
          )*ds(CAVEUP) \
        - ( -un*inner(dot(elem_op(Dn, sigma(u)), n), vp)
            + inner(dot(sigma(u),grad(un)-inner(grad(un),n)*n), vp)
          )*ds(CAVEMID) \
        - ( -un*inner(dot(elem_op(Dn, sigma(u)), n), vp)
            + inner(dot(sigma(u),grad(un)-inner(grad(un),n)*n), vp)
          )*ds(CAVEBOTTOM)
     
    a, L = lhs(F), rhs(F)

    print ("resolviendo la derivada para caso=",caso)
    
    # Assemble system, applying boundary conditions and preserving
    # symmetry)
    A, b = assemble_system(a, L, bcs)

    # Create solution function
    du = Function(V, name="du")

    # Set matrix operator
    solver.set_operator(A);
    
    # Compute solution
    solver.solve(du.vector(), b)
    
    Eder = project(lambda_*div(u)*div(du) + mu2*inner(epsilon(u),epsilon(du)), Q)
    dEpoten = project(-inner(f,du), Q)
    dEdiv = project( (lambda_ + mu2/3.0)*div(u)*div(du), Q)
    dEse = project( mu2*inner(Sepsilon(u),Sepsilon(du)), Q)
           
    print ("Grabando caso=", caso)
           
    if caso == 2:
       File(file_solution_UP.format('du'), "compressed") << du
       File(file_solution_UP.format('un'), "compressed") << un
       File(file_solution_UP.format('Eder'), "compressed") << Eder
       File(file_solution_UP.format('dEpoten'), "compressed") << dEpoten
       File(file_solution_UP.format('dEdiv'), "compressed") << dEdiv
       File(file_solution_UP.format('dEse'), "compressed") << dEse
           
    elif caso == 3:
       File(file_solution_Right.format('du'), "compressed") << du
       File(file_solution_Right.format('un'), "compressed") << un
       File(file_solution_Right.format('Eder'), "compressed") << Eder
       File(file_solution_Right.format('dEpoten'), "compressed") << dEpoten
       File(file_solution_Right.format('dEdiv'), "compressed") << dEdiv
       File(file_solution_Right.format('dEse'), "compressed") << dEse
    
"""
	{ofstream fu1p(pathResult + "u1p_"+caso+".txt");fu1p<<u1p[];}
	{ofstream fu2p(pathResult + "u2p_"+caso+".txt");fu2p<<u2p[];}
	{ofstream fu3p(pathResult + "u3p_"+caso+".txt");fu3p<<u3p[];}
/**/"""