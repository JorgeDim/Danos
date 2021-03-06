# -*- coding: utf-8 -*-
tol = 1E-1
def lados(x, on_boundary):
    return on_boundary and (near(x[0],x0,tol) or near(x[0],x1,tol))
def fondo(x, on_boundary):
    return on_boundary and (near(x[2],z0,tol) )
def frentefondo(x, on_boundary):
    return on_boundary and (near(x[1],y0,tol) or near(x[1],y1,tol))

bc1 = DirichletBC(Vvua.sub(0).sub(0), Constant((0)), lados)
bc2 = DirichletBC(Vvua.sub(0).sub(2), Constant((0)), fondo)
bc3 = DirichletBC(Vvua.sub(0).sub(1), Constant((0)), frentefondo)
bc4 = DirichletBC(Vvua.sub(1).sub(0), Constant((0)), lados)
bc5 = DirichletBC(Vvua.sub(1).sub(2), Constant((0)), fondo)
bc6 = DirichletBC(Vvua.sub(1).sub(1), Constant((0)), frentefondo)
bc=[bc1,bc2,bc3,bc4,bc5,bc6]
