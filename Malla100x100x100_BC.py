# -*- coding: utf-8 -*-
tol = 1E-1
def lados(x, on_boundary):
    return on_boundary and (near(x[0],x0,tol) or near(x[0],x1,tol))
def fondo(x, on_boundary):
    return on_boundary and (near(x[2],z0,tol) )
def arriba(x, on_boundary):
    return on_boundary and (near(x[2],z1,tol) )
def frentefondo(x, on_boundary):
    return on_boundary and (near(x[1],y0,tol) or near(x[1],y1,tol))

bc1 = DirichletBC(Vua.sub(0).sub(0), Constant((0)), lados)
bc2 = DirichletBC(Vua.sub(0).sub(2), Constant((0)), fondo)
bc3 = DirichletBC(Vua.sub(0).sub(1), Constant((0)), frentefondo)
bca = DirichletBC(Vua.sub(1), Constant((0)), arriba)
bc=[bc1,bc2,bc3,bca]

bcalpha = DirichletBC(Va, Constant((0)), arriba)

bc1u = DirichletBC(Vu.sub(0), Constant((0)), lados)
bc2u = DirichletBC(Vu.sub(2), Constant((0)), fondo)
bc3u = DirichletBC(Vu.sub(1), Constant((0)), frentefondo)
bcu_solo=[bc1u,bc2u,bc3u]

