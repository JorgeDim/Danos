from dolfin import *
from mshr import *
h = 0.25
r = 0.3*h
box = Box(Point(0, 0, 0), Point(1, h, h))
s0 = Sphere(Point(0.3, 0.50*h, 0.50*h), r)
s1 = Sphere(Point(0.5, 0.65*h, 0.65*h), r)
s2 = Sphere(Point(0.7, 0.35*h, 0.35*h), r)
domain = box - s0 - s1 - s2
 
# Generate mesh
mesh = generate_mesh(domain, 32)

# Define function space
P2 = VectorElement('P', tetrahedron, 2)
P1 = FiniteElement('P', tetrahedron, 1)
TH = P2 * P1
W = FunctionSpace(mesh, TH)
 
# Define variational problem
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)
a = inner(grad(u), grad(v))*dx - p*div(v)*dx + div(u)*q*dx
L = dot(f, v)*dx
 
# Compute solution
w = Function(W)
solve(a == L, w, [bc1, bc0])



