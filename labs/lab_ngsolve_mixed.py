import ngsolve as ng

mesh = ng.Mesh(ng.unit_square.GenerateMesh(maxh=0.2))

V = ng.HDiv(mesh, order=0, RT=True)
P = ng.L2(mesh, order=0)

space = V * P

v, q = space.TestFunction()
u, p = space.TrialFunction()

normal = ng.specialcf.normal(mesh.dim)

g = ng.LinearForm(space)
g += ng.x * (v.Trace() * normal) * ng.ds
g.Assemble()

A = ng.BilinearForm(space)
A += (u * v) * ng.dx
A += -(p * ng.div(v)) * ng.dx
A += (q * ng.div(u)) * ng.dx
A.Assemble()

sol = ng.GridFunction(space)
sol.vec.data = A.mat.Inverse() * g.vec

sol_u, sol_p = sol.components

sol_u_x = sol_u[0]
sol_u_y = sol_u[1]
zero = ng.CoefficientFunction(0.0)
sol_u_3d = ng.CoefficientFunction((sol_u_x, sol_u_y, zero))

vtk = ng.VTKOutput(mesh, [sol_u_3d, sol_p], ["flux", "pressure"], "darcy")
vtk.Do()
