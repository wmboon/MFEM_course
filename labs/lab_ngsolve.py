import ngsolve as ng
import numpy as np
import scipy.sparse as sps

mesh = ng.Mesh(ng.unit_cube.GenerateMesh(maxh=0.1))
mesh.GetBoundaries

V = ng.HDiv(mesh, order=0, RT=True, dirichlet="back|front")
P = ng.L2(mesh, order=0)

fe_space = V * P

u, p = fe_space.TrialFunction()
v, q = fe_space.TestFunction()

A = ng.BilinearForm(fe_space)
A += u * v * ng.dx
A += - p * ng.div(v) * ng.dx
A += ng.div(u) * q * ng.dx
A.Assemble()

normal = ng.specialcf.normal(mesh.dim)
bc_p = ng.sin(ng.x * 2 * ng.pi) + ng.sin(ng.y * 2 *
                                         ng.pi) + ng.sin(ng.z * 2 * ng.pi)

g = ng.LinearForm(fe_space)
g += - bc_p * (v.Trace() * normal) * ng.ds
g.Assemble()

g_np = np.array(g.vec)

sol = ng.GridFunction(fe_space)

res = g.vec.data - A.mat * sol.vec
sol.vec.data += A.mat.Inverse(freedofs=fe_space.FreeDofs()) * res

flux, pressure = sol.components

# zero = ngs.CoefficientFunction(0.0)
# flux_3d = ngs.CoefficientFunction((flux[0], flux[1], zero))

vtk = ng.VTKOutput(mesh, [flux, pressure], ["flux", "pressure"], "darcy")
vtk.Do()

pass
