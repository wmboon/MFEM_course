import ngsolve as ngs
import numpy as np
import scipy.sparse as sps

mesh = ngs.Mesh(ngs.unit_cube.GenerateMesh(maxh=0.1))

V = ngs.HDiv(mesh, order=0, RT=True, dirichlet="back|front")
P = ngs.L2(mesh, order=0)

fe_space = V * P

u, p = fe_space.TrialFunction()
v, q = fe_space.TestFunction()

A = ngs.BilinearForm(fe_space)
A += u * v * ngs.dx
A += -p * ngs.div(v) * ngs.dx
A += ngs.div(u) * q * ngs.dx
A.Assemble()

normal = ngs.specialcf.normal(mesh.dim)
bc_p = (
    ngs.sin(ngs.x * 2 * ngs.pi)
    + ngs.sin(ngs.y * 2 * ngs.pi)
    + ngs.sin(ngs.z * 2 * ngs.pi)
)

g = ngs.LinearForm(fe_space)
g += -bc_p * (v.Trace() * normal) * ngs.ds
g.Assemble()

# g_np = np.array(g.vec)

sol = ngs.GridFunction(fe_space)
res = g.vec.data - A.mat * sol.vec
sol.vec.data = A.mat.Inverse(freedofs=fe_space.FreeDofs()) * res

flux, pressure = sol.components

# zero = ngs.CoefficientFunction(0.0)
# flux3d = ngs.CoefficientFunction((flux[0], flux[1], flux[2]))

vtk = ngs.VTKOutput(mesh, [flux, pressure], ["flux", "pressure"], "darcy")
vtk.Do()

pass
