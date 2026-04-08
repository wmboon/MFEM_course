
import ngsolve as ng
from netgen.read_gmsh import ReadGmsh
ngmesh = ReadGmsh("geometry.msh")
mesh = ng.Mesh(ngmesh)

print("Elements:", mesh.ne)
print("Boundaries:", mesh.GetBoundaries())

fes = ng.H1(mesh, order=1, dirichlet="wall2|wall3|wall4")
p0 = ng.L2(mesh, order=0)   # piecewise constant, fully discontinuous

sol = ng.GridFunction(fes)
cf = ng.IfPos(5.9-ng.z, 1, 0)

sol.Set(cf, definedon=mesh.Boundaries("wall2|wall3|wall4"))

u = fes.TrialFunction()
v = fes.TestFunction()

dx = ng.dx
f = ng.LinearForm(fes)
# f += v*dx

a = ng.BilinearForm(fes)
a += ng.grad(u)*ng.grad(v)*dx

a.Assemble()
f.Assemble()

r = f.vec - a.mat * sol.vec
sol.vec.data += a.mat.Inverse(fes.FreeDofs()) * r

vtk = ng.VTKOutput(mesh, coefs=[sol, -ng.grad(sol)], names=["sol", "flux"], filename="bifurc")
vtk.Do()