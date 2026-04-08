
import ngsolve as ng
from netgen.read_gmsh import ReadGmsh
from pathlib import Path

# Import the grid
base_dir = Path(__file__).resolve().parent
mesh_file = str(base_dir / "geometry.msh")
ngmesh = ReadGmsh(mesh_file)
mesh = ng.Mesh(ngmesh)

# Output some grid information
print("Number of elements:", mesh.ne)
print("Number of vertices:", mesh.nv)
print("Boundary tags:", mesh.GetBoundaries())

# Set up the finite element space
Lagrange1 = ng.H1(mesh, order=1, dirichlet="bdry2|bdry3|bdry4")

# The trial and test functions
u = Lagrange1.TrialFunction()
v = Lagrange1.TestFunction()

# Assemble the system matrix
a = ng.BilinearForm(Lagrange1)
a += ng.grad(u)*ng.grad(v)*ng.dx
a.Assemble()

# Preallocate the solution
sol = ng.GridFunction(Lagrange1)

# Set the essential boundary conditions
coef_func = ng.IfPos(5.9-ng.z, 1, 0)
sol.Set(coef_func, definedon=mesh.Boundaries("bdry2|bdry3|bdry4"))

# Assemble the source term (zero in this case)
f = ng.LinearForm(Lagrange1)
f.Assemble()

# Solve the problem and update the remaining degrees of freedom
rhs = f.vec - a.mat * sol.vec
sol.vec.data += a.mat.Inverse(Lagrange1.FreeDofs()) * rhs

# Output to vtk format for Paraview
vtk = ng.VTKOutput(mesh, coefs=[sol, -ng.grad(sol)],
                   names=["pressure", "flux"], filename="bifurc")
vtk.Do()
