import ngsolve as ng
from netgen.read_gmsh import ReadGmsh
from pathlib import Path
import pyvista as pv
import numpy as np

# Import the grid
base_dir = Path(__file__).resolve().parent
mesh_file = str(base_dir / "data.msh")

ngmesh = ReadGmsh(mesh_file)
mesh = ng.Mesh(ngmesh)

vtk_file = str(base_dir / "data.vtk")
vtk_grid = pv.read(vtk_file)

# Read the subdomain data from the VTK file
subdomains = vtk_grid.cell_data["subdomains"]
subdomain_array = np.array(subdomains)

P0 = ng.L2(mesh, order=0)
source_p0 = ng.GridFunction(P0)

# Set the source term based on the subdomain array
source = subdomain_array == 3
source_p0.vec.data[:] = source.astype(float)

# Output some grid information
print("Number of elements:", mesh.ne)
print("Number of vertices:", mesh.nv)
print("Boundary tags:", mesh.GetBoundaries())

# Set up the finite element space
Lagrange1 = ng.H1(mesh, order=1, dirichlet=".*")

# The trial and test functions
u = Lagrange1.TrialFunction()
v = Lagrange1.TestFunction()

# Assemble the system matrix
a = ng.BilinearForm(Lagrange1)
a += ng.grad(u) * ng.grad(v) * ng.dx
a.Assemble()

# Preallocate the solution
sol = ng.GridFunction(Lagrange1)

# Assemble the source term (zero in this case)
f = ng.LinearForm(Lagrange1)
f += source_p0 * v * ng.dx
f.Assemble()

# Solve the problem and update the remaining degrees of freedom
rhs = f.vec - a.mat * sol.vec
sol.vec.data += a.mat.Inverse(Lagrange1.FreeDofs()) * rhs

# Output to vtk format for Paraview
filename = str(base_dir / "braaaaain")
vtk = ng.VTKOutput(
    mesh, coefs=[sol, -ng.grad(sol)], names=["pressure", "flux"], filename=filename
)
vtk.Do()
