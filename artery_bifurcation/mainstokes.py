from pathlib import Path

import ngsolve as ng
from netgen.read_gmsh import ReadGmsh

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
Velocity = ng.VectorH1(mesh, order=1, dirichlet="bdry1")
Velocity.SetOrder(ng.TET, 4)
Velocity.Update()
Pressure = ng.H1(mesh, order=1)
fespace = Velocity * Pressure
print("ndof:", fespace.ndof)

# The trial and test functions
u, p = fespace.TrialFunction()
v, q = fespace.TestFunction()

dt = 0.1

# Assemble the system matrix
a = ng.BilinearForm(fespace)
a += ng.InnerProduct(ng.Sym(ng.grad(u)), ng.Sym(ng.grad(v))) * ng.dx
a += -p * ng.div(v) * ng.dx
a += q * ng.div(u) * ng.dx
a += u * v / dt * ng.dx
a.Assemble()


# Output to vtk format for Paraview
sol = ng.GridFunction(fespace)
velocity, pressure = sol.components

filename = str(base_dir / "bifurc_stokes")
vtk = ng.VTKOutput(
    mesh, coefs=[pressure, velocity], names=["pressure", "flux"], filename=filename
)
# vtk.Do(time=0)
# Preallocate u_prev
u_prev = ng.GridFunction(Velocity)
print(type(velocity))
for t in range(1, 30):
    # Set the natural boundary conditions
    normal = ng.specialcf.normal(3)
    f = ng.LinearForm(fespace)
    f += -v * normal * (1 + ng.sin(t * dt)) * ng.ds(definedon="bdry4")
    f += u_prev * v / dt * ng.dx
    f.Assemble()

    # Solve the problem and update the remaining degrees of freedom
    rhs = f.vec - a.mat * sol.vec
    sol.vec.data += a.mat.Inverse(fespace.FreeDofs(), inverse="pardiso") * rhs

    vtk.Do(time=t)
    u_prev = velocity
    print(t)
