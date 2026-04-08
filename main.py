from ngsolve import *
import meshio
import numpy as np

mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))

fes = L2(mesh, order=0)
gf = GridFunction(fes)
gf.Set(x + y)

# Extract mesh nodes and connectivity
points = np.array([mesh[v].point for v in mesh.vertices])
cells = np.array([[mesh[v].nr for v in el.vertices] for el in mesh.Elements()])

cell_values = np.array(gf.vec)

meshio.write("solution.vtu", meshio.Mesh(
    points=points,
    cells=[("triangle", cells)],
    cell_data={"u": [cell_values]}
))