import pyvista as pv
import meshio
from pathlib import Path

# 1. Load the tet mesh
base_dir = Path(__file__).resolve().parent
vtk_file = str(base_dir / "data.vtk")
grid = pv.read(vtk_file)

# 2. Extract the boundary surface, keeping track of original point IDs
surface = grid.extract_surface(pass_pointid=True, pass_cellid=True, algorithm=None)

# 3. Pull out the triangle connectivity (surface.faces is VTK's flat format: [3, i0,i1,i2, 3, i0,i1,i2, ...])
faces = surface.faces.reshape(-1, 4)[:, 1:4]

# 4. Map local surface point indices back to original volume mesh point indices
orig_ids = surface["vtkOriginalPointIds"]
boundary_tris = orig_ids[faces]

# 5. Get the original tet cells
tet_mesh = meshio.read(vtk_file)
tet_cells = tet_mesh.cells_dict["tetra"]

# 6. Build a combined mesh: volume tets + boundary triangles
combined = meshio.Mesh(
    points=tet_mesh.points,
    cells=[("tetra", tet_cells), ("triangle", boundary_tris)],
)

meshio.write(
    str(base_dir / "data.msh"),
    combined,
    file_format="gmsh22",
    binary=False,
)
