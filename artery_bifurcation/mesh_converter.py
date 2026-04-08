import vtk
import meshio
import numpy as np
from scipy.spatial import cKDTree
from pathlib import Path


def convert_to_msh(vtu_file, vtp_file, output_file):
    # Convert vtp to vtu using vtk
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(vtp_file)
    reader.Update()
    converter = vtk.vtkAppendFilter()
    converter.AddInputData(reader.GetOutput())
    converter.Update()
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName("_surface_tmp.vtu")
    writer.SetInputData(converter.GetOutput())
    writer.Write()

    # Read both meshes
    vol = meshio.read(vtu_file)
    surf = meshio.read("_surface_tmp.vtu")

    # Map surface points to volume points
    tree = cKDTree(vol.points)
    _, mapping = tree.query(surf.points)
    surf_triangles = mapping[surf.cells[0].data]
    face_ids = surf.cell_data["ModelFaceID"][0]
    unique_faces = np.unique(face_ids)

    # Write msh file manually
    with open(output_file, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")

        # Physical names
        f.write(f"$PhysicalNames\n{1 + len(unique_faces)}\n")
        f.write('3 1 "volume"\n')
        for fid in unique_faces:
            f.write(f'2 {fid} "bdry{fid}"\n')
        f.write("$EndPhysicalNames\n")

        # Nodes
        points = vol.points
        f.write(f"$Nodes\n{len(points)}\n")
        for i, p in enumerate(points):
            f.write(f"{i+1} {p[0]} {p[1]} {p[2]}\n")
        f.write("$EndNodes\n")

        # Elements
        n_tri = len(surf_triangles)
        n_tet = len(vol.cells[0].data)
        f.write(f"$Elements\n{n_tri + n_tet}\n")
        for i, (tri, fid) in enumerate(zip(surf_triangles, face_ids)):
            nodes = " ".join(str(n+1) for n in tri)
            f.write(f"{i+1} 2 2 {fid} {fid} {nodes}\n")
        for i, tet in enumerate(vol.cells[0].data):
            nodes = " ".join(str(n+1) for n in tet)
            f.write(f"{n_tri+i+1} 4 2 1 1 {nodes}\n")
        f.write("$EndElements\n")

    print(f"Written to {output_file}")
    print(f"Boundaries: {['bdry'+str(fid) for fid in unique_faces]}")


base_dir = Path(__file__).resolve().parent
vtu_file = str(base_dir / "ANR_106cgs.vtu")
vtp_file = str(base_dir / "ANR_106cgs.vtp")

convert_to_msh(mesh_file, mesh_file, "geometry.msh")
