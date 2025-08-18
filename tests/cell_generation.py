from Homogeneous_Analysis.cell_mesh_generation import StrutCell, CellSectionType
import os
from tests.cell_data import (
    get_cube_data, get_octahedron_data, get_octet_truss_data,
    get_tetrahedron_data, get_cross_lattice_data, get_star_cube_data,
    get_enhanced_voxel_data, get_tesseract_data, get_vintiles_data,
    get_face_centered_cubic_data, get_volume_center_data
)


os.chdir(os.path.abspath(os.path.dirname(__file__)))

MESH_SIZE = 0.03    # **the result will be incorrect either increase or decrease this value**

def test_cube():
    """
    Test the StrutCell class with a cube.
    """
    points, grid, diameter = get_cube_data()
    voxel = StrutCell(points, grid, section=CellSectionType.CIRCLE, diameter=diameter)
    voxel.check_valid()
    voxel.create(mesh_size=0.1, filename="cells/cube.msh")

def test_octahedron():
    """
    Test the StrutCell class with an octahedron.
    """
    points, grid, diameter = get_octahedron_data()
    voxel = StrutCell(points, grid, section=CellSectionType.CIRCLE, diameter=diameter)
    voxel.check_valid()
    voxel.create(mesh_size=0.1, filename="cells/octahedron.msh")

def test_octet_truss():
    """
    Test the StrutCell class with an octet truss.
    """
    points, grid, diameter = get_octet_truss_data()
    voxel = StrutCell(points, grid, section=CellSectionType.CIRCLE, diameter=diameter)
    voxel.check_valid()
    voxel.create(mesh_size=0.05, filename="cells/octet_truss.msh")

# asymmetric,  not test
def test_tetrahedron():
    """
    Test the StrutCell class with a tetrahedron.
    """
    points, grid, diameter = get_tetrahedron_data()
    voxel = StrutCell(points, grid, section=CellSectionType.CIRCLE, diameter=diameter)
    voxel.check_valid()
    voxel.create(mesh_size=0.05, filename="cells/tetrahedron.msh")

def test_cross_lattice():
    """
    Test the StrutCell class with a cross lattice.
    """
    points, grid, diameter = get_cross_lattice_data()
    voxel = StrutCell(points, grid, section=CellSectionType.CIRCLE, diameter=diameter)
    voxel.check_valid()
    voxel.create(mesh_size=0.1, filename="cells/cross_lattice.msh")

def test_star_cube():
    """
    Test the StrutCell class with a star cube.
    """
    points, grid, diameter = get_star_cube_data()
    voxel = StrutCell(points, grid, section=CellSectionType.CIRCLE, diameter=0.1)
    voxel.check_valid()
    voxel.create(mesh_size=0.1, filename="cells/star_cube.msh")

def test_enhanced_voxel():
    """
    Test the StrutCell class with an enhanced voxel.
    """
    points, grid, diameter = get_enhanced_voxel_data()
    voxel = StrutCell(points, grid, section=CellSectionType.CIRCLE, diameter=diameter)
    voxel.check_valid()
    voxel.create(mesh_size=0.05, filename="cells/enhanced_voxel.msh")

def test_tesseract_structure():
    """
    Test the StrutCell class with a tesseract structure.
    """
    points, grid, diameter = get_tesseract_data()
    voxel = StrutCell(points, grid, section=CellSectionType.CIRCLE, diameter=diameter)
    voxel.check_valid()
    voxel.create(mesh_size=0.03, filename="cells/tesseract.msh")

def test_vintiles_structure():
    """
    Test the StrutCell class with a vintiles structure.
    """
    points, grid, diameter = get_vintiles_data()
    voxel = StrutCell(points, grid, section=CellSectionType.CIRCLE, diameter=diameter)
    voxel.check_valid()
    voxel.create(mesh_size=0.03, filename="cells/vintiles.msh")

def test_volume_center():
    """
    Test the StrutCell class with a star structure.
    """
    points, grid, diameter = get_volume_center_data()
    voxel = StrutCell(points, grid, section=CellSectionType.CIRCLE, diameter=diameter)
    voxel.check_valid()
    voxel.create(mesh_size=0.1, filename="cells/volume_center.msh")

def test_face_centered_cubic():
    """
    Test the StrutCell class with a star structure.
    :return:
    """
    points, grid, diameter = get_face_centered_cubic_data()
    voxel = StrutCell(points, grid, section=CellSectionType.CIRCLE, diameter=diameter)
    voxel.check_valid()
    voxel.create(mesh_size=0.03, filename="cells/face_centered_cubic.msh")

if __name__ == "__main__":
    test_cube()
    test_octahedron()
    test_octet_truss()
    test_tetrahedron()
    test_cross_lattice()
    test_star_cube()
    test_enhanced_voxel()
    test_tesseract_structure()
    test_vintiles_structure()
    test_volume_center()
    test_face_centered_cubic()
