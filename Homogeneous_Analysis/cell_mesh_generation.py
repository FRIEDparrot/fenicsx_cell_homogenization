import warnings 
import numpy as np
import dolfinx as dfx
from dolfinx.io import gmshio
from mpi4py import MPI
import ufl
import gmsh
from enum import Enum

class CellSectionType(Enum):
    SQUARE = 1
    CIRCLE = 2
    CUSTOM = 3

occ = gmsh.model.occ   # short write
geo = gmsh.model.geo   # short write

EPS = 1e-3

class StrutCell:
    def __init__(self,
                 points : np.ndarray,
                 grid: np.ndarray,
                 section: CellSectionType = CellSectionType.CIRCLE,
                 diameter: float = 0.2):
        """
        Initialize a StrutVoxel object.
        :param points:
        :param grid:
        :param diameter:
        """
        assert diameter > 0, "diameter should be positive"
        self.n = len(points)
        self.diameter = diameter
        self.points = points
        self.grid = grid
        self.section = section
        self.check_valid()

        gmsh.initialize()    # initialize gmsh
        gmsh.model.add("StrutCell")   # add a new model
        gmsh.option.setNumber("Geometry.Tolerance", EPS)
        gmsh.option.setNumber("Geometry.ToleranceBoolean", EPS)
        gmsh.option.setNumber("Geometry.OCCBoundsUseStl", 1) #  use STL for OCC bounds to compute more accurate bounding boxes

    def check_valid(self):
        """
        Make simple check to see if the grid is valid.
        :return:
        """
        assert 0 < self.diameter < 1, "diameter should be positive and less than 1"
        assert self.grid.dtype in (np.int64, np.int32), "Grid should be an array of integers"
        assert self.grid.ndim == 2 and self.grid.shape[1] == 2, "grid should be an array with each row containing two integers"
        assert self.points.dtype == np.float64 and self.points.ndim == 2 and self.points.shape[1] == 3
        all_connections = np.unique(self.grid.reshape(-1))
        assert len(self.points) == len(np.unique(self.points, axis=0)), \
            "some points are duplicated"
        assert np.all(np.sort(all_connections) == np.arange(self.n)), \
            "Invalid grid, some points are not connected or out of range"
        assert len(np.unique(self.grid, axis=0)) == len(self.grid), \
            "Invalid grid, some connections are duplicated"

    @staticmethod
    def locate_box_surface(dim: int, tag: int, p0: np.ndarray, p1: np.ndarray):
        assert dim == 3, "dim should be 3, but got {}".format(dim)
        assert p0.shape == (3,) and p1.shape == (3,), "x0 and x1 should be 3D vectors"

        surfaces = gmsh.model.getBoundary([(dim, tag)], oriented=False)
        if not surfaces:
            raise ValueError(f"No boundaries found for entity (dim={dim}, tag={tag})")
        axis_map = {
            0: ("left", "right"),
            1: ("front", "back"),
            2: ("bottom", "top")
        }
        faces = { "left": [], "right": [], "bottom": [], "top": [], "front": [], "back": [] }
        for dim_srf, tag_srf in surfaces:    # debug
            xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim_srf, tag_srf)
            coords = [(xmin, xmax), (ymin, ymax), (zmin, zmax)]
            for axis, (low_face, high_face) in axis_map.items():
                low, high = coords[axis]
                if abs(low - high) < EPS:
                    if abs(low - p0[axis]) < EPS:
                        faces[low_face].append((dim_srf, tag_srf))
                    elif abs(low - p1[axis]) < EPS:
                        faces[high_face].append((dim_srf, tag_srf))
                    break

        # check the symmetry and
        for axis, (low_face, high_face) in axis_map.items():
            if len(faces[low_face]) != len(faces[high_face]):
                raise ValueError(f"The cell is not symmetric, Number of {low_face} and {high_face} faces should be the same")
            faces[low_face].sort(key=lambda x: gmsh.model.getBoundingBox(*x)[axis])
            faces[high_face].sort(key=lambda x: gmsh.model.getBoundingBox(*x)[axis])
        return faces

    @staticmethod
    def build_affine_mat(dx: float, dy: float, dz: float):
        affine_mat = np.array(
            [[1.0,0,0,dx],
             [0,1.0,0,dy],
             [0,0,1.0,dz],
             [0,0,0, 1.0]],
            dtype=np.float64
        )
        return affine_mat.reshape(-1)

    def create(self,
               add_vertex: bool = True,
               vertex_diameter: float = None,
               mesh_size: float = None,
               filename: str = "voxel.msh"):
        """
        Build the geometry of the strut cell.
        :param add_vertex: if true, add a vertex sphere to points
        :param vertex_diameter: the diameter of the vortex sphere, default is diameter / 2
        :param mesh_size: the size of the mesh, default is mesh_size = diameter / 2
        :param filename: the name of the output file
        :return:
        """
        r = self.diameter / 2.0
        p0 = np.array([0, 0, 0], dtype=np.float64)
        p1 = np.array([1, 1, 1], dtype=np.float64)

        object_volumes = []  # list of (3, tag) to record the volume tags
        if mesh_size is None:
            mesh_size = r
        gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size) # max mesh size

        bnd_box = occ.addBox(*p0,*(p1-p0))  # add a unit cube as boundary box -> return tag 1
        occ.synchronize()
        # create main struct of the cell
        for i,j in self.grid :
            if self.section == CellSectionType.CIRCLE:
                start = self.points[i].astype(float)
                vec   = (self.points[j] - self.points[i]).astype(float)
                cyl_tag = occ.addCylinder(*start, *vec, r=r)
                object_volumes.append((3, cyl_tag))
            elif self.section == CellSectionType.SQUARE:
                raise NotImplementedError
            elif self.section == CellSectionType.CUSTOM:
                raise NotImplementedError

        if add_vertex :
            if vertex_diameter is None:
                vertex_diameter = self.diameter
            if add_vertex :
                if vertex_diameter is None:
                    vertex_diameter = self.diameter
                vertices_tag = [occ.addSphere(*p, vertex_diameter/2) for p in self.points]
                for v_tag in   vertices_tag:
                    object_volumes.append((3, v_tag))
        occ.synchronize()

        if len(object_volumes) > 0:
            fused, _ =  occ.fuse(object_volumes[0:1], object_volumes[1:], -1,  True, True)   # fuse all volumes
            occ.synchronize()
            # the needed part is part of the fused volume inside the box, not the box inside the fused volume
            vol_res, _ = occ.intersect(
                fused,
                [(3, bnd_box)],
                -1,
                True,
                True)
            vol_tag = vol_res[0][1]
        else:
            raise ValueError("No object volumes found")

        occ.removeAllDuplicates()  # remove duplicates
        occ.synchronize()

        # compute the outer faces of the cell
        try:
            # add physical groups after the last synchronization !!!
            sym_faces = self.locate_box_surface(dim=3, tag=vol_tag, p0=p0, p1=p1)
            gmsh.model.mesh.setPeriodic(2,
                                        [tag for (_, tag) in sym_faces["right"]], [tag for (_, tag) in sym_faces["left"]],
                                        self.build_affine_mat(1, 0, 0))
            gmsh.model.mesh.setPeriodic(2,
                                        [tag for (_, tag) in sym_faces["back"]], [tag for (_, tag) in sym_faces["front"]],
                                        self.build_affine_mat(0, 1, 0))
            gmsh.model.mesh.setPeriodic(2,
                                        [tag for (_, tag) in sym_faces["top"]], [tag for (_, tag) in sym_faces["bottom"]],
                                        self.build_affine_mat(0, 0, 1))
        except Exception as e:
            warnings.warn(f"Failed to generate periodic mesh for {filename}: {e}")
        occ.synchronize()  # synchronize again after setting periodic

        pg1 = gmsh.model.addPhysicalGroup(3, [vol_tag], -1, name="cell")  # add physical group to the fused volume
        all_faces = gmsh.model.getEntities(2)
        gmsh.model.addPhysicalGroup(2, [tag for (_, tag) in all_faces], -1, name="faces")
        all_edges = gmsh.model.getEntities(1)
        gmsh.model.addPhysicalGroup(1, [tag for (_, tag) in all_edges], -1, name="edges")
        all_vertices = gmsh.model.getEntities(0)
        gmsh.model.addPhysicalGroup(0, [tag for (_, tag) in all_vertices], -1, name="vertices")
        gmsh.model.mesh.setOrder(1)  # set the order of the mesh to 1
        gmsh.model.mesh.generate(3)
        gmsh.write(filename)
        gmsh.finalize()

# Below is partial code of E:/works/PythonFEMAnalysis/Homogeneous_Analysis/cell_mesh_generation.py:
def get_cell_density(filename) -> float:
    """
    read the .msh file and calculate the relative density
    :param filename:
    :return:
    """
    domain, _, _ = gmshio.read_from_msh(filename,comm=MPI.COMM_WORLD, gdim=3)
    # calculate the volume integral

    # calculate the volume
    vol = dfx.fem.assemble_scalar(dfx.fem.form(1.0 * ufl.dx(domain)))  # assemble the volume
    return vol / 1.0

# ===================== test part ======================

def main():
    from tests.cell_data import get_cube_data
    points, grid, diameter = get_cube_data()
    cell = StrutCell(points, grid, section=CellSectionType.CIRCLE, diameter=0.2)
    cell.create(add_vertex=True, vertex_diameter=0.25, mesh_size=None, filename="cube.msh")
    volume = get_cell_density("cube.msh")
    print("cell density : ", volume)

if __name__ == '__main__':
    main()
