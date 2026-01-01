import warnings
import pathlib
from Homogeneous_Analysis import (
    CellHomogenization,
    HomogenizationVisualizer
)
import os
import pickle
import numpy as np

os.chdir(os.path.dirname(os.path.abspath(__file__)))

voxel_base_path = pathlib.Path("cells")
voxel_img_path = pathlib.Path("../img")
all_voxels = [
    "cross_lattice.msh",
    "cube.msh",
    "enhanced_voxel.msh",
    "face_centered_cubic.msh",
    "octahedron.msh",
    "octet_truss.msh",
    "star_cube.msh",
    "tesseract.msh",
    "tetrahedron.msh",
    "vintiles.msh",
    "volume_center.msh"
]

def main():
    # Define material properties -> structural steel
    E = 200e9  # Young's modulus in Pa
    nu = 0.3  # Poisson's ratio

    results_txt = voxel_img_path / "C_H_results.txt"
    results_pickle = voxel_img_path / "C_H_results.pkl"

    # initialize the file text to blank
    with open(results_txt, "w") as f:
        f.write("# Homogenization Results\n\n")

    all_C_H = {}  # initialize dictionary to store C_H of each cell

    for voxel in all_voxels:
        try:
            cell_name = voxel.split(".")[0]  # remove the extension
            full_path = voxel_base_path / voxel

            print("solving cell homogenization :")
            homo = CellHomogenization(full_path, E, nu)
            C_H = homo.solve_homogenization()
            all_C_H[cell_name] = C_H

            with open(results_txt, "a") as f:
                f.write(f"Cell: {cell_name}\n")
                np.savetxt(f, C_H, fmt="%.6e")
                f.write("\n\n")

            visualizer = HomogenizationVisualizer(C_H, cell_name=cell_name)
            visualizer.plot_voigt_matrix(save_path=str(voxel_img_path / f"{cell_name}_voigt.png"))
            visualizer.plot_directional_young_modulus_3d(resolution=100, save_path=str(voxel_img_path / f"{cell_name}_directional_young_modulus_3d.png"))
            visualizer.plot_directional_young_modulus_2d(plane='xy', save_path=str(voxel_img_path / f"{cell_name}_directional_young_modulus_xy.png"))

            with open(results_pickle, "wb") as f:
                pickle.dump(all_C_H, f)

        except Exception as e:
            warnings.warn(f"Error processing {voxel}: {e}")

    with open(results_pickle, "wb") as f:
        pickle.dump(all_C_H, f)

if __name__ == '__main__':
    main()
