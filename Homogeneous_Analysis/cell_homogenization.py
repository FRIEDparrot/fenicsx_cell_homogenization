import numpy as np
import dolfinx as dfx
import ufl
from dolfinx.io import gmshio
from mpi4py import MPI
from dolfinx.fem import functionspace, locate_dofs_geometrical
from mechanics.iso import (
    lame_parameters, sigma_voigt, constitutive_mat,
    unit_strain_tensor_voigt, epsilon_voigt
)
from dolfinx_mpc import MultiPointConstraint, LinearProblem as MPCProblem
from .cell_visual_toolkit import HomogenizationVisualizer
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))

class CellHomogenization:
    def __init__(self, filename, E, nu, dim=3):
        """
        Reading .msh file and making homogenization of cells
        """
        self.E = E
        self.nu = nu
        self.dim = dim
        self.lmbda, self.mu = lame_parameters(E, nu)
        self.C4, self.C_voigt = constitutive_mat(self.lmbda, self.mu, self.dim)

        msh, _, _ = gmshio.read_from_msh(filename, MPI.COMM_WORLD, gdim=3)
        self.domain = msh
        self.vx = 1.0
        self.vy = 1.0
        self.vz = 1.0
        self.vol = self.vx * self.vy * self.vz  # Volume of the domain

    @staticmethod
    def get_periodic_dofs(V, direction = 0):
        """
        Get the dofs of the periodic boundary conditions
        :param V:
        :param direction:
        :return:
        """
        dofs_slave = locate_dofs_geometrical(V, marker=lambda x: np.isclose(x[direction], 0.0))  # locate master dofs
        dofs_master = locate_dofs_geometrical(V, marker=lambda x: np.isclose(x[direction], 1.0))  # locate slave dofs
        return dofs_master, dofs_slave

    @staticmethod
    def periodic_indicator(x, direction=0):
        """
        :param x:
        :param direction:
        :return:
        """
        return np.isclose(x[direction], 1.0)

    @staticmethod
    def periodic_relation(x, direction=0):
        """
        Check if the mapping is close to the identity mapping
        :return:
        """
        out_x = x.copy()
        out_x[direction] = x[direction] - 1 # from slave surface to master surface
        return out_x

    def apply_periodic_bc(self, V, bcs):
        # Define boundary conditions -> periodic boundary conditions
        mpc = MultiPointConstraint(V = V)  # Define MPC object
        for direction in range(self.dim):
            mpc.create_periodic_constraint_geometrical(
                V,
                lambda x, d=direction:self.periodic_indicator(x, d),
                lambda x, d=direction:self.periodic_relation(x, d),
                bcs
            )
        mpc.finalize()  # finalize creation of MPC constraints
        return mpc

    @staticmethod
    def close_to_center(x):
        """
        :param x:
        :return:
        """
        return np.all(np.isclose(x, 0.5, atol=1e-6),  axis=0)

    def solve_homogenization(self):
        """
        1. Solve the displacement under each load case
        2. Apply the periodic boundary conditions
        3. Calculate the effective stiffness tensor C_H
        :return:
        """
        # Define function space
        V = functionspace(self.domain, ("CG", 2, (self.dim,)))

        # voigt index mapping
        u_solu = []
        voigt_index = np.array([[0,0], [1,1], [2,2], [1,2], [0,2], [0,1]])
        l = len(voigt_index)
        for idx in range(l):
            uh = ufl.TrialFunction(V)  # Displacement Function u^{ij}(x,y,z)
            v = ufl.TestFunction(V)   # Test function v(x,y,z)
            # get the unit strain tensor eps0_ij
            eps0_tensor = unit_strain_tensor_voigt(self.dim, return_numpy=True)[idx]  # unit strain tensor
            eps0 = dfx.fem.Constant(self.domain, eps0_tensor) # Constant strain tensor

            # Define variational problem
            a = ufl.inner(sigma_voigt(epsilon_voigt(uh), self.lmbda, self.mu, dim=self.dim), epsilon_voigt(v)) * ufl.dx(self.domain)
            L = ufl.inner(sigma_voigt(eps0, self.lmbda, self.mu, dim=self.dim), epsilon_voigt(v)) * ufl.dx(self.domain)

            # Dirichlet boundary condition
            dofs = dfx.fem.locate_dofs_geometrical(V, marker=self.close_to_center)
            bc = dfx.fem.dirichletbc(value=np.zeros(3), dofs=dofs, V=V)
            bcs = [bc]

            mpc = self.apply_periodic_bc(V, bcs)
            PETSC_OPTIONS = {"ksp_type": "cg", "pc_type": "hypre"}
            problem = MPCProblem(a, L, mpc, bcs, petsc_options=PETSC_OPTIONS)

            uh = problem.solve() # Solve problem
            u_solu.append(uh)    # get 6 u solutions at each load case

        C_H = np.zeros((l, l))  # Effective stiffness tensor
        C = self.C_voigt  # Constitutive tensor in Voigt notation
        # do integration over the domain to get the effective stiffness tensor

        for I in range(l):
            eps0_I = unit_strain_tensor_voigt(self.dim)[I]  # 1x6 strain tensor
            u_I =  u_solu[I]
            for J in range(I, l):
                eps0_J = unit_strain_tensor_voigt(self.dim)[J]  # 1x6 strain tensor
                u_J = u_solu[J]  # u^{ij} (x,y,z) and u^{kl} (x,y,z)

                # calculate the voigt strain vector
                eps_I_voigt = eps0_I - epsilon_voigt(u_I)   # ij,  I is voigt notation of load case
                eps_J_voigt = eps0_J - epsilon_voigt(u_J)   # kl,  J is voigt notation of load case

                # calculate the sum of central term
                inner_term = ufl.dot(eps_I_voigt, ufl.dot(C, eps_J_voigt))
                r"""
                Effective stiffness tensor
                C_{ijkl}^{H} = \frac{1}{|Y|} \int_{Y} C_{pqrs} (
                    \varepsilon_{pq}^{0(ij)} -\varepsilon_{pq}^{(ij)}) (\varepsilon_{rs}^{0(kl)} - \varepsilon_{rs}^{(kl)}
                ) d\Omega 
                """
                C_H[I, J] = dfx.fem.assemble_scalar(dfx.fem.form(1.0 * inner_term * ufl.dx(self.domain))) / self.vol
                C_H[J, I] = C_H[I, J]  # Symmetric matrix
        return C_H

def main():
    # only use cube cell for test
    filename  = "../voxels/cube.msh"
    E = 200e9 # Young's modulus in Pa
    nu = 0.3  # Poisson's ratio
    homo = CellHomogenization(filename, E, nu)
    C_H = homo.solve_homogenization()
    print(C_H)  # Effective stiffness tensor (in GPa)

    visualizer = HomogenizationVisualizer(C_H)
    visualizer.plot_directional_young_modulus_3d(resolution=100)

if __name__ == "__main__":
    main()
