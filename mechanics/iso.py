import numpy as np
import ufl

def get_voigt_pairs(dim):
    return [(0, 0), (1, 1), (0, 1)] if dim == 2 else \
        [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]

def lame_parameters(E, nu):
    # Material properties
    lmbda = (E * nu) / ((1 + nu) * (1 - 2 * nu))  # Lame's first parameter
    mu = E / (2 * (1 + nu))  # Lame's second parameter
    return lmbda, mu

def constitutive_mat(lmbda, mu, dim, return_numpy=False):
    """
    Compute the constitutive matrix

    Parameters:
    dim (int): Dimension (2 or 3)
    lmbda (float): Lamé's first parameter
    mu (float): Shear modulus (Lamé's second parameter)

    Returns:
    C4 (ndarray): 4th order constitutive tensor
    C_voigt (ndarray): Voigt notation matrix representation
    """
    assert dim in (2, 3), "Dimension must be 2 or 3"

    lam = float(lmbda)
    mu = float(mu)

    # Dimension-dependent parameters
    voigt_pairs = get_voigt_pairs(dim)
    n_voigt = len(voigt_pairs)
    delta = np.eye(dim, dtype=float)

    # Initialize tensors
    C4 = np.zeros((dim, dim, dim, dim), dtype=float)
    C_voigt = np.zeros((n_voigt, n_voigt), dtype=float)

    # Compute 4th order tensor
    for i in range(dim):
        for j in range(dim):
            for k in range(dim):
                for l in range(dim):
                    C4[i, j, k, l] = lam * delta[i, j] * delta[k, l] \
                                     + mu * (delta[i, k] * delta[j, l] + delta[i, l] * delta[j, k])

    # Convert to Voigt notation
    for I, (p, q) in enumerate(voigt_pairs):
        for J, (r, s) in enumerate(voigt_pairs):
            C_voigt[I, J] = C4[p, q, r, s]
    if return_numpy:
        return C4, C_voigt
    return ufl.as_tensor(C4), ufl.as_tensor(C_voigt)

def epsilon(u):
    """
    Small strain tensor
    :param u: ufl function or vector
    """
    return ufl.sym(ufl.grad(u))

def epsilon_voigt(u):
    """Return engineering strain in Voigt notation (vector)"""
    eps = epsilon(u)
    dim = u.ufl_shape[0]
    if dim == 2:
        # (xx, yy, 2*xy)
        idx_pairs = [(0, 0), (1, 1), (0, 1)]
        factors = [1, 1, 2]
    elif dim == 3:
        # (xx, yy, zz, 2*yz, 2*xz, 2*xy)
        idx_pairs = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
        factors = [1, 1, 1, 2, 2, 2]
    else:
        raise ValueError("Dimension must be 2 or 3")
    return ufl.as_vector([f * eps[i, j] for (i, j), f in zip(idx_pairs, factors)])

def sigma(eps, lmbda, mu):
    r"""
        Compute stress tensor for isotropic linear elasticity.
        Properly handles both normal and shear components.
        .. math::
            \sigma = \lambda \mathrm{tr}(\varepsilon)I + 2\mu\varepsilon

        :param eps: Strain tensor (mathematical form, ε_ij for shear, NOT engineering γ_ij)
        :param lmbda: First Lamé parameter
        :param mu: Shear modulus
        :param dim: Dimension (2 or 3)
        :return: Stress tensor
        """
    # eps is mathematical strain (ε_ij), so shear components are NOT pre-multiplied by 2
    dim = eps.ufl_shape[0]
    return lmbda * ufl.tr(eps) * ufl.Identity(dim) + 2 * mu * eps

def sigma_voigt(eps_voigt, lmbda, mu, dim=3):
    r"""
    Calculation for stress tensor
    .. math::
        σ =  \lambda  (\nabla \cdot u ) I  + 2 \mu \varepsilon
    see https://en.wikipedia.org/wiki/Hooke%27s_law
    :param eps_voigt: strain tensor in Voigt notation (1x6 vector)
    :param lmbda:
    :param mu:
    :param dim: dimension of the problem
    :return:
    """
    if dim == 2:
        C = np.array([
            [2 * mu + lmbda, lmbda, 0],  # C11, C12, C13, C14
            [lmbda, 2 * mu + lmbda, 0],  # C21, C22, C23, C24
            [0, 0, mu]
        ])  # C31, C32, C33, C34
    else :
        C = np.array([
            [2 * mu + lmbda, lmbda, lmbda, 0, 0, 0],
            [lmbda, 2 * mu + lmbda, lmbda, 0, 0, 0],
            [lmbda, lmbda, 2 * mu + lmbda, 0, 0, 0],
            [0, 0, 0, mu, 0, 0],
            [0, 0, 0, 0, mu, 0],
            [0, 0, 0, 0, 0, mu]
        ], dtype=float)
    D = ufl.as_tensor(C)
    return  ufl.dot(D,eps_voigt)


def unit_strain_tensor(dim=3, return_numpy=False) -> np.ndarray:
    """
    Returns the unit strain tensor for a given index.
    :param dim: dimension of the problem
    :param return_numpy: if True, returns a numpy array, else returns a UFL tensor
    :return:
    """
    assert dim in (2, 3), "Dimension must be 2 or 3"
    load_case_num = 3 if dim == 2 else 6   # number of load cases
    unit_tensors = np.zeros((load_case_num, dim, dim), dtype=float)
    voigt_pairs = get_voigt_pairs(dim) # voigt notation pairs
    for I in range(load_case_num):
        p, q = voigt_pairs[I]
        if p == q:
            # Normal strain
            unit_tensors[I][p][q] = 1.0
        else:
            # Shear strain -> symmetric, engineering shear = 1
            unit_tensors[I][p][q] = 0.5
            unit_tensors[I][q][p] = 0.5
    if return_numpy:
        return unit_tensors
    return ufl.as_tensor(unit_tensors)


def unit_strain_tensor_voigt(dim=3, return_numpy=False):
    """
    Returns the unit strain tensors in Voigt notation (1x6 vectors)

    [I][eps_x, eps_y, eps_z, 2eps_xy, 2eps_yz, 2eps_xz], where I is the voigt load direction
    :param dim: dimension of the problem
    :param return_numpy: if True, returns a numpy array, else returns a UFL tensor
    :return: array of shape (load_case_num, voigt_dim)
    """
    assert dim in (2, 3), "Dimension must be 2 or 3"
    load_case_num = 3 if dim == 2 else 6  # number of load cases
    voigt_dim = 3 if dim == 2 else 6  # dimension in Voigt notation

    unit_vectors = np.zeros((load_case_num, voigt_dim), dtype=float)

    for I in range(load_case_num):
        unit_vectors[I][I] = 1.0  # All unit vectors (engineering shear strains already include factor of 2)
    if return_numpy:
        return unit_vectors
    return ufl.as_tensor(unit_vectors)  # or convert to appropriate UFL format if needed


def main():
    E, nu = 200e9, 0.3
    lmbda, mu = lame_parameters(E, nu)
    # print(f"Lame's first parameter: {lmbda/1e9:.2f} GPa")  # 115.3 GPa
    # print(f"Lame's second parameter: {mu/1e9:.2f} GPa")    # 76.9 GPa
    tsr = unit_strain_tensor_voigt(dim=3)
    print(tsr)
    stress = sigma_voigt(tsr[3], lmbda, mu, dim=3)
    print(stress)

if __name__ == '__main__':
    main()

