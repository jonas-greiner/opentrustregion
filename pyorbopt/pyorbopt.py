from __future__ import annotations


import numpy as np
import scipy as sc
from typing import List, Tuple, Callable
from pyscf import gto, lo


CONS_TRUST_RADIUS = 0.4
CONV_TOL = 1e-5
N_MACRO = 100
N_MICRO = 50
GLOBAL_RED_FACTOR = 0.05
LOCAL_RED_FACTOR = 0.005

INV_GOLDEN_RATIO = (np.sqrt(5) - 1) / 2
MAX_LINE_SEARCH_DEPTH = 25


def solver(
    func: Callable[[np.ndarray], float],
    grad: Callable[[np.ndarray], np.ndarray],
    hess_diag: Callable[[np.ndarray], np.ndarray],
    hess_x: Callable[[np.ndarray, np.ndarray], np.ndarray],
    norb: int,
    direction: str,
) -> np.ndarray:
    # initialize orbital rotation matrices
    u = np.eye(norb, dtype=np.float64)

    print("-------------------------------------------")
    print(" Iteration | Cost function | Gradient norm ")
    print("-------------------------------------------")

    for imacro in range(N_MACRO):
        # calculate cost function
        f = func(u)

        # calculate gradient
        g = grad(u)

        # calculate gradient norm
        g_norm = np.linalg.norm(g)

        # log results
        print(f"     {imacro:>3d}   |   {f:>1.2e}    |    {g_norm:>1.2e}   ")

        # check for convergence
        if g_norm < CONV_TOL:
            break

        # calculate Hessian diagonal
        hdiag = hess_diag(u)

        # generate trial vectors
        red_space_basis = [g / g_norm]
        trial = -g - hess_x(u, g) / g_norm
        red_space_basis.append(gram_schmidt(trial, red_space_basis))
        trial = np.zeros_like(g, dtype=np.float64)
        min_idx = np.argmin(hdiag)
        if hdiag[min_idx] < 0.0:
            trial[min_idx] = 1.0
            red_space_basis.append(gram_schmidt(trial, red_space_basis))

        # initialize boolean for the calculation of the initial residual
        initial_residual = True

        # calculate linear transformation of basis vectors
        h_basis = [hess_x(u, vec) for vec in red_space_basis]

        # get array of basis vectors
        basis_arr = np.vstack(red_space_basis)

        # construct augmented Hessian in reduced space
        aug_hess = np.zeros(2 * (len(red_space_basis) + 1,), dtype=np.float64)
        for j, vec in enumerate(h_basis, start=1):
            aug_hess[1:, j] = basis_arr @ vec

        # starting trust radius
        trust_radius = CONS_TRUST_RADIUS

        # decrease trust radius until micro iterations converge
        micro_converged = False
        while not micro_converged:
            # loop over micro iterations
            for imicro in range(N_MICRO):
                # perform bisection
                solution, red_solution, mu = bisection(
                    aug_hess, g, trust_radius, basis_arr
                )

                # calculate residual
                residual = g + red_solution.T @ np.vstack(h_basis) - mu * solution

                # get residual norm
                residual_norm = np.linalg.norm(residual)

                # save initial residual
                if initial_residual:
                    initial_residual_norm = residual_norm
                    initial_residual = False

                # determine reduction factor depending on whether local region is
                # reached
                red_factor = GLOBAL_RED_FACTOR if mu != 0.0 else LOCAL_RED_FACTOR

                # check if micro iterations have converged
                if (
                    residual_norm < red_factor * initial_residual_norm
                    or initial_residual_norm < 1e-10
                ):
                    micro_converged = True
                    break
                elif imicro >= 5 and residual_norm > 0.1 * initial_residual_norm:
                    break

                # precondition residual
                precond = hdiag - mu
                precond[abs(precond) < 1e-10] = 1e-10
                residual = np.divide(residual, precond)

                # orthogonalize to previous trial vectors
                red_space_basis.append(gram_schmidt(residual, red_space_basis))

                # add linear transformation of new basis vector
                h_basis.append(hess_x(u, red_space_basis[-1]))

                # get array of basis vectors
                basis_arr = np.vstack(red_space_basis)

                # construct new augmented Hessian
                new_aug_hess = np.zeros(
                    2 * (len(red_space_basis) + 1,), dtype=np.float64
                )
                new_aug_hess[:-1, :-1] = aug_hess
                new_aug_hess[1:, -1] = new_aug_hess[-1, 1:] = basis_arr @ h_basis[-1]
                aug_hess = new_aug_hess

            # decreast trust radius if micro iterations are unable to converge
            if not micro_converged:
                print("Micro iterations are unable to converge. Halving trust radius.")
                trust_radius /= 2

        # get new orbital transformation
        kappa = solution

        # prepare functions for line search
        if direction == "min":
            func_eval = lambda n_kappa: func(
                u @ sc.linalg.expm(unpack(n_kappa * kappa, norb))
            )
        elif direction == "max":
            func_eval = lambda n_kappa: -func(
                u @ sc.linalg.expm(unpack(n_kappa * kappa, norb))
            )

        # conduct a line search along the direction given by kappa
        n_kappa = sc.optimize.golden(func_eval, brack=(0.0, 1.0), tol=1e-3)

        # get orbital rotation matrix
        u @= sc.linalg.expm(unpack(n_kappa * kappa, norb))

    print("-------------------------------------------")

    if np.min(hdiag) < -CONV_TOL:
        raise RuntimeError("Orbital optimization has converged to saddle point!")

    return u


def gram_schmidt(vector: np.ndarray, space: List[np.ndarray]) -> np.ndarray:
    """
    this function orthonormalizes a vector with respect to a vector space
    """
    orth_vector = vector
    for space_vec in space:
        orth_vector -= (
            np.dot(vector, space_vec) / np.dot(space_vec, space_vec) * space_vec
        )

    return orth_vector / np.linalg.norm(orth_vector)


def bisection(
    aug_hess: np.ndarray,
    grad: np.ndarray,
    trust_radius: float,
    red_space_basis: List[np.ndarray],
) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    this function performs bisection to find the parameter alpha that matches the
    desired trust radius
    """

    def get_ah_lowest_eigenvec(
        aug_hess: np.ndarray,
        grad: np.ndarray,
        alpha: float,
        red_space_basis: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray, float]:
        """
        this function returns the lowest eigenvector for an augmented Hessian
        """
        # finish construction of augmented Hessian
        aug_hess[0, 1] = aug_hess[1, 0] = alpha * np.linalg.norm(grad)

        # diagonalize augmented Hessian
        eigenvalues, eigenvecs = np.linalg.eig(aug_hess)

        # get index of lowest eigenvalue
        idx = np.argmin(eigenvalues)

        # scale lowest eigenvector so that its first element is one
        lowest_eigenvalue = eigenvalues[idx]
        lowest_eigenvec = eigenvecs[:, idx] / eigenvecs[0, idx]

        # get solution vector
        solution = lowest_eigenvec[1:].T / alpha @ red_space_basis

        return solution, lowest_eigenvec[1:] / alpha, lowest_eigenvalue

    def get_ah_lowest_eigenvec_norm_trust_radius_diff(
        alpha: float,
        aug_hess: np.ndarray,
        grad: np.ndarray,
        red_space_basis: np.ndarray,
    ) -> float:
        """
        this function returns the norm of the lowest eigenvector for an augmented
        Hessian
        """
        solution, _, _ = get_ah_lowest_eigenvec(aug_hess, grad, alpha, red_space_basis)

        return np.linalg.norm(solution) - trust_radius

    # find root through bisection
    try:
        alpha = sc.optimize.bisect(
            get_ah_lowest_eigenvec_norm_trust_radius_diff,
            1e-4,
            1e6,
            args=(aug_hess, grad, red_space_basis),
        )
    # targeted trust radius is likely outside the bracketing range
    except ValueError:
        # solve Newton equations in reduced space without level shift
        red_grad = np.zeros(len(red_space_basis))
        red_grad[0] = np.linalg.norm(grad)
        red_solution = np.linalg.solve(aug_hess[1:, 1:], -red_grad)
        solution = red_solution.T @ red_space_basis
        return solution, red_solution, 0.0

    # return level-shifted solution
    return get_ah_lowest_eigenvec(aug_hess, grad, alpha, red_space_basis)


def unpack(vector: np.ndarray, norb: int):
    """
    this function unpacks a vector describing the unique parameters of an antihermitian
    matrix
    """
    matrix = np.zeros(2 * (norb,))
    idx = np.tril_indices(norb, -1)
    matrix[idx] = vector

    return matrix - matrix.conj().T


if __name__ == "__main__":
    # prepare molecule
    mol = gto.M(
        verbose=1,
        atom="""
        O	 0.0000000	 0.0000000	 0.0000000
        H	 0.7569685	 0.0000000	-0.5858752
        H	-0.7569685	 0.0000000	-0.5858752
        """,
        unit="Angstrom",
        basis="cc-pvdz",
    )

    # run HF calculation
    mf = mol.RHF().run()

    # canonical initial guess
    mo_coeff = mf.mo_coeff[:, : min(mol.nelec)]

    # atomic initial guess
    # u = lo.boys.atomic_init_guess(mol, mo_coeff)
    # mo_coeff = mo_coeff @ u

    # generate random shift
    # dr = np.cos(np.arange((norb-1)*norb//2)) * 1e-3
    # idx = np.tril_indices(norb, -1)
    # mat = np.zeros((norb,norb))
    # mat[idx] = dr
    # u = sc.linalg.expm(mat - mat.conj().T)
    # mo_coeff = mo_coeff @ u

    # initialize pyscf localization object
    fb = lo.Boys(mol, mo_coeff)

    # cost function
    func = fb.cost_function

    # gradient function
    def grad(u):
        return fb.get_grad(u)

    # hessian diagonal function
    def hess_diag(u):
        return fb.gen_g_hop(u)[2]

    # hessian linear transformation function
    def hess_x(u, x):
        return fb.gen_g_hop(u)[1](x)

    # decide whether to maximize or minimize
    direction = "max" if isinstance(fb, lo.PM) else "min"

    # call solver
    u = solver(func, grad, hess_diag, hess_x, min(mol.nelec), direction)

    # call pyscf solver
    fb.verbose = 4
    loc_orb = fb.kernel(mo_coeff)
