from __future__ import annotations


import numpy as np
import scipy as sc
from typing import List, Tuple, Callable, Optional, Union


CONS_TRUST_RADIUS = 0.4
CONV_TOL = 1e-5
N_MACRO = 150
N_MICRO = 100
GLOBAL_RED_FACTOR = 0.05
LOCAL_RED_FACTOR = 0.005


def solver(
    func: Callable[[np.ndarray], float],
    grad: Callable[[np.ndarray], np.ndarray],
    hess_diag: Callable[[np.ndarray], np.ndarray],
    hess_x: Callable[[np.ndarray, np.ndarray], np.ndarray],
    norb: int,
    direction: str,
    line_search: str = "brent",
    line_search_orb_order: int = 4,
) -> np.ndarray:
    # initialize orbital rotation matrices
    u = np.eye(norb, dtype=np.float64)

    # initialize total number of Hessian linear transformations
    tot_hx = 0

    # initialize list for function evaluations per macro iteration
    func_its = []

    # starting trust radius
    trust_radius = CONS_TRUST_RADIUS

    # inititalize boolean for convergence of macro iterations
    macro_converged = False

    # calculate initial cost function
    f = func(u)

    # calculate initial gradient
    g = grad(u)

    # print header
    print(97 * "-")
    print(
        " Iteration | Cost function | Gradient norm | Level shift "
        "|   Micro    | Trust radius | Step size "
    )
    print(
        "           |               |               |             "
        "| iterations |              |           "
    )
    print(97 * "-")
    mu: float
    imicro: int
    n_kappa: float

    for imacro in range(N_MACRO):
        # add cost function
        func_its.append(f)

        # calculate gradient norm
        g_norm = np.linalg.norm(g)

        # calculate Hessian diagonal
        hdiag = hess_diag(u)

        # check for saddle point
        if g_norm / norb < CONV_TOL and min(hess_diag(u)) < -CONV_TOL:
            print("Reached saddle point. Applying small random perturbation.")

            # generate small, random rotation
            kappa = np.random.random((norb - 1) * norb // 2) * 1e-3
            u = sc.linalg.expm(unpack(kappa, norb))

            # recalculate gradient
            g = grad(u)

            # recalculate gradient norm
            g_norm = np.linalg.norm(g)

            # recalculate Hessian diagonal
            hdiag = hess_diag(u)

        # log results
        if imacro == 0:
            mu_str = f"{'-':^9}"
            imicro_str = f"{'-':>3}"
            trust_radius_str = f"{'-':^8}"
            stepsize_str = f"{'-':^8}"

        else:
            mu_str = f"{mu:>9.2e}"
            imicro_str = f"{imicro:>3d}"
            trust_radius_str = f"{trust_radius:>8.2e}"
            stepsize_str = f"{n_kappa * np.linalg.norm(kappa):>8.2e}"
        print(
            f"     {imacro:>3d}   |    {f:>8.2e}   |    {g_norm  / norb:>8.2e}   "
            f"|  {mu_str}  |      {imicro_str}   |    {trust_radius_str}  "
            f"|  {stepsize_str} "
        )

        # check for convergence
        if g_norm / norb < CONV_TOL:
            macro_converged = True
            break

        # generate trial vectors
        red_space_basis = [g / g_norm]
        trial = -g - hess_x(u, g) / g_norm
        red_space_basis.append(gram_schmidt(trial, red_space_basis))
        ntrial_start = 2
        trial = np.zeros_like(g, dtype=np.float64)
        min_idx = np.argmin(hdiag)
        if hdiag[min_idx] < 0.0:
            trial[min_idx] = 1.0
            red_space_basis.append(gram_schmidt(trial, red_space_basis))
            ntrial_start += 1

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
                elif imicro >= 5 and residual_norm > 0.9 * initial_residual_norm:
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

            # decrease trust radius if micro iterations are unable to converge
            if not micro_converged:
                aug_hess = aug_hess[: ntrial_start + 1, : ntrial_start + 1]
                red_space_basis = red_space_basis[:ntrial_start]
                h_basis = h_basis[:ntrial_start]
                basis_arr = basis_arr[:ntrial_start]
                initial_residual = True
                trust_radius /= 2
                tot_hx += len(red_space_basis)
                if trust_radius < 0.001:
                    print("Trust radius too small.")
                    if (
                        imacro > 4
                        and abs(func_its[-1] - func_its[-5]) < 1e-4 * func_its[-1]
                    ):
                        print("Cost function has converged.")
                        macro_converged = True
                    break
            else:
                trust_radius = min(2 * trust_radius, CONS_TRUST_RADIUS)

        # check if macroiterations have converged
        if macro_converged or not micro_converged:
            break

        # increment total number of Hessian linear transformations
        tot_hx += len(red_space_basis)

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
        grad_eval = lambda n_kappa: grad(
            u @ sc.linalg.expm(unpack(n_kappa * kappa, norb))
        )

        # get maximum frequency, period and upper bound for minimum
        omega_max = max(abs(np.linalg.eig(unpack(kappa, norb))[0]))
        period = 2 * np.pi / (line_search_orb_order * omega_max)
        upper_bound = period / 2

        # conduct a line search along the direction given by kappa
        if line_search == "golden":
            # perform golden section line search
            n_kappa, f, g = golden_search((0, upper_bound), func_eval, grad_eval)
        elif line_search == "brent":
            # perform Brent's line search
            n_kappa, f, g = brent_search((0, upper_bound), func_eval, grad_eval)
        elif line_search == "cubic":
            # perform cubic interpolation line search
            n_kappa, f, g = cubic_search(
                (0, upper_bound),
                func_eval,
                grad_eval,
                kappa,
                initial_f=f if direction == "min" else -f,
                initial_g=g,
            )
        else:
            raise ValueError("Unknown line search method.")

        # get correct function value
        if direction == "max":
            f = -f

        # get orbital rotation matrix
        u @= sc.linalg.expm(unpack(n_kappa * kappa, norb))

    if not macro_converged:
        raise RuntimeError("Orbital optimization has not converged!")

    if np.min(hdiag) < -CONV_TOL:
        raise RuntimeError("Orbital optimization has converged to saddle point!")

    print(97 * "-")
    print(f"Total number of Hessian linear transformations: {tot_hx}")

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
    red_space_basis: np.ndarray,
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
    ) -> Union[float, np.floating]:
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


def golden_search(
    brack: Tuple[float, float],
    func_eval: Callable[[float], float],
    grad_eval: Callable[[float], np.ndarray],
) -> Tuple[float, float, np.ndarray]:
    """
    this function performs a golden-section line search
    """
    # perform golden-section line search
    n_kappa, f, _ = sc.optimize.golden(
        func_eval, brack=brack, tol=1e-2, full_output=True
    )

    # calculate gradient
    g = grad_eval(n_kappa)

    return n_kappa, f, g


def brent_search(
    brack: Tuple[float, float],
    func_eval: Callable[[float], float],
    grad_eval: Callable[[float], np.ndarray],
) -> Tuple[float, float, np.ndarray]:
    """
    this function performs a Brent's line search
    """
    # perform golden-section line search
    n_kappa, f, _, _ = sc.optimize.brent(
        func_eval, brack=brack, tol=1e-2, full_output=True
    )

    # calculate gradient
    g = grad_eval(n_kappa)

    return n_kappa, f, g


def cubic_search(
    brack: Tuple[float, float],
    func_eval: Callable[[float], float],
    grad_eval: Callable[[float], np.ndarray],
    kappa: np.ndarray,
    initial_f: Optional[float] = None,
    initial_g: Optional[np.ndarray] = None,
) -> Tuple[float, float, np.ndarray]:
    """
    this function performs a cubic interpolation line search
    """
    # get function value at current point
    f1 = func_eval(brack[0]) if initial_f is None else initial_f

    # get gradient at current point
    g1 = grad_eval(brack[0]) if initial_g is None else initial_g

    # get directional derivative in search direction, factor 2 comes from multiplying
    # the unpacked matrices
    g1 = 2 * g1.T @ kappa

    # calculate cost function at upper bound
    f2 = func_eval(brack[1])

    # calculate gradient in search direction at upper bound
    g2 = 2 * grad_eval(brack[1]).T @ kappa

    # ensure lower and higher function values are assigned correctly
    if f1 < f2:
        n_kappa_low, f_low, g_low = brack[0], f1, g1
        n_kappa_high, f_high, g_high = brack[1], f2, g2
    else:
        n_kappa_low, f_low, g_low = brack[1], f2, g2
        n_kappa_high, f_high, g_high = brack[0], f1, g1

    while True:
        # perform cubic extrapolation
        d1 = g_low + g_high + 3 * (f_low - f_high) / (n_kappa_high - n_kappa_low)
        d2 = np.sign(n_kappa_high - n_kappa_low) * np.sqrt(d1**2 - g_low * g_high)
        n_kappa = n_kappa_high - (n_kappa_high - n_kappa_low) * (g_high + d2 - d1) / (
            g_high - g_low + 2 * d2
        )

        # calculate cost function at new point
        f = func_eval(n_kappa)

        # calculate gradient at new point
        g = grad_eval(n_kappa)

        # get directional derivative in search direction
        g_kappa = 2 * g.T @ kappa

        # new point is higher than lowest function evaluation,
        if f >= f_low:
            # set new point to higher function evaluation
            n_kappa_high, f_high, g_high = n_kappa, f, g_kappa
        # new point is lowest function evaluation
        else:
            # check convergence
            if abs(g_kappa) <= 1e-2 * -g1:
                break
            # set higher function evaluation to previous lowest function evaluation if
            # minimum is still bracketed
            if g_kappa * (n_kappa_high - n_kappa_low) >= 0:
                n_kappa_high, f_high, g_high = n_kappa_low, f_low, g_low
            # set new point to lowest function evaluation
            n_kappa_low, f_low, g_low = n_kappa, f, g_kappa

        # check if last step is valid
        if n_kappa_low == n_kappa_high:
            break

    return n_kappa, f, g
