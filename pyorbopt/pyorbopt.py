from __future__ import annotations


import logging
import numpy as np
import scipy as sc
from typing import TYPE_CHECKING, cast

if TYPE_CHECKING:
    from typing import List, Tuple, Callable, Optional, Union

    any_float = Union[float, np.floating]


CONS_TRUST_RADIUS = 0.4
N_MACRO = 150
N_MICRO = 100
GLOBAL_RED_FACTOR = 0.05
LOCAL_RED_FACTOR = 0.005
N_STAB_DAVIDSON = 100
SEED = 42


rng = np.random.default_rng(seed=SEED)


class OrbOpt:
    def __init__(
        self,
        conv_tol: float = 1e-5,
        trial_vectors: Union[str, int] = "orig",
        verbose: int = 1,
    ):
        # initialize total number of Hessian linear transformations
        self.tot_hx: int = 0
        # initialize list for function evaluations per macro iteration
        self.func_its: List[float] = []
        # initialize list for gradient RMS per macro iteration
        self.g_rms_its: List[any_float] = []
        # initialize convergence tolerance
        self.conv_tol = conv_tol
        # initialize trial vector choice
        self.trial_vectors = trial_vectors
        # initialize verbosity setting
        self.logger = _create_logger(verbose)

    def solver(
        self,
        unpack: Callable[[np.ndarray], np.ndarray],
        func: Callable[[np.ndarray], float],
        func_grad: Callable[[np.ndarray], Tuple[float, np.ndarray]],
        func_grad_hdiag_hess_x: Callable[
            [np.ndarray],
            Tuple[float, np.ndarray, np.ndarray, Callable[[np.ndarray], np.ndarray]],
        ],
        n_param: int,
        direction: str,
        u: Optional[np.ndarray] = None,
        line_search: str = "brent",
        line_search_orb_order: int = 4,
    ) -> np.ndarray:
        # initialize orbital rotation matrices
        if u is None:
            kappa = np.zeros(n_param, dtype=np.float64)
            u = cast(np.ndarray, sc.linalg.expm(unpack(kappa)))

        # starting trust radius
        trust_radius = CONS_TRUST_RADIUS

        # inititalize boolean for convergence of macro iterations
        macro_converged = False

        # print header
        self.logger.info(101 * "-")
        self.logger.info(
            " Iteration | Objective function | Gradient RMS | Level shift |   Micro    "
            "| Trust radius | Step size "
        )
        self.logger.info(
            "           |                    |              |             | iterations "
            "|              |           "
        )
        self.logger.info(101 * "-")
        mu: float
        imicro: int
        n_kappa: float

        for imacro in range(N_MACRO):
            # calculate cost function, gradient and Hessian diagonal
            f, g, hdiag, hess_x = func_grad_hdiag_hess_x(u)

            # add cost function
            self.func_its.append(float(f))
            self.logger.info(self.func_its[-1])

            # calculate gradient norm
            g_norm = np.linalg.norm(g)

            # RMS gradient
            self.g_rms_its.append(float(g_norm / np.sqrt(n_param)))

            # check for saddle point
            if self.g_rms_its[-1] < self.conv_tol and min(hdiag) < -self.conv_tol:
                self.logger.warning(
                    "Reached saddle point. Applying small random perturbation."
                )

                # generate small, random rotation
                kappa = rng.random(size=n_param) * 1e-3
                u += sc.linalg.expm(unpack(kappa))

                # recalculate cost function, gradient and Hessian diagonal
                f, g, hdiag, hess_x = func_grad_hdiag_hess_x(u)
                self.func_its[-1] = f

                # recalculate gradient norm
                g_norm = np.linalg.norm(g)

                # recalculate RMS gradient
                self.g_rms_its[-1] = g_norm / np.sqrt(n_param)

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
            self.logger.info(
                f"     {imacro:>3d}   |    {f:> 8.6e}   "
                f"|   {self.g_rms_its[-1]:>8.2e}   |  {mu_str}  |      {imicro_str}   "
                f"|    {trust_radius_str}  |  {stepsize_str} "
            )

            # check for convergence
            if self.g_rms_its[-1] < self.conv_tol:
                macro_converged = True
                break

            # generate trial vectors
            red_space_basis, ntrial_start = _generate_trial_vectors(
                self.trial_vectors, g, g_norm, hdiag, hess_x
            )

            # calculate linear transformation of basis vectors
            h_basis = [hess_x(vec) for vec in red_space_basis]

            # get array of basis vectors
            basis_arr = np.vstack(red_space_basis)

            # construct augmented Hessian in reduced space
            aug_hess = np.zeros(2 * (len(red_space_basis) + 1,), dtype=np.float64)
            for j, vec in enumerate(h_basis, start=1):
                aug_hess[1:, j] = basis_arr @ vec

            # decrease trust radius until micro iterations converge
            last_solution_normalized = np.zeros_like(g)
            micro_converged = False
            while not micro_converged:
                # loop over micro iterations
                for imicro in range(N_MICRO):
                    # perform bisection
                    solution, red_solution, mu = _bisection(
                        aug_hess, g, trust_radius, basis_arr
                    )

                    # calculate residual
                    residual = g + red_solution.T @ np.vstack(h_basis) - mu * solution

                    # get residual norm
                    residual_norm = np.linalg.norm(residual)

                    # determine reduction factor depending on whether local region is
                    # reached
                    red_factor = GLOBAL_RED_FACTOR if mu != 0.0 else LOCAL_RED_FACTOR

                    # get norm of solution vector
                    solution_normalized = solution / np.linalg.norm(solution)

                    # reset initial residual norm if state changes
                    if (last_solution_normalized.T @ solution_normalized) ** 2 < 0.5:
                        initial_imicro = imicro
                        initial_residual_norm = residual_norm

                    # check if micro iterations have converged
                    if (
                        residual_norm < red_factor * initial_residual_norm
                        or initial_residual_norm < 1e-10
                    ):
                        micro_converged = True
                        break
                    elif (
                        imicro - initial_imicro >= 5
                        and residual_norm > 0.9 * initial_residual_norm
                    ):
                        break

                    # save level shift
                    last_solution_normalized = solution_normalized.copy()

                    # precondition residual
                    precond = hdiag - mu
                    precond[abs(precond) < 1e-10] = 1e-10
                    residual = np.divide(residual, precond)

                    # orthogonalize to previous trial vectors
                    red_space_basis.append(_gram_schmidt(residual, red_space_basis))

                    # add linear transformation of new basis vector
                    h_basis.append(hess_x(red_space_basis[-1]))

                    # get array of basis vectors
                    basis_arr = np.vstack(red_space_basis)

                    # construct new augmented Hessian
                    aug_hess = _extend_symm_mat(
                        aug_hess, np.insert(basis_arr @ h_basis[-1], 0, 0)
                    )

                # decrease trust radius if micro iterations are unable to converge
                if not micro_converged:
                    aug_hess = aug_hess[: ntrial_start + 1, : ntrial_start + 1]
                    red_space_basis = red_space_basis[:ntrial_start]
                    h_basis = h_basis[:ntrial_start]
                    basis_arr = basis_arr[:ntrial_start]
                    trust_radius /= 2
                    self.tot_hx += len(red_space_basis)
                    if trust_radius < 1e-10:
                        self.logger.error("Trust radius too small.")
                        break
                else:
                    trust_radius = min(2 * trust_radius, CONS_TRUST_RADIUS)

            # check if macroiterations have converged
            if macro_converged or not micro_converged:
                break

            # increment total number of Hessian linear transformations
            self.tot_hx += len(red_space_basis)

            # get new orbital transformation
            kappa = solution

            # prepare functions for line search
            if direction == "min":
                func_eval = lambda n_kappa: func(
                    u @ sc.linalg.expm(unpack(n_kappa * kappa))
                )

                def func_grad_eval(n_kappa):
                    f, g = func_grad(u @ sc.linalg.expm(unpack(n_kappa * kappa)))
                    return f, g

            elif direction == "max":
                func_eval = lambda n_kappa: -func(
                    u @ sc.linalg.expm(unpack(n_kappa * kappa))
                )

                def func_grad_eval(n_kappa):
                    f, g = func_grad(u @ sc.linalg.expm(unpack(n_kappa * kappa)))
                    return -f, g

            # get maximum frequency, period and upper bound for minimum
            omega_max = np.max(np.abs(np.linalg.eig(unpack(kappa))[0]))
            period = 2 * np.pi / (line_search_orb_order * omega_max)
            upper_bound = period / 2

            # conduct a line search along the direction given by kappa
            if line_search == "golden":
                # perform golden section line search
                n_kappa, f = _golden_search((0.0, upper_bound), func_eval)
            elif line_search == "brent":
                # perform Brent's line search
                n_kappa, f = _brent_search((0.0, upper_bound), func_eval)
            elif line_search == "cubic":
                # perform cubic interpolation line search
                n_kappa, f, g = _cubic_search(
                    (0.0, upper_bound),
                    func_grad_eval,
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
            u @= sc.linalg.expm(unpack(n_kappa * kappa))

        if not macro_converged:
            raise RuntimeError("Orbital optimization has not converged!")

        if np.min(hdiag) < -self.conv_tol:
            raise RuntimeError("Orbital optimization has converged to saddle point!")

        self.logger.info(101 * "-")
        self.logger.info(
            f"Total number of Hessian linear transformations: {self.tot_hx}"
        )

        return u

    def stability(
        self,
        func_grad_hdiag_hess_x: Callable[
            [np.ndarray],
            Tuple[float, np.ndarray, np.ndarray, Callable[[np.ndarray], np.ndarray]],
        ],
        unpack: Callable[[np.ndarray], np.ndarray],
        u: np.ndarray,
    ):
        """
        this function performs a stability check by diagonalizing the provided Hessian
        """
        # evaluate gradient and Hessian at current point
        _, g, hdiag, hess_x = func_grad_hdiag_hess_x(u)

        # get quantities for full Hessian
        g *= 2
        hdiag *= 2

        def full_hess_x(x):
            return hess_x(x).real * 2

        g_norm = np.linalg.norm(g)

        # generate trial vectors
        red_space_basis = _generate_trial_vectors(
            self.trial_vectors, g, g_norm, hdiag, hess_x
        )[0]

        # calculate linear transformation of basis vectors
        h_basis = [full_hess_x(vec) for vec in red_space_basis]

        # get array of basis vectors
        basis_arr = np.vstack(red_space_basis)

        # construct reduced space Hessian
        red_space_hess = np.empty(2 * (len(red_space_basis),), dtype=np.float64)
        for j, vec in enumerate(h_basis):
            red_space_hess[:, j] = basis_arr @ vec

        # loop over micro iterations
        for _ in range(N_STAB_DAVIDSON):
            # solve reduced space problem
            eigenvalues, eigenvecs = np.linalg.eig(red_space_hess)

            # get index of lowest eigenvalue
            idx = np.argmin(eigenvalues)

            # get full space solution
            solution = eigenvecs[:, idx].T @ red_space_basis

            # calculate residual
            residual = (
                eigenvecs[:, idx].T @ np.vstack(h_basis) - eigenvalues[idx] * solution
            )

            # check convergence
            if np.linalg.norm(residual) < 1e-4:
                break

            # precondition residual
            precond = hdiag - eigenvalues[idx]
            precond[abs(precond) < 1e-10] = 1e-10
            residual = np.divide(residual, precond)

            # orthogonalize to previous trial vectors
            red_space_basis.append(_gram_schmidt(residual, red_space_basis))

            # add linear transformation of new basis vector
            h_basis.append(full_hess_x(red_space_basis[-1]))

            # get array of basis vectors
            basis_arr = np.vstack(red_space_basis)

            # construct new reduced space Hessian
            red_space_hess = _extend_symm_mat(red_space_hess, basis_arr @ h_basis[-1])

        # determine if saddle point
        stable = not eigenvalues[idx] < -1e-3
        if not stable:
            self.logger.warning(
                f"Solution not stable. Lowest eigenvalue: {eigenvalues[idx]}"
            )

        return stable, u @ sc.linalg.expm(unpack(solution))


def _generate_trial_vectors(
    trial_vectors: Union[str, int],
    g: np.ndarray,
    g_norm: any_float,
    hdiag: np.ndarray,
    hess_x: Callable[[np.ndarray], np.ndarray],
) -> Tuple[List[np.ndarray], int]:
    """
    this function generates trial vectors
    """
    red_space_basis = [g / g_norm]
    trial = -g - hess_x(g) / g_norm
    red_space_basis.append(_gram_schmidt(trial, red_space_basis))
    ntrial_start = 2
    if isinstance(trial_vectors, str):
        if trial_vectors == "orig":
            trial = np.zeros_like(g, dtype=np.float64)
            min_idx = np.argmin(hdiag)
            if hdiag[min_idx] < 0.0:
                trial[min_idx] = 1.0
                red_space_basis.append(_gram_schmidt(trial, red_space_basis))
                ntrial_start += 1
        elif trial_vectors == "lin_comb":
            trial = 1 / (hdiag - min(hdiag) + 0.1)
            trial /= np.linalg.norm(trial)
            red_space_basis.append(_gram_schmidt(trial, red_space_basis))
            ntrial_start += 1
    else:
        min_idx = np.argmin(hdiag)
        for min_idx in np.argsort(hdiag)[: trial_vectors - 2]:
            trial = np.zeros_like(g, dtype=np.float64)
            trial[min_idx] = 1.0
            red_space_basis.append(_gram_schmidt(trial, red_space_basis))
            ntrial_start += 1

    return red_space_basis, ntrial_start


def _extend_symm_mat(old_mat: np.ndarray, vector: np.ndarray) -> np.ndarray:
    """
    this function extends a symmetric matrix by a vector
    """
    new_mat = np.empty(2 * (old_mat.shape[0] + 1,), dtype=np.float64)
    new_mat[:-1, :-1] = old_mat
    new_mat[:, -1] = new_mat[-1] = vector
    return new_mat


def _gram_schmidt(vector: np.ndarray, space: List[np.ndarray]) -> np.ndarray:
    """
    this function orthonormalizes a vector with respect to a vector space
    """
    orth_vector = vector
    for space_vec in space:
        orth_vector -= (
            np.dot(vector, space_vec) / np.dot(space_vec, space_vec) * space_vec
        )

    return orth_vector / np.linalg.norm(orth_vector)


def _bisection(
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


def _golden_search(
    brack: Tuple[float, any_float], func_eval: Callable[[float], float]
) -> Tuple[float, float]:
    """
    this function performs a golden-section line search
    """
    # perform golden-section line search
    n_kappa, f, _ = sc.optimize.golden(
        func_eval, brack=brack, tol=1e-2, full_output=True
    )

    return n_kappa, f


def _brent_search(
    brack: Tuple[float, any_float], func_eval: Callable[[float], float]
) -> Tuple[float, float]:
    """
    this function performs a Brent's line search
    """
    # perform golden-section line search
    n_kappa, f, _, _ = sc.optimize.brent(
        func_eval, brack=brack, tol=1e-2, full_output=True
    )

    return n_kappa, f


def _cubic_search(
    brack: Tuple[any_float, any_float],
    func_grad_eval: Callable[[any_float], Tuple[float, np.ndarray]],
    kappa: np.ndarray,
    initial_f: Optional[float] = None,
    initial_g: Optional[np.ndarray] = None,
) -> Tuple[float, float, np.ndarray]:
    """
    this function performs a cubic interpolation line search
    """
    # get function value and gradient at current point
    if initial_f is None or initial_g is None:
        f1, g1 = func_grad_eval(brack[0])
    else:
        f1 = initial_f
        g1 = initial_g

    # calculate cost function and gradient at upper bound
    f2, g2 = func_grad_eval(brack[1])

    # get directional derivative in search direction, factor 2 comes from multiplying
    # the unpacked matrices
    g1 = 2 * g1.T @ kappa
    g2 = 2 * g2.T @ kappa

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
        radicand = d1**2 - g_low * g_high
        if radicand < 0.0:
            if radicand > -1e-16:
                radicand = np.abs(radicand)
            else:
                raise RuntimeError("Error in cubic extrapolation!")
        d2 = np.sign(n_kappa_high - n_kappa_low) * np.sqrt(radicand)
        n_kappa = n_kappa_high - (n_kappa_high - n_kappa_low) * (g_high + d2 - d1) / (
            g_high - g_low + 2 * d2
        )

        # calculate cost function and gradient at new point
        f, g = func_grad_eval(n_kappa)

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


def _create_logger(verbosity: int):
    """
    Create a logger with a given verbosity level.
    """
    # basic config
    logging.basicConfig(format="%(message)s")

    # map verbosity level to logging levels
    log_levels = {
        3: logging.DEBUG,
        2: logging.INFO,
        1: logging.WARNING,
        0: logging.ERROR,
        0: logging.CRITICAL,
    }

    # Create a logger instance
    logger = logging.getLogger("pyorbopt_logger")
    logger.setLevel(log_levels[verbosity])

    return logger
