from __future__ import annotations


import logging
import numpy as np
import scipy as sc
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import List, Tuple, Callable, Union

    any_float = Union[float, np.floating]


class OrbOpt:
    def __init__(self, verbose: int = 1, random_seed: int = 42):
        # initialize total number of Hessian linear transformations
        self.tot_hx: int = 0
        # initialize list for function evaluations per macro iteration
        self.func_its: List[float] = []
        # initialize list for gradient RMS per macro iteration
        self.g_rms_its: List[any_float] = []
        # initialize verbosity setting
        self.logger = _create_logger(verbose)
        # initialize random number generator
        self.rng = np.random.default_rng(seed=random_seed)

    def solver(
        self,
        func: Callable[[np.ndarray], float],
        update_orbs: Callable[
            [np.ndarray],
            Tuple[float, np.ndarray, np.ndarray, Callable[[np.ndarray], np.ndarray]],
        ],
        n_param: int,
        stability: bool = True,
        line_search: bool = False,
        conv_tol: float = 1e-5,
        n_random_trial_vectors: int = 1,
        start_trust_radius: float = 0.4,
        n_macro: int = 150,
        n_micro: int = 100,
        global_red_factor: float = 1e-3,
        local_red_factor: float = 1e-4,
    ):
        # initialize orbital rotation matrix
        kappa = np.zeros(n_param, dtype=np.float64)

        # check that number of random trial vectors is below number of parameters
        if n_random_trial_vectors > n_param // 2:
            n_random_trial_vectors = n_param // 2
            self.logger.warning(
                "Number of random trial vectors should be smaller than half the number "
                f"of parameters. Setting to {n_random_trial_vectors}."
            )

        # initialize starting trust radius
        trust_radius = start_trust_radius

        # inititalize boolean for convergence of macro iterations
        macro_converged = False

        # initialize boolean for stability
        stable = True

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

        for imacro in range(n_macro):
            # calculate cost function, gradient and Hessian diagonal
            f, g, hdiag, hess_x = update_orbs(kappa)

            # sanity check for array size
            if imacro == 0 and g.size != n_param:
                raise ValueError(
                    "Size of gradient array returned by function update_orbs does not "
                    "equal number of parameters"
                )

            # add cost function
            self.func_its.append(float(f))

            # calculate gradient norm
            g_norm = np.linalg.norm(g)

            # RMS gradient
            self.g_rms_its.append(float(g_norm / np.sqrt(n_param)))

            # log results
            if imacro == 0:
                mu_str = f"{'-':^9}"
                imicro_str = f"{'-':>3}"
                trust_radius_str = f"{'-':^8}"
                stepsize_str = f"{'-':^8}"
            elif not stable:
                mu_str = f"{'-':^9}"
                imicro_str = f"{'-':>3}"
                trust_radius_str = f"{'-':^8}"
                stepsize_str = f"{np.linalg.norm(kappa):>8.2e}"
                stable = True
            else:
                mu_str = f"{mu:>9.2e}"
                imicro_str = f"{imicro:>3d}"
                trust_radius_str = f"{trust_radius:>8.2e}"
                stepsize_str = f"{np.linalg.norm(kappa):>8.2e}"
            self.logger.info(
                f"     {imacro:>3d}   |    {f:> 8.6e}   "
                f"|   {self.g_rms_its[-1]:>8.2e}   |  {mu_str}  |      {imicro_str}   "
                f"|    {trust_radius_str}  |  {stepsize_str} "
            )

            # check for convergence and stability
            if self.g_rms_its[-1] < conv_tol:
                if stability:
                    stable, kappa = self.stability(g, hdiag, hess_x)
                    if not stable:
                        self.logger.warning(
                            "Reached saddle point. This is likely due to symmetry and "
                            "can be avoided by increasing the number of random trial "
                            "vectors. The algorithm will continue by moving along "
                            "eigenvector direction corresponding to negative "
                            "eigenvalue."
                        )
                        continue
                    else:
                        macro_converged = True
                        break
                else:
                    macro_converged = True
                    break

            # generate trial vectors
            red_space_basis, ntrial_start = self._generate_trial_vectors(
                n_random_trial_vectors, g, g_norm, hdiag
            )

            # calculate linear transformation of basis vectors
            h_basis = [hess_x(vec) for vec in red_space_basis]

            # get array of basis vectors
            basis_arr = np.vstack(red_space_basis)

            # construct augmented Hessian in reduced space
            aug_hess = np.zeros(2 * (len(red_space_basis) + 1,), dtype=np.float64)
            for j, vec in enumerate(h_basis, start=1):
                aug_hess[1:, j] = basis_arr @ vec

            # decrease trust radius until micro iterations converge and step is accepted
            last_solution_normalized = np.zeros_like(g)
            accept_step = False
            while not accept_step:
                micro_converged = False
                # loop over micro iterations
                for imicro in range(n_micro):
                    # do a Newton step if the model is positive definite and the step is within the trust region
                    newton_step = False
                    eigvals = np.linalg.eigvalsh(aug_hess[1:, 1:])
                    if min(eigvals) > -1e-5:
                        solution, red_solution = _newton_step(aug_hess, g, basis_arr)
                        mu = 0.0
                        if np.linalg.norm(solution) < trust_radius:
                            newton_step = True

                    # otherwise perform bisection to find the level shift
                    if not newton_step:
                        solution, red_solution, mu = _bisection(
                            aug_hess, g, trust_radius, basis_arr
                        )

                    # calculate residual
                    residual = g + red_solution.T @ np.vstack(h_basis) - mu * solution

                    # get residual norm
                    residual_norm = np.linalg.norm(residual)

                    # determine reduction factor depending on whether local region is
                    # reached
                    red_factor = global_red_factor if mu != 0.0 else local_red_factor

                    # get norm of solution vector
                    solution_normalized = solution / np.linalg.norm(solution)

                    # reset initial residual norm if state changes
                    if (last_solution_normalized.T @ solution_normalized) ** 2 < 0.5:
                        initial_imicro = imicro
                        initial_residual_norm = residual_norm

                    # check if micro iterations have converged
                    if residual_norm < max(red_factor * g_norm, 1e-12):
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

                # evaluate function at predicted point
                new_f = func(solution)

                # calculate ratio of evaluated function and predicted function
                ratio = (new_f - f) / (
                    solution.T @ g
                    + 0.5
                    * (red_solution.T @ basis_arr).T
                    @ (red_solution @ np.vstack(h_basis))
                )

                # decrease trust radius if micro iterations are unable to converge or
                # if function value has not decreased or if individual orbitals change too much
                if (
                    not micro_converged
                    or ratio < 0.0
                    or np.any(np.abs(solution) > np.pi / 4)
                ):
                    aug_hess = aug_hess[: ntrial_start + 1, : ntrial_start + 1]
                    red_space_basis = red_space_basis[:ntrial_start]
                    h_basis = h_basis[:ntrial_start]
                    basis_arr = basis_arr[:ntrial_start]
                    trust_radius *= 0.7
                    self.tot_hx += len(red_space_basis)
                    accept_step = False
                    if trust_radius < 1e-10:
                        self.logger.error("Trust radius too small.")
                        raise RuntimeError
                # check if step is too long
                elif ratio < 0.25:
                    trust_radius *= 0.7
                    accept_step = True
                # check if quadratic approximation is valid
                elif ratio < 0.75:
                    accept_step = True
                # check if step is potentially too short
                else:
                    trust_radius *= 1.2
                    accept_step = True

            # check if step is not accepted
            if not accept_step:
                continue

            # increment total number of Hessian linear transformations
            self.tot_hx += len(red_space_basis)

            # turn off in local region
            if line_search:
                # prepare functions for line search
                func_eval = lambda n_kappa: func(n_kappa * solution)

                # conduct a bracket search along the direction given by kappa
                _, n_kappa, _, _, f, _, _ = sc.optimize.bracket(
                    func_eval, xa=0.0, xb=1.0
                )

            else:
                n_kappa = 1.0

            # set orbital rotation
            kappa = n_kappa * solution

        if not macro_converged:
            raise RuntimeError("Orbital optimization has not converged!")

        self.logger.info(101 * "-")
        self.logger.info(
            f"Total number of Hessian linear transformations: {self.tot_hx}"
        )

        return

    def stability(
        self,
        g: np.ndarray,
        hdiag: np.ndarray,
        hess_x: Callable[[np.ndarray], np.ndarray],
        perturb: float = 1e-2,
        conv_tol: float = 1e-4,
        n_random_trial_vectors: int = 20,
        n_iter: int = 100,
    ):
        """
        this function performs a stability check by diagonalizing the provided Hessian
        """
        # get quantities for full Hessian
        g *= 2
        hdiag *= 2

        def full_hess_x(x):
            return hess_x(x).real * 2

        g_norm = np.linalg.norm(g)

        # generate trial vectors
        red_space_basis = self._generate_trial_vectors(
            n_random_trial_vectors, g, g_norm, hdiag
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
        for _ in range(n_iter):
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
            if np.linalg.norm(residual) < conv_tol:
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

        if stable:
            return stable, np.zeros_like(solution)
        else:
            return stable, perturb * solution

    def _generate_trial_vectors(
        self,
        n_random_trial_vectors: int,
        g: np.ndarray,
        g_norm: any_float,
        hdiag: np.ndarray,
    ) -> Tuple[List[np.ndarray], int]:
        """
        this function generates trial vectors
        """
        red_space_basis = [g / g_norm]
        min_idx = np.argmin(hdiag)
        if hdiag[min_idx] < 0.0:
            trial = np.zeros_like(g, dtype=np.float64)
            trial[min_idx] = 1.0
            red_space_basis.append(_gram_schmidt(trial, red_space_basis))
        for _ in range(n_random_trial_vectors):
            trial = 2 * self.rng.random(size=g.size) - 1
            while np.linalg.norm(trial) < 1e-3:
                trial = 2 * self.rng.random(size=g.size) - 1
            red_space_basis.append(_gram_schmidt(trial, red_space_basis))

        return red_space_basis, len(red_space_basis)


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
        eigenvalues, eigenvecs = np.linalg.eigh(aug_hess)

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
    alpha = sc.optimize.bisect(
        get_ah_lowest_eigenvec_norm_trust_radius_diff,
        1e-4,
        1e6,
        args=(aug_hess, grad, red_space_basis),
    )

    # return level-shifted solution
    return get_ah_lowest_eigenvec(aug_hess, grad, alpha, red_space_basis)


def _newton_step(
    aug_hess: np.ndarray, grad: np.ndarray, red_space_basis: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    this function performs bisection to find the parameter alpha that matches the
    desired trust radius
    """
    # solve Newton equations in reduced space without level shift
    red_grad = np.zeros(len(red_space_basis))
    red_grad[0] = np.linalg.norm(grad)
    red_solution = np.linalg.solve(aug_hess[1:, 1:], -red_grad)
    solution = red_solution.T @ red_space_basis
    return solution, red_solution


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
