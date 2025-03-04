import sys
from ctypes import CDLL, c_bool
from unittest.mock import patch
import numpy as np
from pytrustorbopt import solver_py_interface, stability_check_py_interface

# load calling executable
ext = "dylib" if sys.platform == "darwin" else "so"
libtests = CDLL(sys.argv[1])

# set signature for functions which extracts the test results
libtests.test_solver_result.restypes = c_bool
libtests.test_stability_check_result.restypes = c_bool


# replace original library with mock library
@patch("pytrustorbopt.python_interface.libtrustorbopt.solver", libtests.mock_solver)
def test_solver_py_interface():
    """
    this function tests the solver python interface
    """

    def mock_obj_func(kappa):
        """
        this function is a mock function for the objective function
        """
        return np.sum(kappa)

    def mock_update_orbs(kappa):
        """
        this function is a mock function for the orbital update function
        """
        func = np.sum(kappa)
        grad = 2 * kappa
        h_diag = 3 * kappa

        def hess_x(x):
            return 4 * x

        return func, grad, h_diag, hess_x

    # call solver python interface without optional arguments
    solver_py_interface(mock_obj_func, mock_update_orbs, 3)

    # call solver python interface with optional arguments
    solver_py_interface(
        mock_obj_func,
        mock_update_orbs,
        3,
        stability=False,
        line_search=True,
        conv_tol=1e-3,
        n_random_trial_vectors=5,
        start_trust_radius=0.2,
        n_macro=300,
        n_micro=200,
        global_red_factor=1e-2,
        local_red_factor=1e-3,
        verbose=3,
        seed=33,
    )


# replace original library with mock library
@patch(
    "pytrustorbopt.python_interface.libtrustorbopt.stability_check",
    libtests.mock_stability_check,
)
def test_stability_check_py_interface():
    """
    this function tests the stability check python interface
    """
    grad = np.full(3, 2.0, dtype=np.float64)
    h_diag = np.full(3, 3.0, dtype=np.float64)

    def hess_x(x):
        return 4 * x

    # call stability check python interface without optional arguments
    stable, kappa = stability_check_py_interface(grad, h_diag, hess_x, 3)

    # check if returned variables are correct
    test_stability_check_interface_output = not stable and np.allclose(
        kappa, np.full(3, 1.0, dtype=np.float64)
    )

    # call stability check python interface with optional arguments
    stable, kappa = stability_check_py_interface(
        grad,
        h_diag,
        hess_x,
        3,
        conv_tol=1e-3,
        n_random_trial_vectors=3,
        n_iter=50,
        verbose=3,
    )

    # check if returned variables are correct
    test_stability_check_interface_output = (
        test_stability_check_interface_output
        and not stable
        and np.allclose(kappa, np.full(3, 1.0, dtype=np.float64))
    )
    if stable:
        print(
            "test_stability_check_py_interface failed: Returned stability boolean "
            "wrong."
        )
    if not np.allclose(kappa, np.full(3, 1.0, dtype=np.float64)):
        print("test_stability_check_py_interface failed: Returned direction wrong.")

    return test_stability_check_interface_output


if __name__ == "__main__":

    # run tests and extract results
    test_solver_py_interface()
    test_solver_interface = libtests.test_solver_result()
    test_stability_check_interface_output = test_stability_check_py_interface()
    test_stability_check_interface_input = libtests.test_stability_check_result()

    # print results
    n_test_results = 0
    if test_solver_interface:
        print(" solver_py_interface PASSED")
    else:
        n_test_results += 1
        print(" solver_py_interface FAILED")
    if test_stability_check_interface_input and test_stability_check_interface_output:
        print(" stability_check_py_interface PASSED")
    else:
        n_test_results += 1
        print(" stability_check_py_interface FAILED")

    # exit and pass number failed tests
    sys.exit(n_test_results)
