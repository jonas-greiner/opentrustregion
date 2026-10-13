# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import unittest
from ctypes import c_void_p, POINTER, cast

from pyopentrustregion.tests import (
    lib,
    n_param,
    NUMPY_AVAILABLE,
    add_tests,
    print_separator,
)
from pyopentrustregion.python_interface import c_real, c_int, SolverSettings

if NUMPY_AVAILABLE:
    import numpy as np
    from pyopentrustregion.extensions.common.python_interface import (
        EvaluateDMInterface,
    )


# define all tests in alphabetical order
fortran_tests = {
    "common_tests": [
        "approx_hess_x_callback",
        "channel_rows",
        "check_unwired_extension",
        "compute_sqrt_and_inv_sqrt",
        "conv_check_callback",
        "conv_check_stability_callback",
        "extension_level",
        "fill_extra_trial_vectors_orbital_basis",
        "get_extra_trial_vectors_callback",
        "get_extra_trial_vectors_extension",
        "hess_x_callback",
        "hess_x_extension",
        "init_orbital_settings",
        "level_shifted_divisors",
        "logger_callback",
        "matrix_exponential",
        "obj_func_callback",
        "obj_func_extension",
        "positive_definite_divisors",
        "precond_callback",
        "precond_extension",
        "precond_orbital_basis",
        "precond_pd_callback",
        "precond_pd_extension",
        "precond_pd_orbital_basis",
        "project_callback",
        "project_extension",
        "refresh_response_orbital_basis",
        "stability_hess_x_callback",
        "update_orbs_callback",
        "update_orbs_extension",
        "wire_solver_settings_extension",
    ],
    "common_c_interface_tests": [
        "approx_hess_x_c_wrapper",
        "conv_check_c_wrapper",
        "conv_check_stability_c_wrapper",
        "get_extra_trial_vectors_c_wrapper",
        "hess_x_c_wrapper",
        "holder_level",
        "logger_c_wrapper",
        "obj_func_c_wrapper",
        "precond_c_wrapper",
        "precond_pd_c_wrapper",
        "project_c_wrapper",
        "setup_extension_holder",
        "stability_hess_x_c_wrapper",
        "update_orbs_c_wrapper",
        "wire_solver_settings_c",
    ],
}

# the shared extension routines are only built together with an extension which uses
# them, so raise an AttributeError on import otherwise
getattr(lib, "test_" + fortran_tests["common_tests"][0])

# number of AOs
n_ao = c_int.in_dll(lib, "test_n_ao").value


def setup_wired_solver_settings():
    """
    this function returns solver settings for the Python test of a factory that wires
    them, whose projection and logging function and the convergence check of their
    stability check are built from Python functions, which the factory has to be handed
    without, while their preconditioner is set directly, which it has to be handed,
    together with these Python functions
    """
    solver_settings = SolverSettings()
    solver_settings.project = lambda vector: None
    solver_settings.logger = lambda message: None
    solver_settings.stability_settings.conv_check = lambda residual, eigval: False
    solver_settings.set_optional_callbacks(n_param, {})
    solver_settings.settings_c.precond = cast(lib.mock_level_precond, c_void_p).value
    python_functions = (
        solver_settings.project,
        solver_settings.logger,
        solver_settings.stability_settings.conv_check,
    )

    return solver_settings, python_functions


def check_wired_solver_settings(solver_settings, python_functions, test_name):
    """
    this function checks that solver settings a mock factory wired still hold their
    Python functions, from which the solver builds the callback functions anew, and
    that the preconditioners the mock factory wired into both levels are called with
    the contexts of their levels
    """
    test_passed = True

    if (
        solver_settings.project,
        solver_settings.logger,
        solver_settings.stability_settings.conv_check,
    ) != python_functions:
        print(f" test_{test_name} failed: Python functions of the settings not kept.")
        test_passed = False

    residual = np.full(n_param, 1.0, dtype=np.float64)
    precond_residual = np.empty(n_param, dtype=np.float64)
    for level, settings in (
        ("solver", solver_settings),
        ("stability check", solver_settings.stability_settings),
    ):
        precond_residual[:] = 0.0
        try:
            settings.precond(residual, 2.0, precond_residual)
        except RuntimeError:
            print(
                f" test_{test_name} failed: Wired preconditioner of the {level} "
                "raises error."
            )
            test_passed = False
            continue
        if not np.allclose(precond_residual, 2.0):
            print(
                f" test_{test_name} failed: Wired preconditioner of the {level} "
                "returns wrong result."
            )
            test_passed = False

    return test_passed


@add_tests
class CommonTests(unittest.TestCase):
    """
    this class contains unit tests for the routines shared between extensions
    """

    tests = fortran_tests["common_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for shared extension routines...")
        return super().setUpClass()


@add_tests
class CommonCInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the C interface shared between extensions
    """

    tests = fortran_tests["common_c_interface_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for shared extension C interface...")
        return super().setUpClass()


@unittest.skipUnless(NUMPY_AVAILABLE, "NumPy not available.")
class CommonPyInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the Python interface shared between extensions
    """

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for shared extension Python interface...")
        return super().setUpClass()

    def test_evaluate_dm_py_interface(self):
        """
        this function tests the Python interface to the density matrix evaluating
        function
        """
        # initialize test flag
        test_passed = True

        # check if a density matrix evaluating function which returns no response
        # function although one is requested produces an error
        def evaluate_dm_without_response(dm_ao, fock, get_response):
            return np.sum(dm_ao), None

        exception = {}
        evaluate_dm = EvaluateDMInterface(
            evaluate_dm_without_response, n_ao, 1, True, exception
        )
        dm_ao = np.full(2 * (n_ao,), 1.0, dtype=np.float64)
        energy = (c_real * 1)()
        get_response_funptr = (c_void_p * 1)()
        error = evaluate_dm(
            dm_ao.ctypes.data_as(POINTER(c_real)),
            energy,
            None,
            get_response_funptr,
            None,
        )
        if error == 0 or not isinstance(exception.get("exc"), RuntimeError):
            print(
                " test_evaluate_dm_py_interface failed: Missing response function "
                "requested from density matrix evaluating function does not produce "
                "an error."
            )
            test_passed = False

        self.assertTrue(test_passed, "test_evaluate_dm_py_interface failed")
        print(" test_evaluate_dm_py_interface PASSED")
