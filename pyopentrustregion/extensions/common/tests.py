# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import functools
import gc
import unittest
from ctypes import POINTER, addressof, c_void_p
from unittest.mock import patch

import numpy as np

from pyopentrustregion.tests import (
    lib,
    n_param,
    add_tests,
    print_separator,
)
from pyopentrustregion.python_interface import (
    c_real,
    c_int,
    hess_x_interface_type,
    SolverSettings,
    SolverState,
)

from pyopentrustregion.extensions.common.python_interface import (
    EvaluateDMInterface,
    ExtensionOwner,
    ObjFuncPyInterface,
    UpdateOrbsPyInterface,
    check_wrapped_functions,
)

# define all tests in alphabetical order
fortran_tests = {
    "common_tests": [
        "approx_hess_x_callback",
        "channel_rows",
        "check_evaluate_dm_orbital_basis",
        "compute_sqrt_and_inv_sqrt",
        "conv_check_callback",
        "conv_check_stability_callback",
        "extension_level",
        "fill_extra_trial_vectors_orbital_basis",
        "get_extra_trial_vectors_callback",
        "get_extra_trial_vectors_extension",
        "hess_x_callback",
        "hess_x_extension",
        "init_extension_solver_settings",
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
        "wire_extension_solver_settings",
        "wire_solver_settings_extension",
    ],
    "common_c_interface_tests": [
        "approx_hess_x_c_wrapper",
        "assign_solver_state_f_c",
        "assign_stability_state_f_c",
        "bundle_settings_c",
        "conv_check_c_wrapper",
        "conv_check_stability_c_wrapper",
        "evaluate_dm_cs_f_wrapper",
        "evaluate_dm_os_f_wrapper",
        "finish_extension_holder",
        "get_extra_trial_vectors_c_wrapper",
        "get_response_cs_f_wrapper",
        "get_response_os_f_wrapper",
        "hess_x_c_wrapper",
        "holder_level",
        "logger_c_wrapper",
        "obj_func_c_wrapper",
        "precond_c_wrapper",
        "precond_pd_c_wrapper",
        "project_c_wrapper",
        "refresh_state",
        "release_extension",
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


def release_mock_extension(extension):
    """
    this function frees an extension a mock factory returned through the mock release
    of the library, so that the owners the Python factories set up in their tests,
    whose finalizers keep the release they were set up with, never call the release of
    the library, also when they are finalized after the test
    """
    lib.mock_release_extension.restype = None
    lib.mock_release_extension.argtypes = [c_void_p]
    lib.mock_release_extension(extension)


def setup_wired_solver_settings():
    """
    this function returns solver settings for the Python test of a factory that wires
    them, with a projection and logging function and a convergence check of their
    stability check, whose C callback functions the factory has to be handed, together
    with these Python functions
    """
    solver_settings = SolverSettings()
    solver_settings.project = lambda vector: None
    solver_settings.logger = lambda message: None
    solver_settings.stability_settings.conv_check = lambda residual, eigval: False
    python_functions = (
        solver_settings.project,
        solver_settings.logger,
        solver_settings.stability_settings.conv_check,
    )

    return solver_settings, python_functions


def check_wired_solver_settings(
    solver_settings, python_functions, test_name, returned, extension_settings
):
    """
    this function checks that solver settings a mock factory wired still hold their
    Python functions, that the preconditioners the mock factory wired into both levels
    are called with the contexts of their levels, that a Python function assigned
    after the wiring is called with the context the mock factory wrote into the
    settings, and that the functions the factory returned keep the solver settings and
    the settings of the extension alive, whose callback functions its object calls
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

    # assign a Hessian linear transformation for the stability check after the wiring
    # and call it through the C settings with the context the mock factory wrote and
    # the number of parameters and exception collector the solver would hand it
    x = np.ones(n_param, dtype=np.float64)
    hess_x = np.zeros(n_param, dtype=np.float64)

    def stability_hess_x(x, hess_x):
        hess_x[:] = 4 * x

    solver_settings.stability_hess_x = stability_hess_x
    solver_settings._n_param = n_param
    solver_settings._exception = {}
    error = hess_x_interface_type(solver_settings.settings_c.stability_hess_x)(
        x.ctypes.data_as(POINTER(c_real)),
        hess_x.ctypes.data_as(POINTER(c_real)),
        solver_settings.settings_c.context,
    )
    if error != 0 or not np.allclose(hess_x, 4.0):
        print(
            f" test_{test_name} failed: Python function assigned after the wiring "
            "not called."
        )
        test_passed = False
    solver_settings.stability_hess_x = None

    # the returned functions keep the solver settings alive and the extension's
    # settings directly or through the owner of the extension, which the solver
    # settings keep alive
    owner = getattr(solver_settings, "_extension_owner", None)
    for function in returned:
        kept_settings = function.saved_objects.get("settings")
        if kept_settings is None and owner is not None:
            kept_settings = owner.saved_objects.get("settings")
        if (
            function.saved_objects.get("solver_settings") is not solver_settings
            or kept_settings is not extension_settings
        ):
            print(
                f" test_{test_name} failed: Returned function does not keep the "
                "settings alive."
            )
            test_passed = False

    return test_passed


def check_extension_owner(call_factory, test_name):
    """
    this function checks that a mock factory, which call_factory calls with solver
    settings and a solver state, received the record of the state and handed the
    extension it returned to an owner, which the solver settings keep alive and which
    keeps the state alive
    """
    test_passed = True
    received_solver_state = c_void_p.in_dll(lib, "received_solver_state")
    returned_extension = c_void_p.in_dll(lib, "returned_extension")

    # call the factory
    solver_settings = setup_wired_solver_settings()[0]
    solver_state = SolverState()
    call_factory(solver_settings, solver_state)

    # the factory received the record of the state, and the solver settings keep the
    # owner of the extension it returned, which keeps the state alive
    if received_solver_state.value != addressof(solver_state.state_c):
        print(f" test_{test_name} failed: Record of the solver state not passed.")
        test_passed = False
    owner = getattr(solver_settings, "_extension_owner", None)
    if (
        owner is None
        or owner.extension != returned_extension.value
        or owner.saved_objects.get("solver_state") is not solver_state
    ):
        print(f" test_{test_name} failed: Returned extension not owned.")
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


class CommonPyInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the Python interface shared between extensions
    """

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for shared extension Python interface...")
        return super().setUpClass()

    def test_check_wrapped_functions(self):
        """
        this function tests the function which checks that a factory wraps the orbital
        updating function and the objective function the factory that wired the solver
        settings last returned, or the host's own for settings no factory wired
        """
        # initialize test flag
        test_passed = True

        # host's functions and functions returned by factories for the levels at two
        # addresses
        def host_update_orbs(kappa, grad, h_diag):
            return 0.0, None

        def host_obj_func(kappa):
            return 0.0

        def marked_wrapper(function):
            @functools.wraps(function)
            def wrapper(*args):
                return function(*args)

            return wrapper

        level, other_level = 0x1000, 0x2000
        returned_update_orbs = {
            address: UpdateOrbsPyInterface(
                update_orbs_funptr=None, saved_objects={}, context=address
            )
            for address in (level, other_level)
        }
        returned_obj_func = {
            address: ObjFuncPyInterface(
                obj_func_funptr=None, saved_objects={}, context=address
            )
            for address in (level, other_level)
        }

        # cases of the context of the settings, None for settings no factory wired,
        # the functions passed and whether they have to be accepted
        cases = (
            (
                "host's functions for settings no factory wired",
                None,
                host_update_orbs,
                host_obj_func,
                True,
            ),
            (
                "returned functions for settings no factory wired",
                None,
                returned_update_orbs[level],
                returned_obj_func[level],
                False,
            ),
            (
                "a returned orbital updating function for settings no factory wired",
                None,
                returned_update_orbs[level],
                host_obj_func,
                False,
            ),
            (
                "a returned objective function for settings no factory wired",
                None,
                host_update_orbs,
                returned_obj_func[level],
                False,
            ),
            (
                "functions returned for the wired settings",
                level,
                returned_update_orbs[level],
                returned_obj_func[level],
                True,
            ),
            (
                "host's functions for wired settings",
                level,
                host_update_orbs,
                host_obj_func,
                False,
            ),
            (
                "the host's orbital updating function for wired settings",
                level,
                host_update_orbs,
                returned_obj_func[level],
                False,
            ),
            (
                "the host's objective function for wired settings",
                level,
                returned_update_orbs[level],
                host_obj_func,
                False,
            ),
            (
                "functions returned for settings wired before",
                level,
                returned_update_orbs[other_level],
                returned_obj_func[other_level],
                False,
            ),
            (
                "wrappers made with functools.wraps of the functions returned for the "
                "wired settings",
                level,
                marked_wrapper(returned_update_orbs[level]),
                marked_wrapper(returned_obj_func[level]),
                True,
            ),
            (
                "wrappers made with functools.wraps of functions returned for settings "
                "wired before",
                level,
                marked_wrapper(returned_update_orbs[other_level]),
                marked_wrapper(returned_obj_func[other_level]),
                False,
            ),
            (
                "unmarked wrappers of the functions returned for the wired settings",
                level,
                lambda *args: returned_update_orbs[level](*args),
                lambda *args: returned_obj_func[level](*args),
                False,
            ),
        )
        for case, context, update_orbs, obj_func, accepted in cases:
            solver_settings = SolverSettings()
            solver_settings.settings_c.context = context
            try:
                check_wrapped_functions(update_orbs, obj_func, solver_settings)
                raised = False
            except ValueError:
                raised = True
            if raised == accepted:
                print(
                    f" test_check_wrapped_functions failed: {case} "
                    f"{'refused' if accepted else 'accepted'}."
                )
                test_passed = False

        self.assertTrue(test_passed, "test_check_wrapped_functions failed")
        print(" test_check_wrapped_functions PASSED")

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

    @patch(
        "pyopentrustregion.python_interface.lib.release_extension",
        lib.mock_release_extension,
    )
    def test_extension_owner(self):
        """
        this function tests the owner of the extension object a factory returned for
        the solver settings it wired
        """
        # initialize test flag
        test_passed = True
        released_extension = c_void_p.in_dll(lib, "released_extension")
        n_released_extensions = c_int.in_dll(lib, "n_released_extensions")
        lib.mock_extension.restype = c_void_p
        lib.mock_extension.argtypes = [c_void_p]

        # let two extensions be owned one after the other for the same solver settings
        solver_settings = SolverSettings()
        kept = object()
        extensions = (lib.mock_extension(None), lib.mock_extension(None))
        first = ExtensionOwner(extensions[0], solver_settings, {"kept": kept})
        second = ExtensionOwner(extensions[1], solver_settings, {})

        # the settings keep the latest owner alive, which keeps the earlier one alive,
        # which keeps what it was given alive
        if (
            getattr(solver_settings, "_extension_owner", None) is not second
            or second.saved_objects.get("wired_before") is not first
            or first.saved_objects.get("wired_before") is not None
            or first.saved_objects.get("kept") is not kept
            or first.extension != extensions[0]
        ):
            print(" test_extension_owner failed: Owners not kept alive.")
            test_passed = False

        # dropping the earlier owner releases nothing, since the latest keeps it alive,
        # after releasing the owners earlier tests left in reference cycles
        gc.collect()
        n_released_extensions.value = 0
        del first
        gc.collect()
        if n_released_extensions.value != 0:
            print(" test_extension_owner failed: Owner wired over released too early.")
            test_passed = False

        # dropping the settings and the latest owner releases both, the latest first
        del second, solver_settings
        gc.collect()
        if (
            n_released_extensions.value != 2
            or released_extension.value != extensions[0]
        ):
            print(
                " test_extension_owner failed: Extensions not released in reverse order "
                "of wiring."
            )
            test_passed = False

        self.assertTrue(test_passed, "test_extension_owner failed")
        print(" test_extension_owner PASSED")
