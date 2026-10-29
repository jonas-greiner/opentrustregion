# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import gc
import os
import sys
import copy
import signal
import threading
import unittest
import warnings
import weakref
from ctypes import (
    c_bool,
    pointer,
    c_char_p,
    c_void_p,
    sizeof,
    memset,
    addressof,
    Array,
    byref,
    Structure,
    POINTER,
    cast,
    string_at,
)
from unittest.mock import patch
from functools import partial
from types import SimpleNamespace
from pathlib import Path
import numpy as np

# when run as a script, import the package this file belongs to rather than another
# one that is installed
if not __package__:
    sys.path.insert(0, str(Path(__file__).parent.parent.absolute()))
from pyopentrustregion.python_interface import (
    SolverSettings,
    StabilitySettings,
    SolverState,
    StabilityState,
    solver,
    stability_check,
    c_int,
    c_real,
    ext,
    _load_lib,
    interruption_state,
    defer_interruption,
    ObjFuncInterface,
    HessXInterface,
    UpdateOrbsInterface,
    c_buffer,
    write_back,
    raise_on_failure,
    raise_callback_error,
    call_library,
    check_struct_sizes,
    hess_x_interface_type,
    precond_interface_type,
    precond_pd_interface_type,
    project_interface_type,
    get_extra_trial_vectors_interface_type,
    conv_check_interface_type,
    logger_interface_type,
)

# load the testsuite library
lib = _load_lib([f"libotrtestsuite.{ext}"])


def print_separator(string: str):
    sep_length = 70
    print(sep_length * "-")
    print(string)
    print(sep_length * "-")


# define all tests in alphabetical order
fortran_tests = {
    "opentrustregion_unit_tests": [
        "accept_trust_region_step",
        "add_column",
        "add_error_origin",
        "add_trial_vector",
        "backward_substitution",
        "bisection_ah",
        "bisection_mu",
        "block_davidson",
        "bracket",
        "check_stationary_point",
        "extend_matrix",
        "find_point_on_boundary",
        "forward_substitution",
        "general_mat_diag",
        "generalized_lanczos_trust_region",
        "generate_random_trial_vectors",
        "generate_trial_vectors",
        "get_preconditioned_residual_dot",
        "get_stability_trial_space",
        "get_tridiagonal_lowest_eigenvalue",
        "get_tridiagonal_lowest_eigenvector",
        "gltr_assemble_solution",
        "gltr_lanczos_process",
        "gltr_precondition_residual",
        "gram_schmidt",
        "inherit_solver_callbacks",
        "init_defaults",
        "init_rng",
        "invoke_get_extra_trial_vectors",
        "invoke_hess_x",
        "invoke_obj_func",
        "invoke_precond",
        "invoke_precond_pd",
        "invoke_project",
        "invoke_update_orbs",
        "jacobi_davidson_correction",
        "level_shifted_davidson",
        "level_shifted_diag_precond",
        "lowest_h_diag_unit_vectors",
        "mat_diag",
        "mat_lowest_eigpairs",
        "mat_min_eig",
        "micro_conv_tol",
        "minres",
        "newton_step",
        "noise_level",
        "option_list",
        "orthogonal_projection",
        "orthogonalize_trial_vectors",
        "perturb_vector",
        "print_message",
        "print_results",
        "quadratic_roots",
        "rel_floor_diag_precond",
        "saddle_point_step",
        "solve_factorized_tridiagonal",
        "solve_lanczos_subproblem",
        "solve_tridiagonal",
        "solve_tridiagonal_subproblem",
        "solver",
        "solver_sanity_check",
        "split_string_by_space",
        "stability_check",
        "stability_sanity_check",
        "string_in",
        "string_to_lowercase",
        "symm_mat_diag",
        "truncated_conjugate_gradient",
        "warn_uninitialized_settings",
    ],
    "c_interface_unit_tests": [
        "approx_hess_x_f_wrapper",
        "assign_solver_c_f",
        "assign_solver_f_c",
        "assign_solver_state_c_f",
        "assign_stability_c_f",
        "assign_stability_f_c",
        "assign_stability_state_c_f",
        "character_from_c",
        "character_to_c",
        "conv_check_f_wrapper",
        "convert_solver_settings_c",
        "copy_state_to_c",
        "get_extra_trial_vectors_f_wrapper",
        "hess_x_f_wrapper",
        "hess_x_f_wrapper_impl",
        "init_solver_settings_c",
        "init_stability_settings_c",
        "logger_f_wrapper",
        "obj_func_f_wrapper",
        "precond_f_wrapper",
        "precond_pd_f_wrapper",
        "project_f_wrapper",
        "solver_c_wrapper",
        "solver_settings_size",
        "solver_state_size",
        "stability_check_c_wrapper",
        "stability_hess_x_f_wrapper",
        "stability_settings_size",
        "stability_state_size",
        "store_optional_c_callbacks",
        "update_orbs_f_wrapper",
    ],
    "opentrustregion_system_tests": [
        "h2o_fb_solver_davidson_ah",
        "h2o_fb_solver_default",
        "h2o_fb_solver_gltr",
        "h2o_fb_solver_hess_symm",
        "h2o_fb_solver_jacobi_davidson_ah",
        "h2o_fb_solver_jacobi_davidson_ls",
        "h2o_fb_solver_line_search",
        "h2o_fb_solver_stability",
        "h2o_fb_solver_tcg",
        "h2o_fb_stability_check_default",
        "h2o_fb_stability_check_hess_symm",
        "h2o_fb_stability_check_jacobi_davidson",
    ],
    "c_integration_tests": [
        "settings_layout",
        "solver_c",
        "solver_settings_init",
        "stability_check_c",
        "stability_settings_init",
        "state_layout",
    ],
}

# define signatures of the Fortran functions providing the reference settings
lib.reference_field.restype = None
lib.get_reference_solver_values.argtypes = [
    POINTER(SolverSettings.c_struct),
    c_char_p,
]
lib.get_reference_solver_values.restype = None

# define signatures of the Fortran functions providing the reference states
lib.reference_state_field.argtypes = [c_char_p]
lib.reference_state_field.restype = c_real
lib.get_reference_solver_state_values.argtypes = [POINTER(SolverState.c_struct)]
lib.get_reference_solver_state_values.restype = None
lib.get_reference_stability_state_values.argtypes = [POINTER(StabilityState.c_struct)]
lib.get_reference_stability_state_values.restype = None


def reference(name, field_type):
    """
    this function returns the reference value of a settings field by its name,
    prefixed by "stability_settings." for the nested settings, which also describe
    standalone stability check settings, from the Fortran test reference module
    """
    value = c_real()
    keyword = dict(SolverSettings.c_struct._fields_)["subsystem_solver"]()
    lib.reference_field(name.encode(), byref(value), keyword)
    if issubclass(field_type, Array):
        return keyword.value.decode()
    if field_type == c_bool:
        return bool(value.value)
    if field_type == c_int:
        return int(value.value)
    return value.value


def reference_state(name, field_type):
    """
    this function returns the reference value of a state field by its name, prefixed
    by "stability_state." for the nested state, which also describes a standalone
    stability check state, from the Fortran test reference module
    """
    value = lib.reference_state_field(name.encode())
    if field_type == c_bool:
        return bool(value)
    if field_type == c_int:
        return int(value)
    return value


def check_callback_interfaces(test_name, interfaces, read_only, input_array):
    """
    this function checks that interfaces to Python callback functions, given by name
    as functions calling the interface with a callback function and an exception
    collector, report an exception or an interruption raised by the callback function
    as an error and keep it, and that those named in read_only hand the callback
    function their input, the given array, read-only as its first argument, so that
    writing into it is reported as an error and leaves the input unchanged
    """
    test_passed = True

    def raising(error_type):
        def callback(*args):
            raise error_type("callback function failure")

        return callback

    for error_type in (ValueError, KeyboardInterrupt):
        for name, call in interfaces.items():
            exception = {}
            if call(raising(error_type), exception) != 1:
                print(
                    f" {test_name} failed: {error_type.__name__} raised by {name} not "
                    "reported as error."
                )
                test_passed = False
            if not isinstance(exception.get("exc"), error_type):
                print(
                    f" {test_name} failed: {error_type.__name__} raised by {name} not "
                    "kept."
                )
                test_passed = False

    def writing_into_input(array, *args):
        array += 1.0

    original = input_array.copy()
    for name in read_only:
        exception = {}
        error = interfaces[name](writing_into_input, exception)
        if error != 1 or not isinstance(exception.get("exc"), ValueError):
            print(
                f" {test_name} failed: Writing into the input of the {name} not "
                "reported as error."
            )
            test_passed = False
        if not np.array_equal(input_array, original):
            print(f" {test_name} failed: Input of the {name} changed.")
            test_passed = False
            input_array[...] = original

    return test_passed


# define function to add tests to test classes
def add_tests(cls):

    def create_test(func_name):
        def test(self):
            func = getattr(lib, func_name)
            func.restype = c_bool
            result = func()
            if result:
                print(f" {func_name} PASSED")
            self.assertTrue(result, f"{func_name} failed")

        return test

    for func_name in cls.tests:
        setattr(cls, f"test_{func_name}", create_test(f"test_{func_name}"))

    return cls


# number of parameters
n_param = c_int.in_dll(lib, "test_n_param").value


@add_tests
class OpenTrustRegionUnitTests(unittest.TestCase):
    """
    this class contains unit tests for opentrustregion
    """

    tests = fortran_tests["opentrustregion_unit_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for OpenTrustRegion...")
        return super().setUpClass()


@add_tests
class CInterfaceUnitTests(unittest.TestCase):
    """
    this class contains unit tests for the C interface
    """

    tests = fortran_tests["c_interface_unit_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for C interface...")
        return super().setUpClass()


class PyInterfaceUnitTests(unittest.TestCase):
    """
    this class contains unit tests for the Python interface
    """

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for Python interface...")

        return super().setUpClass()

    @staticmethod
    def _mock_precond(residual, mu, precond_residual):
        """
        this function is a mock function for the preconditioner function
        """
        precond_residual[:] = mu * residual

    @staticmethod
    def _mock_precond_pd(residual, precond_residual):
        """
        this function is a mock function for the positive-definite preconditioner
        function
        """
        precond_residual[:] = 3 * residual

    @staticmethod
    def _mock_project(vector):
        """
        this function is a mock function for the projection function
        """
        vector[:] = 2 * vector

    @staticmethod
    def _mock_hess_x(x, hess_x):
        """
        this function is a mock function for the Hessian linear transformation
        function
        """
        hess_x[:] = 4 * x

    @staticmethod
    def _mock_get_extra_trial_vectors(trial_vectors):
        """
        this function is a mock function for the extra trial vector function, which
        fills the i-th trial vector with i
        """
        trial_vectors[:] = np.arange(1, len(trial_vectors) + 1).reshape(-1, 1)

    @staticmethod
    def _mock_conv_check():
        """
        this function is a mock function for the convergence check function
        """
        return True

    @staticmethod
    def _raising_logger(message):
        """
        this function is a logging function that raises
        """
        raise ValueError("logging failure")

    # replace original library with mock library
    @patch("pyopentrustregion.python_interface.lib.solver", lib.mock_solver)
    def test_solver_py_interface(self):
        """
        this function tests the solver python interface
        """

        def mock_obj_func(kappa):
            """
            this function is a mock function for the objective function
            """
            return np.sum(kappa)

        def mock_update_orbs(kappa, grad, h_diag):
            """
            this function is a mock function for the orbital update function
            """
            func = np.sum(kappa)
            grad[:] = 2 * kappa
            h_diag[:] = 3 * kappa

            def hess_x(x, hess_x):
                hess_x[:] = 4 * x

            return func, hess_x

        # the callback functions record that they were called and the logging
        # functions record the messages, separately for the settings and the nested
        # settings so that their callback functions can be told apart
        called = set()

        def recording(name, func):
            def callback(*args):
                called.add(name)
                return func(*args)

            return callback

        messages = []
        stability_messages = []
        settings = SolverSettings()
        settings.precond = recording("precond", self._mock_precond)
        settings.precond_pd = recording("precond_pd", self._mock_precond_pd)
        settings.project = recording("project", self._mock_project)
        settings.get_extra_trial_vectors = recording(
            "get_extra_trial_vectors", self._mock_get_extra_trial_vectors
        )
        settings.conv_check = self._mock_conv_check
        settings.stability_hess_x = recording("stability_hess_x", self._mock_hess_x)
        settings.logger = messages.append
        settings.stability_settings.precond = recording(
            "stability_precond", self._mock_precond
        )
        settings.stability_settings.project = recording(
            "stability_project", self._mock_project
        )
        settings.stability_settings.approx_hess_x = recording(
            "stability_approx_hess_x", self._mock_hess_x
        )
        settings.stability_settings.get_extra_trial_vectors = recording(
            "stability_get_extra_trial_vectors", self._mock_get_extra_trial_vectors
        )
        settings.stability_settings.conv_check = recording(
            "stability_conv_check", self._mock_conv_check
        )
        settings.stability_settings.logger = stability_messages.append
        for field_info in settings.c_struct._fields_:
            field_name, field_type = field_info[:2]
            if (
                field_type == c_void_p
                or field_name == "initialized"
                or (isinstance(field_type, type) and issubclass(field_type, Structure))
            ):
                continue
            setattr(settings, field_name, reference(field_name, field_type))

        # set reference values for nested stability check settings
        for field_info in settings.stability_settings.c_struct._fields_:
            field_name, field_type = field_info[:2]
            if field_type == c_void_p or field_name == "initialized":
                continue
            setattr(
                settings.stability_settings,
                field_name,
                reference("stability_settings." + field_name, field_type),
            )

        # call solver python interface with optional arguments, the result of the mock
        # is cleared before and read right after the call, keep the C settings to check
        # that the calls leave them unchanged, and hand the settings the number of
        # parameters and the exception collector of a call this one interrupts, which
        # the call has to hand back when it returns
        c_settings = bytes(settings.settings_c)
        interface_flag = c_bool.in_dll(lib, "test_solver_interface")
        interface_flag.value = False
        state = SolverState()
        interrupted_exception: dict = {}
        settings.prepare_call(n_param + 1, interrupted_exception)
        solver(mock_update_orbs, mock_obj_func, n_param, settings, state=state)
        interface_passed = interface_flag.value
        call_restored = (
            settings._n_param == n_param + 1
            and settings._exception is interrupted_exception
            and settings.stability_settings._n_param == n_param + 1
            and settings.stability_settings._exception is interrupted_exception
        )
        if not call_restored:
            print(
                " test_solver_py_interface failed: Call data of the interrupted call "
                "not restored."
            )

        # check that the state is passed, which the mock fills with the reference
        # values
        state_passed = self._state_holds_reference(
            "test_solver_py_interface", state, ""
        )

        # check if every callback function was called through its own slot
        test_logger = "test" in messages and "test" in stability_messages
        if not test_logger:
            print(" test_solver_py_interface failed: Called logging function wrong.")
        callbacks_called = True
        for name in [
            "precond",
            "precond_pd",
            "project",
            "get_extra_trial_vectors",
            "stability_hess_x",
            "stability_precond",
            "stability_project",
            "stability_approx_hess_x",
            "stability_get_extra_trial_vectors",
            "stability_conv_check",
        ]:
            if name not in called:
                print(
                    f" test_solver_py_interface failed: Callback function {name} was "
                    "not called."
                )
                callbacks_called = False

        # a logging function that raises must not be silent, the exception is reported
        # as the cause of a warning once the solver has returned without an error
        settings.logger = self._raising_logger
        logger_error_reported = False
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                solver(mock_update_orbs, mock_obj_func, n_param, settings)
                logger_error_reported = any(
                    issubclass(w.category, RuntimeWarning)
                    and isinstance(w.message.__cause__, ValueError)
                    for w in caught
                )
            except RuntimeError:
                pass
        if not logger_error_reported:
            print(
                " test_solver_py_interface failed: Exception raised by logging "
                "function was not reported as a warning."
            )

        # a callback function that raises makes the solver fail with the exception as
        # the cause of the reported error
        def raising_obj_func(kappa):
            raise ValueError("objective function failure")

        settings.logger = messages.append
        callback_error_reported = False
        try:
            solver(mock_update_orbs, raising_obj_func, n_param, settings)
        except RuntimeError as e:
            callback_error_reported = isinstance(e.__cause__, ValueError)
        if not callback_error_reported:
            print(
                " test_solver_py_interface failed: Exception raised by objective "
                "function was not reported."
            )

        # the C settings are left unchanged by the calls and by assigning another
        # logging function
        c_settings_kept = bytes(settings.settings_c) == c_settings
        if not c_settings_kept:
            print(" test_solver_py_interface failed: C settings changed by the calls.")
        # an orbital update that raises makes the solver fail with the exception as the
        # cause of the reported error, which a library calling it right away shows
        update_orbs_failure = ValueError("orbital update failure")

        def raising_update_orbs(kappa, grad, h_diag):
            raise update_orbs_failure

        def library_calling_update_orbs(update_orbs_interface, *args):
            vector = np.zeros(n_param)
            vector_ptr = vector.ctypes.data_as(POINTER(c_real))
            return update_orbs_interface(
                vector_ptr,
                pointer(c_real()),
                vector_ptr,
                vector_ptr,
                pointer(hess_x_interface_type()),
                None,
            )

        update_orbs_error_reported = False
        with patch(
            "pyopentrustregion.python_interface.lib.solver",
            library_calling_update_orbs,
        ):
            try:
                solver(raising_update_orbs, mock_obj_func, n_param)
            except RuntimeError as e:
                update_orbs_error_reported = e.__cause__ is update_orbs_failure
        if not update_orbs_error_reported:
            print(
                " test_solver_py_interface failed: Exception raised by orbital update "
                "was not reported."
            )

        # the solver marks the library as innermost while it runs, around the callback
        # function marked within it, and raises an interruption of a callback function
        # as it is
        interruption = KeyboardInterrupt()
        marks = []

        def interrupted_obj_func(kappa):
            marks.append(interruption_state.stack)
            raise interruption

        raised = None
        try:
            solver(mock_update_orbs, interrupted_obj_func, n_param, settings)
        except KeyboardInterrupt as e:
            raised = e
        interruption_raised = raised is interruption and marks == [(False, True)]
        if not interruption_raised:
            print(
                " test_solver_py_interface failed: Library not marked as innermost or "
                "interruption of callback function not raised."
            )

        # call solver python interface without settings, which have to be passed as
        # none so that the library uses its default settings, the results of the mock
        # start opposite to the expected ones
        settings_flag = c_bool.in_dll(lib, "solver_settings_passed")
        settings_flag.value = True
        interface_flag.value = False
        solver(mock_update_orbs, mock_obj_func, n_param)
        no_settings_passed = interface_flag.value and not settings_flag.value
        if settings_flag.value:
            print(
                " test_solver_py_interface failed: Settings passed although none were "
                "given."
            )

        self.assertTrue(
            interface_passed
            and call_restored
            and state_passed
            and test_logger
            and callbacks_called
            and logger_error_reported
            and callback_error_reported
            and c_settings_kept
            and update_orbs_error_reported
            and interruption_raised
            and no_settings_passed,
            "test_solver_py_interface failed",
        )
        print(" test_solver_py_interface PASSED")

    # replace original library with mock library
    @patch(
        "pyopentrustregion.python_interface.lib.stability_check",
        lib.mock_stability_check,
    )
    def test_stability_check_py_interface(self):
        """
        this function tests the stability check python interface
        """
        h_diag = np.full(n_param, 3.0, dtype=np.float64)

        def mock_hess_x(x, hess_x):
            hess_x[:] = 4 * x

        # initialize settings object, the logging function records the messages
        messages = []
        settings = StabilitySettings()
        settings.precond = self._mock_precond
        settings.project = self._mock_project
        settings.approx_hess_x = mock_hess_x
        settings.get_extra_trial_vectors = self._mock_get_extra_trial_vectors
        settings.conv_check = self._mock_conv_check
        settings.logger = messages.append
        for field_info in settings.c_struct._fields_:
            field_name, field_type = field_info[:2]
            if field_type == c_void_p or field_name == "initialized":
                continue
            setattr(
                settings,
                field_name,
                reference("stability_settings." + field_name, field_type),
            )

        # allocate memory for descent direction
        kappa = np.zeros(n_param, dtype=np.float64)

        # call stability check python interface with optional arguments, the result of
        # the mock is cleared before and read right after the call, check that the
        # call leaves the C settings unchanged, and hand the settings the number of
        # parameters and the exception collector of a call this one interrupts, which
        # the call has to hand back when it returns
        interface_flag = c_bool.in_dll(lib, "test_stability_check_interface")
        interface_flag.value = False
        state = StabilityState()
        c_settings = bytes(settings.settings_c)
        interrupted_exception: dict = {}
        settings.prepare_call(n_param + 1, interrupted_exception)
        stable = stability_check(
            h_diag, mock_hess_x, n_param, settings, kappa=kappa, state=state
        )
        interface_passed = interface_flag.value
        c_settings_kept = bytes(settings.settings_c) == c_settings
        if not c_settings_kept:
            print(
                " test_stability_check_py_interface failed: C settings changed by the "
                "call."
            )
        call_restored = (
            settings._n_param == n_param + 1
            and settings._exception is interrupted_exception
        )
        if not call_restored:
            print(
                " test_stability_check_py_interface failed: Call data of the "
                "interrupted call not restored."
            )

        # check that the state is passed, which the mock fills with the reference
        # values
        state_passed = self._state_holds_reference(
            "test_stability_check_py_interface", state, "stability_state."
        )

        # check if logger was called correctly
        test_logger = "test" in messages
        if not test_logger:
            print(
                " test_stability_check_py_interface failed: Called logging function "
                "wrong."
            )

        # check if returned variables are correct
        if not stable:
            print(
                " test_stability_check_py_interface failed: Returned stability boolean "
                "wrong."
            )
        wrong_direction = not np.allclose(
            kappa, np.full(n_param, 1.0, dtype=np.float64)
        )
        if wrong_direction:
            print(
                " test_stability_check_py_interface failed: Returned direction wrong."
            )

        # a logging function that raises must not be silent, the exception is reported
        # as the cause of a warning once the stability check has returned without an
        # error, which still returns its result
        settings.logger = self._raising_logger
        logger_error_reported = False
        stable_with_failed_logger = False
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                stable_with_failed_logger = stability_check(
                    h_diag, mock_hess_x, n_param, settings, kappa=kappa
                )
                logger_error_reported = any(
                    issubclass(w.category, RuntimeWarning)
                    and isinstance(w.message.__cause__, ValueError)
                    for w in caught
                )
            except RuntimeError:
                pass
        if not logger_error_reported:
            print(
                " test_stability_check_py_interface failed: Exception raised by "
                "logging function was not reported as a warning."
            )
        if not stable_with_failed_logger:
            print(
                " test_stability_check_py_interface failed: Returned stability boolean "
                "wrong with failed logging function."
            )

        # call stability check python interface without a returned direction
        settings.logger = messages.append
        stable_without_direction = stability_check(
            h_diag, mock_hess_x, n_param, settings
        )
        if not stable_without_direction:
            print(
                " test_stability_check_py_interface failed: Returned stability boolean "
                "wrong without returned direction."
            )

        # a callback function that raises makes the stability check fail with the
        # exception as the cause of the reported error
        def raising_hess_x(x, hess_x):
            raise ValueError("Hessian linear transformation failure")

        callback_error_reported = False
        try:
            stability_check(h_diag, raising_hess_x, n_param, settings)
        except RuntimeError as e:
            callback_error_reported = isinstance(e.__cause__, ValueError)
        if not callback_error_reported:
            print(
                " test_stability_check_py_interface failed: Exception raised by "
                "Hessian linear transformation was not reported."
            )

        # the stability check marks the library as innermost while it runs, around the
        # callback function marked within it, and raises an interruption of a callback
        # function as it is
        interruption = KeyboardInterrupt()
        marks = []

        def interrupted_hess_x(x, hess_x):
            marks.append(interruption_state.stack)
            raise interruption

        raised = None
        try:
            stability_check(h_diag, interrupted_hess_x, n_param, settings)
        except KeyboardInterrupt as e:
            raised = e
        interruption_raised = raised is interruption and marks == [(False, True)]
        if not interruption_raised:
            print(
                " test_stability_check_py_interface failed: Library not marked as "
                "innermost or interruption of callback function not raised."
            )

        # a Hessian diagonal given as a strided view, of another type or as a list is
        # passed with its values, while a Hessian diagonal of another length and an
        # array for the direction which the library can not write into are refused
        matrix = np.zeros((n_param, n_param))
        np.fill_diagonal(matrix, 3.0)
        h_diag_converted = True
        for h_diag_given, case in [
            (np.diag(matrix), "a strided view"),
            (np.full(n_param, 3.0, dtype=np.float32), "a single precision array"),
            ([3.0] * n_param, "a list"),
        ]:
            interface_flag.value = False
            stability_check(h_diag_given, mock_hess_x, n_param, settings)
            if not interface_flag.value:
                print(
                    " test_stability_check_py_interface failed: Hessian diagonal "
                    f"given as {case} not passed correctly."
                )
                h_diag_converted = False
        read_only = np.zeros(n_param)
        read_only.flags.writeable = False
        wrong_inputs_refused = True
        for h_diag_given, kappa_given, case in [
            (np.full(n_param + 1, 3.0), None, "Hessian diagonal of another length"),
            (h_diag, np.zeros(n_param + 1), "Direction of another length"),
            (h_diag, np.zeros(n_param, dtype=np.float32), "Single precision direction"),
            (h_diag, np.zeros(2 * n_param)[::2], "Strided direction"),
            (h_diag, read_only, "Read-only direction"),
            (h_diag, [0.0] * n_param, "Direction given as a list"),
        ]:
            try:
                stability_check(
                    h_diag_given, mock_hess_x, n_param, settings, kappa=kappa_given
                )
                print(f" test_stability_check_py_interface failed: {case} not refused.")
                wrong_inputs_refused = False
            except ValueError:
                pass

        # call stability check python interface without settings, which have to be
        # passed as none so that the library uses its default settings, the results of
        # the mock start opposite to the expected ones
        settings_flag = c_bool.in_dll(lib, "stability_settings_passed")
        settings_flag.value = True
        interface_flag.value = False
        stable_without_settings = stability_check(h_diag, mock_hess_x, n_param)
        no_settings_passed = (
            interface_flag.value and not settings_flag.value and stable_without_settings
        )
        if settings_flag.value:
            print(
                " test_stability_check_py_interface failed: Settings passed although "
                "none were given."
            )

        self.assertTrue(
            interface_passed
            and c_settings_kept
            and call_restored
            and state_passed
            and test_logger
            and stable
            and not wrong_direction
            and logger_error_reported
            and stable_with_failed_logger
            and stable_without_direction
            and callback_error_reported
            and interruption_raised
            and h_diag_converted
            and wrong_inputs_refused
            and no_settings_passed,
            "test_stability_check_py_interface failed",
        )
        print(" test_stability_check_py_interface PASSED")

    @staticmethod
    def _settings_hold_reference(test_name, settings, prefix, location=""):
        """
        this function checks that the optional callback functions of a settings object
        are unset and that every other field holds its reference value, the names of
        the stability check settings carry the prefix "stability_settings.", the host
        context and the initialization flag are skipped since the Python wrapper
        exposes neither
        """
        test_passed = True
        for field_info in settings.c_struct._fields_:
            field_name, field_type = field_info[:2]
            if field_name in ("context", "initialized"):
                continue
            if field_type == c_void_p:
                if getattr(settings, field_name) is not None:
                    print(
                        f" {test_name} failed: Optional function pointer {field_name} "
                        f"not initialized correctly{location}."
                    )
                    test_passed = False
            elif isinstance(field_type, type) and issubclass(field_type, Structure):
                test_passed &= PyInterfaceUnitTests._settings_hold_reference(
                    test_name,
                    getattr(settings, field_name),
                    "stability_settings.",
                    " for nested stability settings",
                )
            else:
                ref_value = reference(prefix + field_name, field_type)
                value = getattr(settings, field_name)
                if field_type == c_real:
                    match = np.isclose(value, ref_value)
                else:
                    match = value == ref_value
                if not match:
                    print(
                        f" {test_name} failed: Field {field_name} not initialized "
                        f"correctly{location}."
                    )
                    test_passed = False
        return test_passed

    @staticmethod
    def _settings_values(settings):
        """
        this function returns the bytes of the C struct of a settings object without
        the C callback functions, which every settings object builds for itself, also
        those of the nested stability check settings
        """
        settings_c = type(settings.settings_c).from_buffer_copy(settings.settings_c)
        structs = [(settings_c, settings)]
        if hasattr(settings, "stability_settings"):
            structs.append((settings_c.stability_settings, settings.stability_settings))
        for struct, owner in structs:
            for name in owner.callback_fields:
                setattr(struct, name, None)
        return bytes(settings_c)

    def _settings_safeguarded(self, test_name, settings):
        """
        this function checks that a settings object accepts no attribute that is not a
        setting or a callback function and that its copies hold the settings in a C
        struct of their own, a copy sharing the callback functions and a deep copy
        copying them, also for the nested stability check settings
        """
        test_passed = True

        # neither the initialization flag nor the host context nor a misspelled
        # setting can be set
        for name in ["initialized", "context", "conv_tl"]:
            try:
                setattr(settings, name, 0)
                print(f" {test_name} failed: Attribute {name} could be set.")
                test_passed = False
            except AttributeError:
                pass

        # change a setting, so that only a copy of the values holds it, and set a
        # callback function that a deep copy copies, on the nested settings too
        all_settings = [settings]
        if hasattr(settings, "stability_settings"):
            all_settings.append(settings.stability_settings)
        for some_settings in all_settings:
            some_settings.conv_tol = 3 * some_settings.conv_tol
            some_settings.precond = partial(self._mock_precond)

        for copy_function, shared, case in [
            (copy.copy, True, "copy"),
            (copy.deepcopy, False, "deep copy"),
        ]:
            settings_copy = copy_function(settings)
            copies = [settings_copy]
            if hasattr(settings_copy, "stability_settings"):
                copies.append(settings_copy.stability_settings)
            for original, duplicate in zip(all_settings, copies):
                # the copy holds the same settings, a change of which leaves the
                # original unchanged
                if self._settings_values(duplicate) != self._settings_values(original):
                    print(f" {test_name} failed: The {case} holds other settings.")
                    test_passed = False
                conv_tol = original.conv_tol
                duplicate.conv_tol = 2 * conv_tol
                if original.conv_tol != conv_tol:
                    print(
                        f" {test_name} failed: Changing the {case} changed the "
                        "original."
                    )
                    test_passed = False

                # a copy shares the callback functions, a deep copy copies them
                if (duplicate.precond is original.precond) != shared or not callable(
                    duplicate.precond
                ):
                    print(
                        f" {test_name} failed: Callback function wrong for the "
                        f"{case}."
                    )
                    test_passed = False

        # settings with a context, which only the wiring of a factory sets, whose object
        # calls the callback functions of these settings, can not be copied, also when
        # nested
        for some_settings in all_settings:
            some_settings.settings_c.context = 0x1000
            for copy_function, case in [
                (copy.copy, "copy"),
                (copy.deepcopy, "deep copy"),
            ]:
                try:
                    copy_function(some_settings)
                    print(f" {test_name} failed: Wired settings allow a {case}.")
                    test_passed = False
                except TypeError:
                    pass
            some_settings.settings_c.context = None

        return test_passed

    @staticmethod
    def _callbacks_wired(test_name, settings_class):
        """
        this function checks that the C callback functions of a settings object report
        an exception raised by the callable assigned to a callback function during a
        call as an error and keep it with the exceptions of the call, also for the
        nested stability check settings
        """
        test_passed = True
        n_param = 2
        vector = np.zeros(n_param)
        vector_ptr = vector.ctypes.data_as(POINTER(c_real))
        arguments = {
            "precond": (vector_ptr, pointer(c_real(1.0)), vector_ptr),
            "precond_pd": (vector_ptr, vector_ptr),
            "project": (vector_ptr,),
            "get_extra_trial_vectors": (vector_ptr, 1),
            "conv_check": (pointer(c_bool()),),
            "stability_hess_x": (vector_ptr, vector_ptr),
            "approx_hess_x": (vector_ptr, vector_ptr),
        }

        def raising(failure):
            def callback(*args):
                raise failure

            return callback

        # every callback function with an error code raises an exception of its own
        settings = settings_class()
        all_settings = [(settings, "")]
        if hasattr(settings, "stability_settings"):
            all_settings.append((settings.stability_settings, "stability_settings."))
        failures = {}
        for some_settings, prefix in all_settings:
            for name in arguments.keys() & some_settings.callback_fields.keys():
                failures[prefix + name] = ValueError(prefix + name)
                setattr(some_settings, name, raising(failures[prefix + name]))

        # hand the settings the exception collector of a call and call the C callback
        # function of each callback function
        exception = {}
        settings.prepare_call(n_param, exception)
        for some_settings, prefix in all_settings:
            for name in arguments.keys() & some_settings.callback_fields.keys():
                exception.clear()
                interface = cast(
                    c_void_p(getattr(some_settings.settings_c, name)),
                    some_settings.callback_fields[name][1],
                )
                if interface(*arguments[name], None) != 1 or (
                    exception.get("exc") is not failures[prefix + name]
                ):
                    print(
                        f" {test_name} failed: Exception raised by {prefix}{name} not "
                        "reported as error or not kept."
                    )
                    test_passed = False

        return test_passed

    @staticmethod
    def _state_holds_reference(test_name, state, prefix, location=""):
        """
        this function checks that every field of a state object holds its reference
        value, the names of the stability check state carry the prefix
        "stability_state."
        """
        test_passed = True
        for field_info in state.c_struct._fields_:
            field_name, field_type = field_info[:2]
            if isinstance(field_type, type) and issubclass(field_type, Structure):
                test_passed &= PyInterfaceUnitTests._state_holds_reference(
                    test_name,
                    getattr(state, field_name),
                    "stability_state.",
                    " for nested stability check state",
                )
            elif getattr(state, field_name) != reference_state(
                prefix + field_name, field_type
            ):
                print(
                    f" {test_name} failed: Field {field_name} not read correctly"
                    f"{location}."
                )
                test_passed = False
        return test_passed

    def _state_copies_kept(self, test_name, state, prefix):
        """
        this function checks that a copy and a deep copy of a state holding the
        reference values keep them when the original state is reset, since they hold
        them in a C struct of their own
        """
        test_passed = True
        copies = [copy.copy(state), copy.deepcopy(state)]
        memset(addressof(state.state_c), 0, sizeof(state.c_struct))
        for state_copy, case in zip(copies, ["copy", "deep copy"]):
            test_passed &= self._state_holds_reference(
                test_name, state_copy, prefix, f" for the {case} of a reset state"
            )
        return test_passed

    @staticmethod
    def _state_holds_defaults(test_name, state, location=""):
        """
        this function checks that every field of a state object holds its default
        value, zero or false
        """
        test_passed = True
        for field_info in state.c_struct._fields_:
            field_name, field_type = field_info[:2]
            if isinstance(field_type, type) and issubclass(field_type, Structure):
                test_passed &= PyInterfaceUnitTests._state_holds_defaults(
                    test_name,
                    getattr(state, field_name),
                    " for nested stability check state",
                )
            elif getattr(state, field_name) != 0:
                print(
                    f" {test_name} failed: Field {field_name} does not start from "
                    f"its default value{location}."
                )
                test_passed = False
        return test_passed

    def test_defer_interruption(self):
        """
        this function tests that the handler of SIGINT only records an interruption
        while the library runs and raises it otherwise
        """
        test_passed = True
        try:
            for stack, deferred, case in [
                ((False,), True, "while the library is innermost"),
                ((False, True, False), True, "while a nested call is innermost"),
                ((False, True), False, "while a callback function is innermost"),
                ((), False, "outside a call into the library"),
            ]:
                interruption_state.stack = stack
                interruption_state.pending = None
                raised = False
                try:
                    defer_interruption(signal.SIGINT, None)
                except KeyboardInterrupt:
                    raised = True
                if raised == deferred:
                    print(
                        " test_defer_interruption failed: Interruption "
                        f"{'raised' if raised else 'not raised'} {case}."
                    )
                    test_passed = False
                if (
                    isinstance(interruption_state.pending, KeyboardInterrupt)
                    != deferred
                ):
                    print(
                        " test_defer_interruption failed: Interruption "
                        f"{'not recorded' if deferred else 'recorded'} {case}."
                    )
                    test_passed = False
        finally:
            interruption_state.stack = ()
            interruption_state.pending = None

        self.assertTrue(test_passed, "test_defer_interruption failed")
        print(" test_defer_interruption PASSED")

    def test_callback_interfaces(self):
        """
        this function tests that the interfaces to the Python callback functions,
        those of the settings reached through the C callback functions the settings
        build, report an exception or an interruption as an error and keep it, also
        one of the interface itself and a call before the number of parameters is
        known, that they hand the callback functions their inputs read-only and mark
        a callback function as innermost while it runs, that the convergence check
        interfaces pass on the result and that the logging interface keeps its first
        exception and records an interruption
        """
        test_passed = True
        n_param = 2
        vector = np.zeros(n_param)
        vector_ptr = vector.ctypes.data_as(POINTER(c_real))
        output = np.zeros(n_param)
        output_ptr = output.ctypes.data_as(POINTER(c_real))

        # settings whose C callback functions are called below, kept alive until then,
        # since these belong to the settings and a garbage collection frees them
        # together with settings that nothing else refers to
        called_settings = []

        def settings_callback(settings, name, func, interface_type, exception):
            """
            this function assigns a callable to a callback function of settings,
            hands them the number of parameters and the exception collector of a call
            and returns their C callback function
            """
            setattr(settings, name, func)
            settings.prepare_call(n_param, exception)
            called_settings.append(settings)
            return interface_type(getattr(settings.settings_c, name))

        # every interface, called with a callback function and an exception collector,
        # with the vector as its input and the output as its output
        interfaces = {
            "objective function": lambda func, exc: ObjFuncInterface(
                func, n_param, exc
            )(vector_ptr, pointer(c_real()), None),
            "Hessian linear transformation": lambda func, exc: HessXInterface(
                func, n_param, exc
            )(vector_ptr, output_ptr, None),
            "orbital update": lambda func, exc: UpdateOrbsInterface(func, n_param, exc)(
                vector_ptr, pointer(c_real()), output_ptr, output_ptr, None, None
            ),
            "preconditioner": lambda func, exc: settings_callback(
                SolverSettings(), "precond", func, precond_interface_type, exc
            )(vector_ptr, pointer(c_real(1.0)), output_ptr, None),
            "positive-definite preconditioner": lambda func, exc: settings_callback(
                SolverSettings(), "precond_pd", func, precond_pd_interface_type, exc
            )(vector_ptr, output_ptr, None),
            "projection": lambda func, exc: settings_callback(
                SolverSettings(), "project", func, project_interface_type, exc
            )(vector_ptr, None),
            "extra trial vector function": lambda func, exc: settings_callback(
                SolverSettings(),
                "get_extra_trial_vectors",
                func,
                get_extra_trial_vectors_interface_type,
                exc,
            )(output_ptr, 1, None),
            "convergence check": lambda func, exc: settings_callback(
                SolverSettings(), "conv_check", func, conv_check_interface_type, exc
            )(pointer(c_bool()), None),
            "Hessian linear transformation of the stability check": (
                lambda func, exc: settings_callback(
                    SolverSettings(),
                    "stability_hess_x",
                    func,
                    hess_x_interface_type,
                    exc,
                )(vector_ptr, output_ptr, None)
            ),
        }

        # every interface reports an exception or an interruption as an error and keeps
        # it, and those whose callback function receives an input hand it over
        # read-only
        test_passed &= check_callback_interfaces(
            "test_callback_interfaces",
            interfaces,
            [
                "objective function",
                "Hessian linear transformation",
                "orbital update",
                "preconditioner",
                "positive-definite preconditioner",
                "Hessian linear transformation of the stability check",
            ],
            vector,
        )

        # a callback function of settings called before the number of parameters of a
        # call is known reports an error
        settings = SolverSettings()
        settings.project = lambda vector: None
        if project_interface_type(settings.settings_c.project)(vector_ptr, None) != 1:
            print(
                " test_callback_interfaces failed: Callback function called before "
                "the number of parameters is known does not report an error."
            )
            test_passed = False

        # the orbital update interface leaves the Hessian linear transformation unset
        # without an error when the orbital update does not provide one, so that the
        # solver reports the missing one, and reports one that is not callable as an
        # error of the orbital update
        for returned, expected_error, case in [
            (None, 0, "without Hessian linear transformation"),
            (
                np.zeros(n_param),
                1,
                "with Hessian linear transformation that is not callable",
            ),
        ]:
            exception = {}
            hess_x_funptr = pointer(hess_x_interface_type())
            error = UpdateOrbsInterface(
                lambda kappa, grad, h_diag: (0.0, returned), n_param, exception
            ).call(vector_ptr, pointer(c_real()), vector_ptr, vector_ptr, hess_x_funptr)
            if error != expected_error:
                print(
                    f" test_callback_interfaces failed: Error reported wrong for "
                    f"orbital update {case}."
                )
                test_passed = False
            if expected_error and not isinstance(exception.get("exc"), TypeError):
                print(
                    f" test_callback_interfaces failed: Exception not kept for "
                    f"orbital update {case}."
                )
                test_passed = False
            if hess_x_funptr[0]:
                print(
                    f" test_callback_interfaces failed: Hessian linear transformation "
                    f"set for orbital update {case}."
                )
                test_passed = False

        # the Hessian linear transformation an orbital update returns reports an
        # exception it raises as an error and keeps it with those of the orbital update
        def raising_hess_x(x, hess_x):
            raise ValueError("Hessian linear transformation failure")

        exception = {}
        hess_x_funptr = pointer(hess_x_interface_type())
        UpdateOrbsInterface(
            lambda kappa, grad, h_diag: (0.0, raising_hess_x), n_param, exception
        ).call(vector_ptr, pointer(c_real()), vector_ptr, vector_ptr, hess_x_funptr)
        if hess_x_funptr[0](vector_ptr, output_ptr, None) != 1 or not isinstance(
            exception.get("exc"), ValueError
        ):
            print(
                " test_callback_interfaces failed: Exception raised by Hessian linear "
                "transformation returned by orbital update not reported."
            )
            test_passed = False

        # the convergence check interface passes on that the optimization has not
        # converged
        converged = c_bool(True)
        conv_check = settings_callback(
            SolverSettings(), "conv_check", lambda: False, conv_check_interface_type, {}
        )
        if conv_check(pointer(converged), None) != 0:
            print(
                " test_callback_interfaces failed: Error reported for convergence "
                "check."
            )
            test_passed = False
        if converged.value:
            print(
                " test_callback_interfaces failed: Result of convergence check not "
                "passed on."
            )
            test_passed = False

        # the logging interface keeps the first of several exceptions
        failures = iter([ValueError("first"), ValueError("second")])

        def failing_logger(message):
            raise next(failures)

        exception = {}
        logger = settings_callback(
            SolverSettings(), "logger", failing_logger, logger_interface_type, exception
        )
        logger(b"test", None)
        logger(b"test", None)
        if str(exception.get("logger")) != "first":
            print(
                " test_callback_interfaces failed: First exception raised by "
                "logging function not kept."
            )
            test_passed = False

        # an exception raised by an interface itself, as for a pointer it can not
        # write the result through, is reported as an error and kept
        exception = {}
        if ObjFuncInterface(lambda kappa: 0.0, n_param, exception)(
            vector_ptr, None, None
        ) != 1 or not isinstance(exception.get("exc"), TypeError):
            print(
                " test_callback_interfaces failed: Exception raised by interface not "
                "reported as error."
            )
            test_passed = False

        # an interface marks a callback function as innermost while it runs and
        # unmarks it afterwards, also when it raises
        try:
            for error_type, case in [(None, "returns"), (ValueError, "raises")]:
                interruption_state.stack = (False,)
                marks = []

                def marking(kappa):
                    marks.append(interruption_state.stack)
                    if error_type is not None:
                        raise error_type("objective function failure")
                    return 0.0

                ObjFuncInterface(marking, n_param, {})(
                    vector_ptr, pointer(c_real()), None
                )
                if marks != [(False, True)] or interruption_state.stack != (False,):
                    print(
                        " test_callback_interfaces failed: Callback function not "
                        f"marked as innermost while it runs when it {case}."
                    )
                    test_passed = False
        finally:
            interruption_state.stack = ()

        # the logging interface records an interruption, which the next callback
        # function with an error code reports without being called
        try:
            interruption = KeyboardInterrupt()

            def interrupted_logger(message):
                raise interruption

            exception = {}
            settings_callback(
                SolverSettings(),
                "logger",
                interrupted_logger,
                logger_interface_type,
                exception,
            )(b"test", None)
            settings_callback(
                SolverSettings(),
                "logger",
                lambda message: None,
                logger_interface_type,
                exception,
            )(b"test", None)
            if interruption_state.pending is not interruption or "logger" in exception:
                print(
                    " test_callback_interfaces failed: Interruption of logging "
                    "function not recorded."
                )
                test_passed = False
            obj_func_calls = []
            error = ObjFuncInterface(obj_func_calls.append, n_param, exception)(
                vector_ptr, pointer(c_real()), None
            )
            if error != 1 or exception.get("exc") is not interruption or obj_func_calls:
                print(
                    " test_callback_interfaces failed: Recorded interruption not "
                    "reported as error by next callback function."
                )
                test_passed = False
            if interruption_state.pending is not None:
                print(
                    " test_callback_interfaces failed: Recorded interruption kept "
                    "after it was reported."
                )
                test_passed = False
        finally:
            interruption_state.pending = None

        self.assertTrue(test_passed, "test_callback_interfaces failed")
        print(" test_callback_interfaces PASSED")

    def test_c_buffer(self):
        """
        this function tests that an array a C function can read and write in place is
        handed over itself and every other array as a writeable, C-contiguous float64
        copy of its values
        """
        test_passed = True
        read_only = np.random.random(n_param)
        read_only.flags.writeable = False

        for array, in_place, case in [
            (
                np.random.random(n_param),
                True,
                "a writeable, C-contiguous float64 array",
            ),
            (np.random.random((n_param, 2))[:, 1], False, "a non-contiguous view"),
            (np.random.random(n_param).astype(np.float32), False, "a float32 array"),
            (read_only, False, "a read-only array"),
        ]:
            buffer = c_buffer(array)
            if (buffer is array) is not in_place:
                print(
                    f" test_c_buffer failed: Array {'not ' if in_place else ''}handed "
                    f"over itself for {case}."
                )
                test_passed = False
            if not (
                buffer.dtype == np.float64
                and buffer.flags.c_contiguous
                and buffer.flags.writeable
                and np.allclose(buffer, array)
            ):
                print(
                    f" test_c_buffer failed: Buffer not a writeable, C-contiguous "
                    f"float64 array of the values for {case}."
                )
                test_passed = False

        self.assertTrue(test_passed, "test_c_buffer failed")
        print(" test_c_buffer PASSED")

    def test_write_back(self):
        """
        this function tests that the values a C function wrote into the buffer of an
        output are written back into an array whose buffer is a copy, also of another
        kind, and that a read-only array refuses them
        """
        test_passed = True

        for array, case in [
            (np.random.random((n_param, 2))[:, 1], "a non-contiguous view"),
            (np.random.random(n_param).astype(np.float32), "a float32 array"),
            (np.random.random(n_param), "an array handed over itself"),
        ]:
            buffer = c_buffer(array)
            buffer[:] = np.arange(n_param)
            write_back(array, buffer)
            if not np.allclose(array, np.arange(n_param)):
                print(f" test_write_back failed: Values not written back for {case}.")
                test_passed = False

        read_only = np.random.random(n_param)
        read_only.flags.writeable = False
        try:
            write_back(read_only, c_buffer(read_only))
        except ValueError:
            pass
        else:
            print(" test_write_back failed: Values written into a read-only array.")
            test_passed = False

        self.assertTrue(test_passed, "test_write_back failed")
        print(" test_write_back PASSED")

    def test_raise_on_failure(self):
        """
        this function tests that an error returned by the library is raised with its
        code and the exception of a failed callback function as its cause, that a
        failed logging function is reported as a warning with its exception as the
        cause both with and without an error, that nothing is reported without either
        and that an interruption of a callback function is raised as it is
        """
        test_passed = True
        callback_failure = ValueError("callback function failure")
        logger_failure = ValueError("logging failure")

        for error, exception, case in [
            (7, {}, "for an error alone"),
            (
                7,
                {"exc": callback_failure},
                "for an error of a failed callback function",
            ),
            (7, {"logger": logger_failure}, "for an error and a failed logger"),
            (0, {"logger": logger_failure}, "for a failed logger alone"),
            (0, {}, "without failure"),
        ]:
            raised = None
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                try:
                    raise_on_failure(error, exception, "solver")
                except RuntimeError as e:
                    raised = e

            # an error is raised with its code and the exception of a failed callback
            # function as its cause, never the one of a failed logging function
            if error and (raised is None or "(code 7)" not in str(raised)):
                print(
                    f" test_raise_on_failure failed: Error not raised with its code "
                    f"{case}."
                )
                test_passed = False
            if raised is not None and raised.__cause__ is not exception.get("exc"):
                print(
                    f" test_raise_on_failure failed: Error raised with wrong cause "
                    f"{case}."
                )
                test_passed = False
            if not error and raised is not None:
                print(f" test_raise_on_failure failed: Raised without error {case}.")
                test_passed = False

            # a failed logging function is reported as a warning with its exception as
            # the cause and nothing else is warned about
            warned_causes = [
                w.message.__cause__
                for w in caught
                if issubclass(w.category, RuntimeWarning)
            ]
            if "logger" in exception and warned_causes != [logger_failure]:
                print(
                    f" test_raise_on_failure failed: Failed logging function not "
                    f"reported as a warning {case}."
                )
                test_passed = False
            if "logger" not in exception and caught:
                print(
                    f" test_raise_on_failure failed: Warned without failed logger {case}."
                )
                test_passed = False

        # an interruption of a callback function is raised as it is, rather than as an
        # error of the library
        interruption = KeyboardInterrupt()
        raised = None
        try:
            raise_on_failure(7, {"exc": interruption}, "solver")
        except BaseException as e:
            raised = e
        if raised is not interruption:
            print(
                " test_raise_on_failure failed: Interruption of callback function not "
                "raised as it is."
            )
            test_passed = False

        self.assertTrue(test_passed, "test_raise_on_failure failed")
        print(" test_raise_on_failure PASSED")

    def test_raise_callback_error(self):
        """
        this function tests that the error of a C function called from Python raises
        a RuntimeError, with the exception of a failing Python callback function it
        reached as its cause, while an interruption of such a function is raised as it
        is and no error raises nothing
        """
        test_passed = True
        failure = ValueError("callback function failure")
        interruption = KeyboardInterrupt()
        for error, exception, expected, case in [
            (0, {"exc": failure}, None, "without error"),
            (1, {"exc": failure}, failure, "with a failing callback function"),
            (1, None, None, "without exception collector"),
            (1, {"exc": interruption}, interruption, "with an interruption"),
        ]:
            raised = None
            try:
                raise_callback_error(error, exception, "C function raised error.")
            except BaseException as e:
                raised = e
            if expected is interruption:
                correct = raised is interruption
            elif error:
                correct = (
                    isinstance(raised, RuntimeError) and raised.__cause__ is expected
                )
            else:
                correct = raised is None
            if not correct:
                print(f" test_raise_callback_error failed: Wrong exception {case}.")
                test_passed = False

        self.assertTrue(test_passed, "test_raise_callback_error failed")
        print(" test_raise_callback_error PASSED")

    def test_call_library(self):
        """
        this function tests that a call into the library passes on the arguments and
        the error code and marks the library as innermost, that it replaces Python's
        default handler of SIGINT in the main thread only and restores it afterwards,
        also when the call raises, and that it raises an interruption that arrived
        during the call
        """
        test_passed = True
        seen = []

        def library(*args):
            seen.append(
                (args, interruption_state.stack, signal.getsignal(signal.SIGINT))
            )
            return 7

        def check_restored(case):
            if signal.getsignal(signal.SIGINT) is not signal.default_int_handler:
                print(f" test_call_library failed: Handler not restored {case}.")
                return False
            if interruption_state.stack != ():
                print(f" test_call_library failed: Marks not restored {case}.")
                return False
            return True

        original_handler = signal.signal(signal.SIGINT, signal.default_int_handler)
        try:
            # the outermost call from the main thread defers interruptions
            seen.clear()
            if call_library(library, 1, "a") != 7 or seen != [
                ((1, "a"), (False,), defer_interruption)
            ]:
                print(
                    " test_call_library failed: Outermost call from the main thread "
                    "not made correctly."
                )
                test_passed = False
            test_passed &= check_restored("after the outermost call")

            # a nested call keeps the handler and nests its mark
            seen.clear()
            call_library(lambda: call_library(library))
            if seen != [((), (False, False), defer_interruption)]:
                print(" test_call_library failed: Nested call not made correctly.")
                test_passed = False
            test_passed &= check_restored("after a nested call")

            # a call from another thread neither replaces the handler nor marks the
            # main thread
            seen.clear()
            thread = threading.Thread(target=call_library, args=(library,))
            thread.start()
            thread.join()
            if seen != [((), (False,), signal.default_int_handler)]:
                print(
                    " test_call_library failed: Call from another thread not made "
                    "correctly."
                )
                test_passed = False
            test_passed &= check_restored("after a call from another thread")

            # a call raises the exception of the library and restores the handler
            def failing_library():
                raise ValueError("library failure")

            try:
                call_library(failing_library)
                print(" test_call_library failed: Exception of library not raised.")
                test_passed = False
            except ValueError:
                pass
            test_passed &= check_restored("after a call that raises")

            # an interruption the library received is raised when the call returns
            interruption = KeyboardInterrupt()

            def interrupted_library():
                interruption_state.pending = interruption
                return 0

            raised = None
            try:
                call_library(interrupted_library)
            except KeyboardInterrupt as e:
                raised = e
            if raised is not interruption or interruption_state.pending is not None:
                print(
                    " test_call_library failed: Interruption received by the library "
                    "not raised when the call returns."
                )
                test_passed = False
            test_passed &= check_restored("after an interruption")

            # an interruption that arrives while the handler is restored, which Python
            # handles while it switches the handler, is raised when the call returns,
            # after the handler is restored
            real_signal = signal.signal

            def signal_handling_pending(signum, handler):
                if handler is signal.default_int_handler:
                    signal.getsignal(signum)(signum, None)
                return real_signal(signum, handler)

            raised = None
            with patch.object(signal, "signal", signal_handling_pending):
                try:
                    call_library(library)
                except KeyboardInterrupt as e:
                    raised = e
            if raised is None or interruption_state.pending is not None:
                print(
                    " test_call_library failed: Interruption that arrives while the "
                    "handler is restored not raised when the call returns."
                )
                test_passed = False
            test_passed &= check_restored("after an interruption while it is restored")

            # another handler is kept
            def custom_handler(signum, frame):
                pass

            signal.signal(signal.SIGINT, custom_handler)
            seen.clear()
            call_library(library)
            if seen[0][2] is not custom_handler:
                print(" test_call_library failed: Other handler replaced.")
                test_passed = False
            if signal.getsignal(signal.SIGINT) is not custom_handler:
                print(" test_call_library failed: Other handler not kept.")
                test_passed = False
        finally:
            signal.signal(signal.SIGINT, original_handler)
            interruption_state.stack = ()
            interruption_state.pending = None

        self.assertTrue(test_passed, "test_call_library failed")
        print(" test_call_library PASSED")

    def test_check_struct_sizes(self):
        """
        this function tests that the check of the structure sizes passes for a library
        whose structures have the sizes of the Python structures and raises an import
        error for one whose structure differs in size or which exports no sizes
        """
        test_passed = True
        structures = [
            ("solver_settings", SolverSettings.c_struct),
            ("stability_settings", StabilitySettings.c_struct),
            ("solver_state", SolverState.c_struct),
            ("stability_state", StabilityState.c_struct),
        ]
        sizes = {name: sizeof(structure) for name, structure in structures}

        def mock_library(sizes):
            """
            this function returns a mock library which exports the given sizes
            """
            return SimpleNamespace(
                _name="mock library",
                **{
                    f"{name}_size": (lambda size=size: size)
                    for name, size in sizes.items()
                },
            )

        # the check passes for the sizes of the Python structures
        try:
            check_struct_sizes(mock_library(sizes), structures)
        except ImportError:
            print(" test_check_struct_sizes failed: Raised for matching sizes.")
            test_passed = False

        # the check raises for a structure of another size and for a missing size
        for name in sizes:
            for wrong_sizes, case in [
                ({**sizes, name: sizes[name] + 8}, f"for a differing {name} size"),
                (
                    {key: size for key, size in sizes.items() if key != name},
                    f"for a missing {name} size",
                ),
            ]:
                try:
                    check_struct_sizes(mock_library(wrong_sizes), structures)
                    print(f" test_check_struct_sizes failed: Not raised {case}.")
                    test_passed = False
                except ImportError:
                    pass

        self.assertTrue(test_passed, "test_check_struct_sizes failed")
        print(" test_check_struct_sizes PASSED")

    @patch.object(SolverSettings, "init_c_struct", lib.mock_init_solver_settings)
    def test_solver_settings(self):
        """
        this function ensures the SolverSettings object is properly initialized and
        synchronized with the underlying C struct
        """
        settings = SolverSettings()
        test_passed = self._settings_hold_reference(
            "test_solver_settings", settings, ""
        )
        test_passed &= self._settings_safeguarded("test_solver_settings", settings)
        test_passed &= self._callbacks_wired("test_solver_settings", SolverSettings)

        # assigning a callable writes the C callback function the settings built for it
        # into the C settings and reading returns the callable
        calls = []

        def precond(residual, mu, precond_residual):
            calls.append(mu)

        settings.precond = precond
        c_funptr = settings.settings_c.precond
        if not c_funptr or settings.precond is not precond:
            print(
                " test_solver_settings failed: Assigned callback function not set "
                "correctly."
            )
            test_passed = False

        # the C callback function calls the callable with the number of parameters of
        # the call whatever context it receives, as for a callback function assigned
        # after a factory wrote its context into the settings
        vector = np.zeros(n_param)
        vector_ptr = vector.ctypes.data_as(POINTER(c_real))
        settings.prepare_call(n_param, {})
        error = precond_interface_type(c_funptr)(
            vector_ptr, pointer(c_real(2.0)), vector_ptr, 12345
        )
        if error != 0 or calls != [2.0]:
            print(
                " test_solver_settings failed: C callback function does not call the "
                "assigned callback function."
            )
            test_passed = False

        # assigning another callable keeps the C settings
        settings.precond = lambda residual, mu, precond_residual: None
        if settings.settings_c.precond != c_funptr:
            print(
                " test_solver_settings failed: C settings changed by assigning "
                "another callback function."
            )
            test_passed = False

        # an unset callback function leaves its C function pointer null
        settings.precond = None
        if settings.settings_c.precond is not None or settings.precond is not None:
            print(
                " test_solver_settings failed: Unset callback function does not "
                "leave a null function pointer."
            )
            test_passed = False

        # a C function wired into the settings or their nested settings is read as a
        # callable around it, which hands it the arguments, calls it with the context
        # of the settings it was read from while marking the library as innermost,
        # passes on its outputs, raises its error and keeps these settings alive
        settings.settings_c.context = 6789
        settings.settings_c.stability_settings.context = 4321
        x = np.random.random(n_param)
        n_vectors = 2
        received = []
        marks = []
        error_code = 0

        def as_vector(ptr):
            return np.ctypeslib.as_array(ptr, shape=(n_param,))

        def scale_into(factor):
            def body(x_ptr, out_ptr):
                as_vector(out_ptr)[:] = factor * as_vector(x_ptr)

            return body

        def precond_body(residual_ptr, mu_ptr, precond_residual_ptr):
            as_vector(precond_residual_ptr)[:] = mu_ptr[0] * as_vector(residual_ptr)

        def project_body(vector_ptr):
            as_vector(vector_ptr)[:] *= 5.0

        def trial_vectors_body(trial_vectors_ptr, n):
            trial_vectors = np.ctypeslib.as_array(trial_vectors_ptr, shape=(n, n_param))
            trial_vectors[:] = n

        def conv_check_body(converged_ptr):
            converged_ptr[0] = True

        def logger_body(message):
            received.append(string_at(message).decode("utf-8"))

        def project_call(func):
            vector = x.copy()
            func(vector)
            return vector

        def output_call(call, shape=(n_param,)):
            def read_and_call(func):
                out = np.zeros(shape)
                call(func, out)
                return out

            return read_and_call

        nested = settings.stability_settings
        cases = [
            (
                settings,
                "precond",
                precond_body,
                output_call(lambda func, out: func(x, 3.0, out)),
                3.0 * x,
            ),
            (
                settings,
                "precond_pd",
                scale_into(2.0),
                output_call(lambda func, out: func(x, out)),
                2.0 * x,
            ),
            (
                settings,
                "project",
                project_body,
                project_call,
                5.0 * x,
            ),
            (
                settings,
                "get_extra_trial_vectors",
                trial_vectors_body,
                output_call(lambda func, out: func(out), (n_vectors, n_param)),
                np.full((n_vectors, n_param), float(n_vectors)),
            ),
            (settings, "conv_check", conv_check_body, lambda func: func(), True),
            (
                settings,
                "stability_hess_x",
                scale_into(4.0),
                output_call(lambda func, out: func(x, out)),
                4.0 * x,
            ),
            (settings, "logger", logger_body, lambda func: func("message"), None),
            (
                nested,
                "precond",
                precond_body,
                output_call(lambda func, out: func(x, 3.0, out)),
                3.0 * x,
            ),
            (
                nested,
                "project",
                project_body,
                project_call,
                5.0 * x,
            ),
            (
                nested,
                "approx_hess_x",
                scale_into(6.0),
                output_call(lambda func, out: func(x, out)),
                6.0 * x,
            ),
            (
                nested,
                "get_extra_trial_vectors",
                trial_vectors_body,
                output_call(lambda func, out: func(out), (n_vectors, n_param)),
                np.full((n_vectors, n_param), float(n_vectors)),
            ),
            (nested, "conv_check", conv_check_body, lambda func: func(), True),
            (nested, "logger", logger_body, lambda func: func("message"), None),
        ]
        c_functions = []
        for case_settings, name, body, read_and_call, expected in cases:
            case = (
                f"{name} of the "
                f"{'nested ' if case_settings is nested else ''}settings"
            )
            interface_type = case_settings.callback_fields[name][1]

            def c_function(*args, body=body, name=name):
                received.append((name, args[-1]))
                marks.append(interruption_state.stack)
                body(*args[:-1])
                return None if name == "logger" else error_code

            c_functions.append(interface_type(c_function))
            setattr(
                case_settings.settings_c, name, cast(c_functions[-1], c_void_p).value
            )
            received.clear()
            marks.clear()
            error_code = 0
            func = getattr(case_settings, name)
            result = read_and_call(func)
            if not received or received[0] != (name, case_settings.settings_c.context):
                print(
                    f" test_solver_settings failed: C function wired as {case} not "
                    "called with the context of these settings."
                )
                test_passed = False
            if marks != [(False,)]:
                print(
                    f" test_solver_settings failed: Library not marked as innermost "
                    f"while the C function wired as {case} runs."
                )
                test_passed = False
            if name == "logger":
                passed_on = received[1:] == ["message"]
            else:
                passed_on = np.allclose(result, expected)
            if not passed_on or func._owner is not case_settings:
                print(
                    f" test_solver_settings failed: C function wired as {case} not "
                    "read correctly."
                )
                test_passed = False

            # the error of the C function is raised
            if name != "logger":
                error_code = 3
                try:
                    read_and_call(func)
                except RuntimeError:
                    pass
                else:
                    print(
                        f" test_solver_settings failed: Error of the C function "
                        f"wired as {case} not raised."
                    )
                    test_passed = False
            setattr(case_settings.settings_c, name, None)

        # the number of parameters and the exception collector of a call are held for
        # the calling thread, so that a call in another thread does not see them, and
        # restoring the call another call interrupted hands back its number of
        # parameters and exception collector, also to the nested stability check
        # settings
        interrupted_exception: dict = {}
        settings.prepare_call(n_param, interrupted_exception)
        seen_in_thread = []
        thread = threading.Thread(
            target=lambda: seen_in_thread.append(settings._n_param)
        )
        thread.start()
        thread.join()
        if seen_in_thread != [None]:
            print(
                " test_solver_settings failed: Number of parameters of a call seen in "
                "another thread."
            )
            test_passed = False
        settings.restore_call(settings.prepare_call(n_param + 1, {}))
        if (
            settings._n_param != n_param
            or settings._exception is not interrupted_exception
            or settings.stability_settings._n_param != n_param
            or settings.stability_settings._exception is not interrupted_exception
        ):
            print(
                " test_solver_settings failed: Call data of the interrupted call not "
                "restored."
            )
            test_passed = False

        # a keyword given as bytes, as read from a C struct, is stored as well
        settings.subsystem_solver = b"davidson_ah"
        if settings.subsystem_solver != "davidson_ah":
            print(
                " test_solver_settings failed: Keyword given as bytes not stored "
                "correctly."
            )
            test_passed = False

        # the nested stability check settings keep the settings alive, which own the
        # extension objects wired into both, once nothing else refers to them
        stability_settings = settings.stability_settings
        settings_ref = weakref.ref(settings)
        del settings, func, cases
        gc.collect()
        if settings_ref() is None:
            print(
                " test_solver_settings failed: Nested stability check settings do not "
                "keep the settings alive."
            )
            test_passed = False

        self.assertTrue(test_passed, "test_solver_settings failed")
        print(" test_solver_settings PASSED")

    @patch.object(StabilitySettings, "init_c_struct", lib.mock_init_stability_settings)
    def test_stability_settings(self):
        """
        this function ensures the StabilitySettings object is properly initialized and
        synchronized with the underlying C struct
        """
        settings = StabilitySettings()
        test_passed = self._settings_hold_reference(
            "test_stability_settings", settings, "stability_settings."
        )
        test_passed &= self._settings_safeguarded("test_stability_settings", settings)
        test_passed &= self._callbacks_wired(
            "test_stability_settings", StabilitySettings
        )
        self.assertTrue(test_passed, "test_stability_settings failed")
        print(" test_stability_settings PASSED")

    def test_solver_state(self):
        """
        this function ensures the SolverState object starts from the default values
        and reads every field of the underlying C struct, including those of the
        nested stability check state, which it does not let be set
        """
        state = SolverState()
        test_passed = self._state_holds_defaults("test_solver_state", state)
        lib.get_reference_solver_state_values(byref(state.state_c))
        test_passed &= self._state_holds_reference("test_solver_state", state, "")
        test_passed &= self._state_copies_kept("test_solver_state", state, "")
        try:
            state.n_update_orbs = 1
            print(" test_solver_state failed: Field could be set.")
            test_passed = False
        except AttributeError:
            pass
        self.assertTrue(test_passed, "test_solver_state failed")
        print(" test_solver_state PASSED")

    def test_stability_state(self):
        """
        this function ensures the StabilityState object starts from the default values
        and reads every field of the underlying C struct, which it does not let be set
        """
        state = StabilityState()
        test_passed = self._state_holds_defaults("test_stability_state", state)
        lib.get_reference_stability_state_values(byref(state.state_c))
        test_passed &= self._state_holds_reference(
            "test_stability_state", state, "stability_state."
        )
        test_passed &= self._state_copies_kept(
            "test_stability_state", state, "stability_state."
        )
        try:
            state.n_hess_x = 1
            print(" test_stability_state failed: Field could be set.")
            test_passed = False
        except AttributeError:
            pass
        self.assertTrue(test_passed, "test_stability_state failed")
        print(" test_stability_state PASSED")


@add_tests
class OpenTrustRegionSystemTests(unittest.TestCase):
    """
    this class contains system tests for opentrustregion
    """

    tests = fortran_tests["opentrustregion_system_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running system tests for OpenTrustRegion...")
        test_data = Path(__file__).parent / "test_data"
        if not os.path.isdir(test_data):
            raise RuntimeError(
                "test_data directory does not exist in same directory as testsuite.py."
            )
        lib.set_test_data_path(str(test_data).encode("utf-8"))
        return super().setUpClass()


@add_tests
class CIntegrationTests(unittest.TestCase):
    """
    this class contains integration tests that drive the real library through the C
    header
    """

    tests = fortran_tests["c_integration_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running integration tests for C interface...")
        return super().setUpClass()


class PyIntegrationTests(unittest.TestCase):
    """
    this class contains integration tests that drive the real library (not the mock
    library) through the Python interface, with the Hartmann 6D problem as a small
    workload that exercises every callback, settings field and return value crossing
    the interface
    """

    @classmethod
    def setUpClass(cls):
        print_separator("Running integration tests for Python interface...")

        # read the Hartmann 6D problem from the library
        def read_array(name, *shape):
            size = int(np.prod(shape))
            array = (c_real * size).in_dll(lib, name)
            return np.frombuffer(array, dtype=np.dtype(c_real), count=size).reshape(
                shape, order="F"
            )

        cls.n_param = c_int.in_dll(lib, "hartmann6d_n_param").value
        cls.n_terms = c_int.in_dll(lib, "hartmann6d_n_terms").value
        cls.alpha = read_array("hartmann6d_alpha", cls.n_terms)
        cls.A = read_array("hartmann6d_A", cls.n_terms, cls.n_param)
        cls.P = read_array("hartmann6d_P", cls.n_terms, cls.n_param)
        cls.minimum1 = read_array("hartmann6d_minimum1", cls.n_param)
        cls.near_minimum = read_array("hartmann6d_near_minimum", cls.n_param)
        cls.saddle_point = read_array("hartmann6d_saddle_point", cls.n_param)

        return super().setUpClass()

    # Hartmann 6D primitives shared by the callbacks

    @classmethod
    def _exp_terms(cls, x):
        return np.exp(-np.sum(cls.A * (x - cls.P) ** 2, axis=1))

    @classmethod
    def _func(cls, x):
        return -np.dot(cls.alpha, cls._exp_terms(x))

    @classmethod
    def _grad(cls, x):
        e = cls._exp_terms(x)
        return np.array(
            [
                np.sum(2.0 * cls.alpha * cls.A[:, j] * (x[j] - cls.P[:, j]) * e)
                for j in range(cls.n_param)
            ]
        )

    @classmethod
    def _hess(cls, x):
        e = cls._exp_terms(x)
        H = np.zeros((cls.n_param, cls.n_param))
        for i in range(cls.n_param):
            H[i, i] = 2.0 * np.sum(
                cls.alpha
                * cls.A[:, i]
                * e
                * (1.0 - 2.0 * cls.A[:, i] * (x[i] - cls.P[:, i]) ** 2)
            )
            for j in range(i):
                H[i, j] = -4.0 * np.sum(
                    cls.alpha
                    * cls.A[:, i]
                    * cls.A[:, j]
                    * (x[i] - cls.P[:, i])
                    * (x[j] - cls.P[:, j])
                    * e
                )
                H[j, i] = H[i, j]
        return H

    # Tests

    def test_settings_layout(self):
        """
        this function checks that every field of the Python settings structures,
        including the nested stability settings, is read back under its own name from
        settings whose fields Fortran sets one by one to the reference values
        """
        test_passed = True

        # names of the logicals, prefixed for the nested settings
        logicals = [
            name
            for name, field_type in SolverSettings.c_struct._fields_
            if field_type == c_bool
        ] + [
            "stability_settings." + name
            for name, field_type in StabilitySettings.c_struct._fields_
            if field_type == c_bool
        ]

        def read(settings_c, name):
            for part in name.split("."):
                settings_c = getattr(settings_c, part)
            return settings_c

        # check that every logical is read back under its own name, only one is set at
        # a time so that swapped logicals can be told apart
        for true_logical in logicals:
            settings_c = SolverSettings.c_struct()
            lib.get_reference_solver_values(byref(settings_c), true_logical.encode())
            for name in logicals:
                if read(settings_c, name) != (name == true_logical):
                    print(f" test_settings_layout failed: Field {name} misplaced.")
                    test_passed = False

        # check that every other field is read back under its own name
        settings_c = SolverSettings.c_struct()
        lib.get_reference_solver_values(byref(settings_c), None)
        for struct, prefix in [
            (settings_c, ""),
            (settings_c.stability_settings, "stability_settings."),
        ]:
            for name, field_type in struct._fields_:
                if field_type == c_bool or field_type == StabilitySettings.c_struct:
                    continue
                value = getattr(struct, name)
                if field_type == c_void_p:
                    value = value or 0
                elif issubclass(field_type, Array):
                    value = value.decode()
                if value != reference(prefix + name, field_type):
                    print(
                        f" test_settings_layout failed: Field {prefix + name} "
                        "misplaced."
                    )
                    test_passed = False
        self.assertTrue(test_passed, "test_settings_layout failed")
        print(" test_settings_layout PASSED")

    def test_state_layout(self):
        """
        this function checks that every field of the Python state structures,
        including the nested stability check state, is read back under its own name
        from a state whose fields Fortran sets one by one to the reference values
        """
        test_passed = True

        # check that every field is read back under its own name
        state_c = SolverState.c_struct()
        lib.get_reference_solver_state_values(byref(state_c))
        for struct, prefix in [
            (state_c, ""),
            (state_c.stability_state, "stability_state."),
        ]:
            for name, field_type in struct._fields_:
                if field_type == StabilityState.c_struct:
                    continue
                if getattr(struct, name) != reference_state(prefix + name, field_type):
                    print(
                        f" test_state_layout failed: Field {prefix + name} misplaced."
                    )
                    test_passed = False
        self.assertTrue(test_passed, "test_state_layout failed")
        print(" test_state_layout PASSED")

    def test_solver_py(self):
        """
        this function drives solver() end-to-end from near a minimum with every
        optional callback set
        """
        test_passed = True

        # mutable closure state holding the current point, starting in the quadratic
        # region near the first minimum, its Hessian and the names of the callbacks
        # that were reached
        state = {
            "curr": self.near_minimum.copy(),
            "hess": None,
            "called": set(),
            "n_update_orbs_calls": 0,
            "state_not_current": False,
            "nested_n_hess_x_seen": None,
            "residual_rms_seen": None,
        }

        # the solver state record, which the callbacks read during the solve
        solver_state = SolverState()

        def update_orbs(delta_vars, grad, h_diag):
            # the state is copied before every callback, so it counts the updates
            # before this one
            if solver_state.n_update_orbs != state["n_update_orbs_calls"]:
                state["state_not_current"] = True
            state["n_update_orbs_calls"] += 1
            state["curr"] = state["curr"] + delta_vars
            x = state["curr"]
            state["hess"] = self._hess(x)
            grad[:] = self._grad(x)
            h_diag[:] = np.diag(state["hess"])

            def hess_x(v, hv):
                # the internal stability check calls this function too when the
                # settings provide no Hessian linear transformation of the stability
                # check, and its state then counts the calls before this one
                state["nested_n_hess_x_seen"] = solver_state.stability_state.n_hess_x
                hv[:] = state["hess"] @ v

            return self._func(x), hess_x

        def obj_func(delta_vars):
            return self._func(state["curr"] + delta_vars)

        # identity preconditioners and projections and loggers that record their
        # names, the stability variants are set only on the nested stability settings,
        # so the test can tell whether the internal stability check reached its own
        # callback slots rather than the solver's, the identity preconditioners
        # exercise the callback without producing a zero vector when mu=0 (which would
        # trip the Gram-Schmidt zero-vector guard)
        def recording_precond(name):
            def precond(residual, mu, out):
                state["called"].add(name)
                out[:] = residual

            return precond

        def recording_callback(name):
            def callback(*args):
                state["called"].add(name)

            return callback

        def conv_check():
            # report convergence by the gradient criterion with the default threshold
            # it replaces, so the solve is unaffected
            state["called"].add("conv_check")
            return (
                np.sqrt(np.mean(self._grad(state["curr"]) ** 2))
                < SolverSettings().conv_tol
            )

        # a positive-definite preconditioner for the truncated conjugate gradient and
        # generalized Lanczos subsystem solvers, which the Davidson subsystem solver of
        # this test does not use
        def precond_pd(residual, out):
            out[:] = residual

        # extra trial vector functions that contribute no direction, so that the solve
        # is unaffected
        def recording_get_extra_trial_vectors(name):
            def get_extra_trial_vectors(trial_vectors):
                state["called"].add(name)
                trial_vectors[:] = 0.0

            return get_extra_trial_vectors

        # the Hessian linear transformation of the internal stability check, which it
        # uses instead of the one the orbital update returns
        def stability_hess_x(v, hv):
            state["called"].add("stability_hess_x")
            state["nested_n_hess_x_seen"] = solver_state.stability_state.n_hess_x
            hv[:] = state["hess"] @ v

        # an approximate Hessian linear transformation, which the stability check only
        # uses for extra trial vectors when no extra trial vector function is provided
        def approx_hess_x(v, hv):
            state["called"].add("stability_approx_hess_x")
            hv[:] = state["hess"] @ v

        # convergence check of the internal stability check, which reads the RMS
        # residual of the lowest eigenpair from the nested state, records it and
        # accepts a converged residual
        def stability_conv_check():
            state["called"].add("stability_conv_check")
            state["residual_rms_seen"] = solver_state.stability_state.residual_rms
            return solver_state.stability_state.residual_rms < 1e-9

        settings = SolverSettings()
        settings.precond = recording_precond("precond")
        settings.precond_pd = precond_pd
        settings.project = recording_callback("project")
        settings.get_extra_trial_vectors = recording_get_extra_trial_vectors(
            "get_extra_trial_vectors"
        )
        settings.conv_check = conv_check
        settings.stability_hess_x = stability_hess_x
        settings.logger = recording_callback("logger")
        settings.stability = True
        settings.stability_settings.precond = recording_precond("stability_precond")
        settings.stability_settings.project = recording_callback("stability_project")
        settings.stability_settings.approx_hess_x = approx_hess_x
        settings.stability_settings.get_extra_trial_vectors = (
            recording_get_extra_trial_vectors("stability_get_extra_trial_vectors")
        )
        settings.stability_settings.conv_check = stability_conv_check
        settings.stability_settings.logger = recording_callback("stability_logger")
        settings.verbose = 3  # ensure the loggers are exercised

        # the maximum precision flag starts opposite to the one the converging solve
        # returns, so that its write-back is detected
        solver_state.state_c.max_precision_reached = True

        solver(update_orbs, obj_func, self.n_param, settings, state=solver_state)
        if solver_state.max_precision_reached:
            print(
                " test_solver_py failed: Maximum precision reached flag was not "
                "returned."
            )
            test_passed = False
        for name, description in [
            ("precond", "Preconditioner"),
            ("project", "Projection"),
            ("get_extra_trial_vectors", "Extra trial vector function"),
            ("conv_check", "Convergence check"),
            ("logger", "Logger"),
        ]:
            if name not in state["called"]:
                print(f" test_solver_py failed: {description} was not called.")
                test_passed = False
        if "stability_hess_x" not in state["called"]:
            print(
                " test_solver_py failed: Hessian linear transformation of the "
                "stability check was not called by the internal stability check."
            )
            test_passed = False
        if solver_state.n_update_orbs <= 0:
            print(" test_solver_py failed: Orbital update counter was not populated.")
            test_passed = False
        if solver_state.n_obj_func <= 0:
            print(
                " test_solver_py failed: Objective function counter was not populated."
            )
            test_passed = False
        if solver_state.n_hess_x <= 0:
            print(
                " test_solver_py failed: Hessian linear transformation counter was not "
                "populated."
            )
            test_passed = False
        if solver_state.stability_state.n_hess_x <= 0:
            print(
                " test_solver_py failed: Hessian linear transformation counter of the "
                "internal stability check was not populated."
            )
            test_passed = False

        # every orbital update read the updates before it, the last Hessian linear
        # transformation of the internal stability check the ones before it and its
        # last convergence check the residual of its last iteration
        if (
            state["state_not_current"]
            or state["n_update_orbs_calls"] != solver_state.n_update_orbs
            or state["nested_n_hess_x_seen"]
            != solver_state.stability_state.n_hess_x - 1
            or state["residual_rms_seen"] != solver_state.stability_state.residual_rms
        ):
            print(" test_solver_py failed: State was not current in the callbacks.")
            test_passed = False
        for name, description in [
            ("stability_precond", "Preconditioner"),
            ("stability_project", "Projection"),
            ("stability_get_extra_trial_vectors", "Extra trial vector function"),
            ("stability_conv_check", "Convergence check"),
            ("stability_logger", "Logger"),
        ]:
            if name not in state["called"]:
                print(
                    f" test_solver_py failed: {description} set on the nested "
                    "stability settings was not called by the internal stability check."
                )
                test_passed = False

        # solve again without an extra trial vector function of the nested settings,
        # which the internal stability check does not inherit, so that it reaches the
        # approximate Hessian linear transformation of the nested settings instead, and
        # without a Hessian linear transformation of the stability check, so that the
        # internal stability check calls the one the orbital update returns, which
        # finds the state of the internal stability check current
        state["curr"] = self.near_minimum.copy()
        state["called"].clear()
        state["nested_n_hess_x_seen"] = None
        settings.stability_settings.get_extra_trial_vectors = None
        settings.stability_hess_x = None
        solver(update_orbs, obj_func, self.n_param, settings, state=solver_state)
        if (
            "stability_hess_x" in state["called"]
            or state["nested_n_hess_x_seen"] is None
            or state["nested_n_hess_x_seen"]
            != solver_state.stability_state.n_hess_x - 1
        ):
            print(
                " test_solver_py failed: Hessian linear transformation of the orbital "
                "update was not called by the internal stability check with its state "
                "current."
            )
            test_passed = False
        if "stability_approx_hess_x" not in state["called"]:
            print(
                " test_solver_py failed: Approximate Hessian linear transformation set "
                "on the nested stability settings was not called by the internal "
                "stability check."
            )
            test_passed = False
        self.assertTrue(test_passed, "test_solver_py failed")
        print(" test_solver_py PASSED")

    def test_stability_check_py(self):
        """
        this function drives stability_check() end-to-end at a minimum and at a saddle
        point with every optional callback set
        """
        test_passed = True
        called = set()

        # identity preconditioner, identity projection and logger that record their
        # names, the identity preconditioner exercises the callback without producing
        # a zero vector when mu=0 (which would trip the Gram-Schmidt zero-vector guard)
        def precond(residual, mu, out):
            called.add("precond")
            out[:] = residual

        def project(vector):
            called.add("project")

        def logger(msg):
            called.add("logger")

        # an approximate Hessian linear transformation, which the stability check only
        # uses for extra trial vectors when no extra trial vector function is provided
        def approx_hess_x(v, hv):
            called.add("approx_hess_x")
            hv[:] = H @ v

        # an extra trial vector function that contributes no direction and a
        # convergence check that reads the RMS residual of the lowest eigenpair from the
        # state, records it and accepts a converged residual
        def get_extra_trial_vectors(trial_vectors):
            called.add("get_extra_trial_vectors")
            trial_vectors[:] = 0.0

        residual_rms_seen = []

        def conv_check():
            called.add("conv_check")
            residual_rms_seen.append(stability_state.residual_rms)
            return stability_state.residual_rms < 1e-9

        settings = StabilitySettings()
        settings.precond = precond
        settings.project = project
        settings.approx_hess_x = approx_hess_x
        settings.get_extra_trial_vectors = get_extra_trial_vectors
        settings.conv_check = conv_check
        settings.logger = logger
        settings.verbose = 3  # ensure logger is exercised

        # at a minimum, expect stable
        H = self._hess(self.minimum1)
        h_diag = np.diag(H).copy()

        def hess_x(v, hv):
            hv[:] = H @ v

        # the stability check state record, which the Hessian linear transformation
        # reads during the first call
        stability_state = StabilityState()
        n_hess_x_seen = []

        def watching_hess_x(v, hv):
            n_hess_x_seen.append(stability_state.n_hess_x)
            hess_x(v, hv)

        kappa = np.zeros(self.n_param)
        stable = stability_check(
            h_diag,
            watching_hess_x,
            self.n_param,
            settings,
            kappa=kappa,
            state=stability_state,
        )
        if not stable:
            print(
                " test_stability_check_py failed: Stability incorrectly classifies "
                "stability of minimum."
            )
            test_passed = False
        if stability_state.iteration <= 0:
            print(" test_stability_check_py failed: State was not returned.")
            test_passed = False

        # the last Hessian linear transformation read the ones before it and the last
        # convergence check the residual of the last iteration
        if (
            not n_hess_x_seen
            or n_hess_x_seen[-1] != stability_state.n_hess_x - 1
            or not residual_rms_seen
            or residual_rms_seen[-1] != stability_state.residual_rms
        ):
            print(
                " test_stability_check_py failed: State was not current in the "
                "callbacks."
            )
            test_passed = False
        for name, description in [
            ("precond", "Preconditioner"),
            ("project", "Projection"),
            ("get_extra_trial_vectors", "Extra trial vector function"),
            ("conv_check", "Convergence check"),
            ("logger", "Logger"),
        ]:
            if name not in called:
                print(f" test_stability_check_py failed: {description} was not called.")
                test_passed = False
        if stability_state.n_hess_x <= 0:
            print(
                " test_stability_check_py failed: Hessian linear transformation "
                "counter was not populated."
            )
            test_passed = False
        if not stability_state.min_eigval > 0.0:
            print(
                " test_stability_check_py failed: State does not hold a positive "
                "lowest eigenvalue for minimum."
            )
            test_passed = False

        # also exercise the call without a direction at the minimum, where the
        # stability flag differs from the False the wrapper starts from, without extra
        # trial vector function so that the extra trial vectors are taken from the
        # approximate Hessian linear transformation, with the state, which the
        # convergence check reads
        settings.get_extra_trial_vectors = None
        called.discard("approx_hess_x")
        stable = stability_check(
            h_diag, hess_x, self.n_param, settings, state=stability_state
        )
        if not stable:
            print(
                " test_stability_check_py failed: Stability incorrectly classifies "
                "stability of minimum when not passing direction."
            )
            test_passed = False
        if "approx_hess_x" not in called:
            print(
                " test_stability_check_py failed: Approximate Hessian linear "
                "transformation was not called without extra trial vector function."
            )
            test_passed = False

        # at a saddle, the stability flag and the sign of the lowest eigenvalue have
        # to differ from the ones returned at the minimum and the returned direction
        # has to be a normalized direction of negative curvature, with the extra trial
        # vector function again, so that the extra trial vectors are not taken from the
        # approximate Hessian linear transformation, which is exact here
        H = self._hess(self.saddle_point)
        h_diag = np.diag(H).copy()
        settings.get_extra_trial_vectors = get_extra_trial_vectors
        stable = stability_check(
            h_diag, hess_x, self.n_param, settings, kappa=kappa, state=stability_state
        )
        if stable:
            print(
                " test_stability_check_py failed: Stability incorrectly classifies "
                "stability of saddle point."
            )
            test_passed = False
        if not stability_state.min_eigval < 0.0:
            print(
                " test_stability_check_py failed: State does not hold a negative "
                "lowest eigenvalue for saddle point."
            )
            test_passed = False
        if not np.isclose(np.linalg.norm(kappa), 1.0, rtol=0.0, atol=1e-6):
            print(
                " test_stability_check_py failed: Stability check does not return a "
                "normalized direction for saddle point."
            )
            test_passed = False
        if not kappa @ H @ kappa < 0.0:
            print(
                " test_stability_check_py failed: Stability check does not return a "
                "direction of negative curvature for saddle point."
            )
            test_passed = False
        self.assertTrue(test_passed, "test_stability_check_py failed")
        print(" test_stability_check_py PASSED")
