# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import annotations

import os
import sys
import unittest
from importlib import resources
from ctypes import (
    CDLL,
    c_bool,
    c_char,
    c_void_p,
    Array,
    byref,
    CFUNCTYPE,
    Structure,
    POINTER,
)
from unittest.mock import patch
from pathlib import Path
from typing import TYPE_CHECKING

# check if numpy is available
try:
    import numpy as np

    NUMPY_AVAILABLE = True
except:
    NUMPY_AVAILABLE = False


# check if pyopentrustregion is installed or import module in same directory
try:
    from pyopentrustregion import (
        SolverSettings,
        StabilitySettings,
        solver,
        stability_check,
        c_int,
        c_real,
    )

    if TYPE_CHECKING:
        from pyopentrustregion.python_interface import Settings
except ImportError:
    sys.path.insert(0, str(Path(__file__).parent.absolute()))
    from pyopentrustregion import (
        SolverSettings,
        StabilitySettings,
        solver,
        stability_check,
        c_int,
        c_real,
    )

    if TYPE_CHECKING:
        from pyopentrustregion.python_interface import Settings

ext = "dylib" if sys.platform == "darwin" else "so"


def load_testsuite_library() -> CDLL:
    # try to load from installed package (site-packages)
    lib_path = resources.files("pyopentrustregion") / f"libtestsuite.{ext}"
    if lib_path.is_file():
        return CDLL(str(lib_path))

    # fallback: try to load from same directory (editable install)
    local_path = os.path.join(os.path.dirname(__file__), f"libtestsuite.{ext}")
    if os.path.exists(local_path):
        return CDLL(local_path)

    # fallback: try to load from ../build (development build)
    build_path = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../build", f"libtestsuite.{ext}")
    )
    if os.path.exists(build_path):
        return CDLL(build_path)

    # if all failed
    raise FileNotFoundError(
        f"Cannot find any of the expected libraries: libtestsuite.{ext}"
    )


# load library
lib = load_testsuite_library()


def print_separator(string: str):
    sep_length = 70
    print(sep_length * "-")
    print(string)
    print(sep_length * "-")


# define all tests in alphabetical order
fortran_tests = {
    "opentrustregion_tests": [
        "accept_trust_region_step",
        "add_column",
        "add_error_origin",
        "bisection",
        "bracket",
        "extend_matrix",
        "general_mat_diag",
        "general_mat_min_eig",
        "generalized_lanczos_trust_region",
        "generate_random_trial_vectors",
        "generate_trial_vectors",
        "get_precond_level_shift",
        "gram_schmidt",
        "init_rng",
        "init_solver_settings",
        "init_stability_settings",
        "jacobi_davidson_correction",
        "level_shifted_davidson",
        "level_shifted_diag_precond",
        "print_message",
        "minres",
        "newton_step",
        "orthogonal_projection",
        "print_results",
        "solver",
        "solver_sanity_check",
        "split_string_by_space",
        "stability_check",
        "stability_sanity_check",
        "string_to_lowercase",
        "symm_mat_diag",
        "symm_mat_min_eig",
        "truncated_conjugate_gradient",
    ],
    "c_interface_tests": [
        "assign_solver_c_f",
        "assign_solver_f_c",
        "assign_stability_c_f",
        "assign_stability_f_c",
        "character_from_c",
        "character_to_c",
        "conv_check_f_wrapper",
        "hess_x_f_wrapper",
        "init_solver_settings_c",
        "init_stability_settings_c",
        "logger_f_wrapper",
        "modify_step_f_wrapper",
        "obj_func_f_wrapper",
        "precond_f_wrapper",
        "project_f_wrapper",
        "solver_c_wrapper",
        "stability_check_c_wrapper",
        "update_orbs_f_wrapper",
    ],
    "system_tests": ["h2o_atomic_fb", "h2o_saddle_fb"],
}

# define return type of Fortran functions
for tests in fortran_tests.values():
    for test in tests:
        getattr(lib, f"test_{test}").restype = c_bool


# define function to add tests to test classes
def add_tests(cls):

    def create_test(func_name):
        def test(self):
            result = getattr(lib, func_name)()
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
class OpenTrustRegionTests(unittest.TestCase):
    """
    this class contains unit tests for opentrustregion
    """

    tests = fortran_tests["opentrustregion_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for OpenTrustRegion...")
        return super().setUpClass()


@add_tests
class CInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the C interface
    """

    tests = fortran_tests["c_interface_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for C interface...")
        return super().setUpClass()


@unittest.skipUnless(NUMPY_AVAILABLE, "NumPy not available.")
class PyInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the Python interface
    """

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for Python interface...")

        # get combined fields
        combined_fields = (
            SolverSettings.c_struct._fields_ + StabilitySettings.c_struct._fields_
        )

        # get non-duplicate fields and allocate memory for reference values
        ref_fields = []
        seen_names = set()
        # loop by type order
        for curr_type in (c_bool, c_real, c_int):
            # solver fields of this type
            for name, t in combined_fields:
                if t == curr_type and name not in seen_names and name != "initialized":
                    ref_fields.append((name + "_ref", t))
                    seen_names.add(name)

        # handle fixed-size c_char arrays (strings)
        for name, t in combined_fields:
            if (
                issubclass(t, Array)
                and getattr(t, "_type_", None) is c_char
                and name not in seen_names
            ):
                ref_fields.append((name + "_ref", t))
                seen_names.add(name)

        # create class to read reference values
        class RefSettingsC(Structure):
            _fields_ = ref_fields

        # create instance
        ref_settings = RefSettingsC()

        # call Fortran to fill values
        lib.get_reference_values.argtypes = [POINTER(RefSettingsC)]
        lib.get_reference_values.restype = None
        lib.get_reference_values(byref(ref_settings))

        # extract Python values
        for name, _ in ref_fields:
            ref_value = getattr(ref_settings, name)
            if isinstance(ref_value, bytes):
                ref_value = ref_value.decode("utf-8")
            setattr(cls, name, ref_value)

        return super().setUpClass()

    def assign_ref_to_settings(self, settings: Settings):
        for field_info in settings.c_struct._fields_:
            field_name, field_type = field_info[:2]
            if field_type == c_void_p:
                setattr(settings, field_name, getattr(self, "mock_" + field_name))
            elif field_name == "initialized":
                continue
            else:
                setattr(settings, field_name, getattr(self, field_name + "_ref"))

    def equal_settings_to_ref(self, settings: Settings) -> bool:
        test_passed = True
        for field_info in settings.c_struct._fields_:
            field_name, field_type = field_info[:2]
            if field_type == c_void_p:
                if (
                    getattr(settings, field_name) is not None
                    or getattr(settings.settings_c, field_name) is not None
                ):
                    print(
                        f" Optional function pointer {field_name} not initialized "
                        "correctly."
                    )
                    test_passed = False
            elif field_name == "initialized":
                if not getattr(settings, field_name):
                    print(" Field initialized not initialized correctly.")
                    test_passed = False
            else:
                ref_value = getattr(self, field_name + "_ref")
                if field_type == c_real:
                    match = np.isclose(getattr(settings, field_name), ref_value)
                else:
                    match = getattr(settings, field_name) == ref_value
                if not match:
                    print(field_name, getattr(settings, field_name), ref_value)
                    print(f" Field {field_name} not initialized correctly.")
                    test_passed = False
        return test_passed

    def mock_obj_func(self, kappa):
        """
        this function is a mock function for the objective function
        """
        return np.sum(kappa)

    def mock_hess_x(self, x, hess_x):
        """
        this function is a mock function for the Hessian linear transformation
        function
        """
        hess_x[:] = 4 * x

    def mock_update_orbs(self, kappa, grad, h_diag):
        """
        this function is a mock function for the orbital update function
        """
        func = np.sum(kappa)
        grad[:] = 2 * kappa
        h_diag[:] = 3 * kappa

        return func, self.mock_hess_x

    def mock_precond(self, residual, mu, precond_residual):
        """
        this function is a mock function for the preconditioner function
        """
        precond_residual[:] = mu * residual

    def mock_project(self, vector):
        """
        this function is a mock function for the projection function
        """
        vector[:] = 2 * vector

    def mock_modify_step(self, kappa):
        """
        this function is a mock function for the step modification function
        """
        kappa[:] = 2 * kappa

    def mock_conv_check(self):
        """
        this function is a mock function for the convergence check function
        """
        return True

    mock_approx_hess_x = mock_hess_x
    mock_stability_hess_x = mock_hess_x

    def mock_logger(self, message):
        """
        this function is a mock function for the logging function
        """
        if message != "test":
            self.test_logger = False
        return

    # replace original library with mock library
    @patch("pyopentrustregion.python_interface.lib.solver", lib.mock_solver)
    def test_solver_py_interface(self):
        """
        this function tests the solver python interface
        """
        # initialize settings object and assign reference values
        settings = SolverSettings()
        self.assign_ref_to_settings(settings)

        # initialize logging boolean
        self.test_logger = True

        # call solver python interface with optional arguments
        solver(self.mock_obj_func, self.mock_update_orbs, n_param, settings)

        # check if logger was called correctly
        if not self.test_logger:
            print(" test_solver_py_interface failed: Called logging function wrong.")

        self.assertTrue(
            c_bool.in_dll(lib, "test_solver_interface").value and self.test_logger,
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

        # initialize settings object and assign reference values
        settings = StabilitySettings()
        self.assign_ref_to_settings(settings)

        # allocate memory for descent direction
        kappa = np.empty(n_param, dtype=np.float64)

        # initialize logging boolean
        self.test_logger = True

        # call stability check python interface with optional arguments
        stable = stability_check(
            h_diag, self.mock_hess_x, n_param, settings, kappa=kappa
        )

        # check if logger was called correctly
        if not self.test_logger:
            print(
                " test_stability_check_py_interface failed: Called logging function "
                "wrong."
            )

        # check if returned variables are correct
        if stable:
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

        self.assertTrue(
            c_bool.in_dll(lib, "test_stability_check_interface").value
            and self.test_logger
            and not stable
            and not wrong_direction,
            "test_stability_check_py_interface failed",
        )
        print(" test_stability_check_py_interface PASSED")

    @patch.object(SolverSettings, "init_c_struct", lib.mock_init_solver_settings)
    def test_solver_settings(self):
        """
        this function ensure the SolverSettings object is properly initialized and
        synchronized with the underlying C struct
        """
        test_passed = True
        settings = SolverSettings()
        if not self.equal_settings_to_ref(settings):
            print(" test_solver_settings failed: Settings not initialized correctly.")
            test_passed = False

        dummy_error_code = 42

        def dummy_precond():
            return dummy_error_code

        settings.set_optional_callback(
            "precond", dummy_precond, lambda x: x, CFUNCTYPE(c_int)
        )

        c_ptr = getattr(settings.settings_c, "precond")
        c_interface = getattr(settings.settings_c, "precond_interface", None)

        if (
            c_ptr is None
            or (isinstance(c_ptr, c_void_p) and c_ptr.value is None)
            or not callable(c_interface)
            or c_interface() != dummy_error_code
        ):
            print(
                " test_solver_settings failed: Optional callbacks are not set "
                "correctly."
            )
            test_passed = False

        self.assertTrue(test_passed, "test_solver_settings failed")
        print(" test_solver_settings PASSED")

    @patch.object(StabilitySettings, "init_c_struct", lib.mock_init_stability_settings)
    def test_stability_settings(self):
        """
        this function ensure the StabilitySettings object is properly initialized and
        synchronized with the underlying C struct
        """
        test_passed = True
        settings = StabilitySettings()
        if not self.equal_settings_to_ref(settings):
            print(
                " test_stability_settings failed: Settings not initialized correctly."
            )
            test_passed = False

        self.assertTrue(test_passed, "test_stability_settings failed")
        print(" test_stability_settings PASSED")


@add_tests
class SystemTests(unittest.TestCase):
    """
    this class contains system tests for opentrustregion
    """

    tests = fortran_tests["system_tests"]

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
