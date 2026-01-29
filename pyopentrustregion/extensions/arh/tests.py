# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import unittest
from ctypes import (
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

from pyopentrustregion.tests import lib, NUMPY_AVAILABLE, add_tests, print_separator
from pyopentrustregion.python_interface import c_real, c_int

from pyopentrustregion.extensions.arh import ARHSettings, arh_factory, arh_deconstructor

if NUMPY_AVAILABLE:
    import numpy as np


# define all tests in alphabetical order
fortran_tests = {
    "arh_c_interface_tests": [
        "arh_deconstructor_c_wrapper",
        "arh_factory_c_wrapper",
        "assign_arh_c_f",
        "assign_arh_f_c",
        "get_energy_f_wrapper",
        "get_fock_f_wrapper",
        "get_fock_jk_f_wrapper",
        "hess_x_arh_c_wrapper",
        "init_arh_settings_c",
        "update_orbs_arh_c_wrapper",
    ],
}

# number of AOs
n_ao = 3


@add_tests
class ARHCInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the ARH C interface
    """

    tests = fortran_tests["arh_c_interface_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for ARH C interface...")
        return super().setUpClass()


@unittest.skipUnless(NUMPY_AVAILABLE, "NumPy not available.")
class ARHPyInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the Python interface
    """

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for ARH Python interface...")

        # get fields
        fields = ARHSettings.c_struct._fields_

        # get fields by order
        ref_fields = []
        # loop by type order
        for curr_type in (c_bool, c_real, c_int):
            # fields of this type
            for name, t in fields:
                if t == curr_type and name != "initialized":
                    ref_fields.append((name + "_ref", t))

        # handle fixed-size c_char arrays (strings)
        for name, t in fields:
            if issubclass(t, Array) and getattr(t, "_type_", None) is c_char:
                ref_fields.append((name + "_ref", t))

        # create class to read reference values
        class RefSettingsC(Structure):
            _fields_ = ref_fields

        # create instance
        ref_settings = RefSettingsC()

        # call Fortran to fill values
        lib.get_reference_arh_values.argtypes = [POINTER(RefSettingsC)]
        lib.get_reference_arh_values.restype = None
        lib.get_reference_arh_values(byref(ref_settings))

        # extract Python values
        for name, _ in ref_fields:
            ref_value = getattr(ref_settings, name)
            if isinstance(ref_value, bytes):
                ref_value = ref_value.decode("utf-8")
            setattr(cls, name, ref_value)

        return super().setUpClass()

    # replace original library with mock library
    @patch("pyopentrustregion.python_interface.lib.arh_factory", lib.mock_arh_factory)
    def test_arh_factory_py_interface(self):
        """
        this function tests the ARH factory python interface (only tests if dm_ao,
        mock_get_energy and mock_get_fock_jk are passed correctly for the open-shell
        case since everything else is the same in the closed-shell case)
        """
        n_param = n_ao * (n_ao - 1) // 2
        ao_overlap = np.full(2 * (n_ao,), 2.0, dtype=np.float64)

        def mock_get_energy(dm_ao):
            """
            this function is a mock function for the energy function
            """
            return np.sum(dm_ao)

        def mock_get_fock(dm_ao, fock):
            """
            this function is a mock function for the Fock matrix function
            """
            fock[:] = 2 * dm_ao

            return np.sum(dm_ao)

        def mock_get_fock_jk(dm_ao, fock, coulomb, exchange):
            """
            this function is a mock function for the Fock matrix function
            """
            fock[:] = 2 * dm_ao
            coulomb[:] = 3 * dm_ao
            exchange[:] = 4 * dm_ao

            return np.sum(dm_ao)

        def mock_logger(message):
            """
            this function is a mock function for the logging function
            """
            nonlocal test_logger
            if message == "test":
                test_logger = True
            return

        # initialize test flag
        test_passed = True

        # initialize settings object
        settings = ARHSettings()
        settings.logger = mock_logger
        for field_info in settings.c_struct._fields_:
            field_name, field_type = field_info[:2]
            if field_type == c_void_p or field_name == "initialized":
                continue
            setattr(settings, field_name, getattr(self, field_name + "_ref"))

        # initialize logging boolean
        test_logger = False

        # number of particles
        n_particle = 1

        # initialize density matrix
        dm_ao = np.full(2 * (n_ao,), 1.0, dtype=np.float64)

        # call ARH factory python interface
        obj_func_arh, update_orbs_arh, precond_arh = arh_factory(
            dm_ao,
            ao_overlap,
            n_particle,
            n_ao,
            mock_get_energy,
            mock_get_fock,
            settings,
        )

        # check if logger was called correctly
        if not test_logger:
            print(
                " test_arh_factory_py_interface failed: Called logging function wrong."
            )
            test_passed = False

        # call returned ARH objective function
        kappa = np.ones(n_param, dtype=np.float64)
        try:
            func = obj_func_arh(kappa)
        except RuntimeError:
            print(
                " test_arh_factory_py_interface failed: Returned ARH objective "
                "function raises error."
            )
            test_passed = False

        # check results
        if func != 3.0:
            print(
                " test_arh_factory_py_interface failed: Returned function value of "
                "returned ARH objective function wrong."
            )
            test_passed = False

        # call returned ARH orbital updating function
        grad = np.empty(n_param, dtype=np.float64)
        h_diag = np.empty(n_param, dtype=np.float64)
        try:
            func, hess_x_funptr = update_orbs_arh(kappa, grad, h_diag)
        except RuntimeError:
            print(
                " test_arh_factory_py_interface failed: Returned ARH orbital updating "
                "function raises error."
            )
            test_passed = False

        # check results
        if func != 3.0:
            print(
                " test_arh_factory_py_interface failed: Returned function value of "
                "returned ARH orbital updating function wrong."
            )
            test_passed = False
        if not np.allclose(grad, np.full(n_param, 2.0, dtype=np.float64)):
            print(
                " test_arh_factory_py_interface failed: Returned gradient of returned "
                "ARH orbital updating function wrong."
            )
            test_passed = False
        if not np.allclose(h_diag, np.full(n_param, 3.0, dtype=np.float64)):
            print(
                " test_arh_factory_py_interface failed: Returned ARH Hessian diagonal "
                "of returned ARH orbital updating function wrong."
            )
            test_passed = False

        # call returned hess_x function
        x = np.ones(n_param, dtype=np.float64)
        hess_x = np.empty(n_param, dtype=np.float64)

        try:
            hess_x_funptr(x, hess_x)
        except RuntimeError:
            print(
                " test_arh_factory_py_interface failed: Returned ARH Hessian linear "
                "transformation function of returned ARH orbital updating function "
                "raises error."
            )
            test_passed = False

        # check results
        if not np.allclose(hess_x, np.full(n_param, 4.0, dtype=np.float64)):
            print(
                " test_arh_factory_py_interface failed: Returned Hessian linear "
                "transformation of returned ARH orbital updating function wrong."
            )
            test_passed = False

        # call returned ARH preconditioning function
        residual = np.full(n_param, 1.0, dtype=np.float64)
        precond_residual = np.empty(n_param, dtype=np.float64)
        try:
            precond_arh(residual, 2.0, precond_residual)
        except RuntimeError:
            print(
                " test_arh_factory_py_interface failed: Returned ARH preconditioning "
                "function raises error."
            )
            test_passed = False

        # check results
        if not np.allclose(precond_residual, np.full(n_param, 2.0, dtype=np.float64)):
            print(
                " test_arh_factory_py_interface failed: Returned preconditioned "
                "residual of returned ARH preconditioning function wrong."
            )
            test_passed = False

        # number of particles
        n_particle = 2

        # initialize density matrix
        dm_ao = np.full((n_particle, n_ao, n_ao), 1.0, dtype=np.float64)

        # call ARH factory python interface
        arh_factory(
            dm_ao,
            ao_overlap,
            n_particle,
            n_ao,
            mock_get_energy,
            mock_get_fock_jk,
            settings,
        )

        self.assertTrue(
            c_bool.in_dll(lib, "test_arh_factory_interface").value and test_passed,
            "test_arh_factory_py_interface failed",
        )
        print(" test_arh_factory_py_interface PASSED")

    def test_arh_deconstructor_py_interface(self):
        """
        this function tests the ARH deconstructor python interface
        """

        @patch(
            "pyopentrustregion.python_interface.lib.arh_deconstructor",
            lib.mock_arh_deconstructor_2d,
        )
        def test_arh_deconstructor_py_interface_2d():
            """
            this function tests the ARH deconstructor python interface for 2D density
            matrices
            """
            # initialize density matrix
            dm_ao = np.empty(2 * (n_ao,), dtype=np.float64)

            # call ARH deconstructor python interface
            arh_deconstructor(dm_ao)

            # check results
            if not np.allclose(dm_ao, np.full(2 * (n_ao,), 1.0, dtype=np.float64)):
                print(
                    " test_arh_deconstructor_py_interface failed: Returned AO density "
                    "matrix wrong for 2D density matrices."
                )
                return False

            return True

        @patch(
            "pyopentrustregion.python_interface.lib.arh_deconstructor",
            lib.mock_arh_deconstructor_3d,
        )
        def test_arh_deconstructor_py_interface_3d():
            """
            this function tests the ARH deconstructor python interface for 3D density
            matrices
            """
            # number of particles
            n_particle = 2

            # initialize density matrix
            dm_ao = np.empty((n_particle, n_ao, n_ao), dtype=np.float64)

            # call ARH deconstructor python interface
            arh_deconstructor(dm_ao)

            # check results
            if not np.allclose(
                dm_ao, np.full((n_particle, n_ao, n_ao), 1.0, dtype=np.float64)
            ):
                print(
                    " test_arh_deconstructor_py_interface failed: Returned AO density "
                    "matrix wrong for 3D density matrices."
                )
                return False

            return True

        # run test for 2D density matrices
        test_passed_2d = test_arh_deconstructor_py_interface_2d()

        # run test for 3D density matrices
        test_passed_3d = test_arh_deconstructor_py_interface_3d()

        self.assertTrue(
            test_passed_2d and test_passed_3d,
            "test_arh_deconstructor_py_interface failed",
        )
        print(" test_arh_deconstructor_py_interface PASSED")

    @patch.object(ARHSettings, "init_c_struct", lib.mock_init_arh_settings)
    def test_arh_settings(self):
        """
        this function ensure the ARHSettings object is properly initialized and
        synchronized with the underlying C struct
        """
        settings = ARHSettings()
        test_passed = True
        for field_info in settings.c_struct._fields_:
            field_name, field_type = field_info[:2]
            if field_type == c_void_p:
                if (
                    getattr(settings, field_name) is not None
                    or getattr(settings.settings_c, field_name) is not None
                ):
                    print(
                        " test_arh_settings failed: Optional function pointer "
                        f"{field_name} not initialized correctly."
                    )
                    test_passed = False
            elif field_name == "initialized":
                if not getattr(settings, field_name):
                    print(
                        " test_arh_settings failed: Field initialized not initialized "
                        "correctly."
                    )
                    test_passed = False
            else:
                ref_value = getattr(self, field_name + "_ref")
                if field_type == c_real:
                    match = np.isclose(getattr(settings, field_name), ref_value)
                else:
                    match = getattr(settings, field_name) == ref_value
                if not match:
                    print(field_name, getattr(settings, field_name), ref_value)
                    print(
                        f" test_arh_settings failed: Field {field_name} not "
                        "initialized correctly."
                    )
                    test_passed = False

        dummy_error_code = 42

        def dummy_logger():
            return dummy_error_code

        settings.set_optional_callback("logger", CFUNCTYPE(c_int)(dummy_logger))

        c_ptr = getattr(settings.settings_c, "logger")
        c_interface = getattr(settings.settings_c, "logger_interface", None)

        if (
            c_ptr is None
            or (isinstance(c_ptr, c_void_p) and c_ptr.value is None)
            or not callable(c_interface)
            or c_interface() != dummy_error_code
        ):
            print(
                " test_arh_settings failed: Optional callbacks are not set correctly."
            )
            test_passed = False

        self.assertTrue(test_passed, "test_arh_settings failed")
        print(" test_arh_settings PASSED")
