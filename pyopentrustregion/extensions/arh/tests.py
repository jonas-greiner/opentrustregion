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

from pyopentrustregion.tests import (
    lib,
    n_param,
    NUMPY_AVAILABLE,
    add_tests,
    print_separator,
    PyInterfaceTests,
)
from pyopentrustregion.python_interface import c_real, c_int
from pyopentrustregion.extensions.arh import ARHSettings, arh_factory, arh_deconstructor
from pyopentrustregion.extensions.oao.tests import n_ao, OAOPyInterfaceTests

if NUMPY_AVAILABLE:
    import numpy as np


# define all tests in alphabetical order
fortran_tests = {
    "arh_tests": [
        "arh_deconstructor",
        "arh_factory_cs",
        "arh_factory_os",
        "arh_sanity_check",
        "cross_symmetrize_weighted",
        "get_arh_metric",
        "get_ms_a_inv_cs",
        "get_ms_a_inv_os_linear",
        "get_ms_a_inv_os_nonlinear",
        "get_response_contribution_cs",
        "get_response_contribution_ms_sr1_nonlinear",
        "get_response_contribution_os_separated",
        "hess_x_arh",
        "init_arh_settings",
        "multiply_with_inverse_metric",
        "noise_threshold",
        "prepend",
        "regularized_eigval_inv",
        "symmetrize_exact",
        "symmetrize_weighted",
        "truncated_eigval_inv",
        "update_orbs_arh_cs",
        "update_orbs_arh_os",
    ],
    "arh_c_interface_tests": [
        "arh_deconstructor_c_wrapper",
        "arh_factory_c_wrapper",
        "assign_arh_c_f",
        "assign_arh_f_c",
        "hess_x_arh_c_wrapper",
        "init_arh_settings_c",
        "update_dm_cs_f_wrapper",
        "update_dm_os_f_wrapper",
        "update_orbs_arh_c_wrapper",
    ],
}


@add_tests
class ARHTests(unittest.TestCase):
    """
    this class contains unit tests for ARH
    """

    tests = fortran_tests["arh_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for ARH...")
        return super().setUpClass()


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

    assign_ref_to_settings = PyInterfaceTests.assign_ref_to_settings
    equal_settings_to_ref = PyInterfaceTests.equal_settings_to_ref

    mock_get_energy = OAOPyInterfaceTests.mock_get_energy

    def mock_update_dm_cs(self, dm_ao, fock, v_nonlinear):
        """
        this function is a mock function for the density matrix updating function with
        a separate non-linear potential contribution for the closed-shell case
        """
        fock[:] = 2 * dm_ao
        v_nonlinear[:] = 3 * dm_ao

        return np.sum(dm_ao)

    def mock_update_dm_os(self, dm_ao, fock, v_same_spin, v_opposite_spin, v_nonlinear):
        """
        this function is a mock function for the density matrix updating function with
        separate same- and opposite spin and non-linear contributions for the
        open-shell case
        """
        fock[:] = 2 * dm_ao
        v_same_spin[:] = 3 * dm_ao
        v_opposite_spin[:] = 4 * dm_ao
        v_nonlinear[:] = 5 * dm_ao

        return np.sum(dm_ao)

    mock_logger = PyInterfaceTests.mock_logger

    # replace original library with mock library
    @patch("pyopentrustregion.python_interface.lib.arh_factory", lib.mock_arh_factory)
    def test_arh_factory_py_interface(self):
        """
        this function tests the ARH factory python interface (only tests if dm_ao,
        mock_get_energy and mock_update_dm_os are passed correctly for the open-shell
        case since everything else is the same in the closed-shell case)
        """
        ao_overlap = np.full(2 * (n_ao,), 2.0, dtype=np.float64)

        # initialize test flag
        test_passed = True

        # initialize settings object
        settings = ARHSettings()
        self.assign_ref_to_settings(settings)

        # initialize logging boolean
        self.test_logger = True

        # number of particles
        n_particle = 1

        # initialize density matrix
        dm_ao = np.full(2 * (n_ao,), 1.0, dtype=np.float64)

        # call ARH factory python interface
        obj_func_arh, update_orbs_arh, precond_arh, precond_pd_arh, project_arh = (
            arh_factory(
                dm_ao,
                ao_overlap,
                n_particle,
                n_ao,
                self.mock_get_energy,
                self.mock_update_dm_cs,
                settings,
            )
        )

        # check if logger was called correctly
        if not self.test_logger:
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

        # check if density matrix was updated
        if not np.allclose(dm_ao, np.full(2 * (n_ao,), 2.0, dtype=np.float64)):
            print(
                " test_arh_factory_py_interface failed: Density matrix not updated "
                "correctly."
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

        # call returned ARH projection function
        vector = np.full(n_param, 1.0, dtype=np.float64)
        try:
            project_arh(vector)
        except RuntimeError:
            print(
                " test_arh_factory_py_interface failed: Returned ARH projection "
                "function raises error."
            )
            test_passed = False

        # check results
        if not np.allclose(vector, np.full(n_param, 2.0, dtype=np.float64)):
            print(
                " test_arh_factory_py_interface failed: Returned projected vector "
                "of returned ARH projection function wrong."
            )
            test_passed = False

        # call returned ARH level-shifted preconditioner function
        residual = np.full(n_param, 1.0, dtype=np.float64)
        precond_residual = np.empty(n_param, dtype=np.float64)
        mu = 5.0
        try:
            precond_arh(residual, mu, precond_residual)
        except RuntimeError:
            print(
                " test_arh_factory_py_interface failed: Returned ARH level-shifted "
                "preconditioner function raises error."
            )
            test_passed = False

        # check results
        if not np.allclose(precond_residual, np.full(n_param, mu, dtype=np.float64)):
            print(
                " test_arh_factory_py_interface failed: Returned preconditioned "
                "residual of returned ARH level-shifted preconditioner function "
                "wrong."
            )
            test_passed = False

        # call returned ARH positive-definite preconditioner function
        precond_pd_residual = np.empty(n_param, dtype=np.float64)
        try:
            precond_pd_arh(residual, precond_pd_residual)
        except RuntimeError:
            print(
                " test_arh_factory_py_interface failed: Returned ARH "
                "positive-definite preconditioner function raises error."
            )
            test_passed = False

        # check results
        if not np.allclose(
            precond_pd_residual, np.full(n_param, 3.0, dtype=np.float64)
        ):
            print(
                " test_arh_factory_py_interface failed: Returned preconditioned "
                "residual of returned ARH positive-definite preconditioner function "
                "wrong."
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
            self.mock_get_energy,
            self.mock_update_dm_os,
            settings,
        )

        self.assertTrue(
            c_bool.in_dll(lib, "test_arh_factory_interface").value and test_passed,
            "test_arh_factory_py_interface failed",
        )
        print(" test_arh_factory_py_interface PASSED")

    # replace original library with mock library
    @patch(
        "pyopentrustregion.python_interface.lib.arh_deconstructor",
        lib.mock_arh_deconstructor,
    )
    def test_arh_deconstructor_py_interface(self):
        """
        this function tests the ARH deconstructor python interface
        """
        # initialize test flag
        test_passed = True

        # call ARH deconstructor python interface
        arh_deconstructor()

        # check if deconstructor was called correctly
        if not c_bool.in_dll(lib, "test_arh_deconstructor_interface").value:
            print(
                " test_arh_deconstructor_py_interface failed: Deconstructor called "
                "wrong."
            )
            test_passed = False

        self.assertTrue(test_passed, "test_arh_deconstructor_py_interface failed")
        print(" test_arh_deconstructor_py_interface PASSED")

    @patch.object(ARHSettings, "init_c_struct", lib.mock_init_arh_settings)
    def test_arh_settings(self):
        """
        this function ensure the ARHSettings object is properly initialized and
        synchronized with the underlying C struct
        """
        test_passed = True
        settings = ARHSettings()
        if not self.equal_settings_to_ref(settings):
            print(" test_arh_settings failed: Settings not initialized correctly.")
            test_passed = False

        dummy_error_code = 42

        def dummy_logger():
            return dummy_error_code

        settings.set_optional_callback(
            "logger", dummy_logger, lambda x: x, CFUNCTYPE(c_int)
        )

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
