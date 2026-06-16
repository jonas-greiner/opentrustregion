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
from pyopentrustregion.extensions.oao import OAOSettings, oao_factory, oao_deconstructor

if NUMPY_AVAILABLE:
    import numpy as np


# define all tests in alphabetical order
fortran_tests = {
    "oao_c_interface_tests": [
        "assign_oao_c_f",
        "assign_oao_f_c",
        "get_energy_f_wrapper",
        "get_response_f_wrapper",
        "init_oao_settings_c",
        "oao_deconstructor_c_wrapper",
        "oao_factory_c_wrapper",
        "obj_func_oao_c_wrapper",
        "project_oao_c_wrapper",
        "update_dm_f_wrapper",
    ],
}

# number of AOs
n_ao = c_int.in_dll(lib, "test_n_ao").value


@add_tests
class OAOCInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the OAO C interface
    """

    tests = fortran_tests["oao_c_interface_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for OAO C interface...")
        return super().setUpClass()


@unittest.skipUnless(NUMPY_AVAILABLE, "NumPy not available.")
class OAOPyInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the Python interface
    """

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for OAO Python interface...")

        # get fields
        fields = OAOSettings.c_struct._fields_

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
        lib.get_reference_oao_values.argtypes = [POINTER(RefSettingsC)]
        lib.get_reference_oao_values.restype = None
        lib.get_reference_oao_values(byref(ref_settings))

        # extract Python values
        for name, _ in ref_fields:
            ref_value = getattr(ref_settings, name)
            if isinstance(ref_value, bytes):
                ref_value = ref_value.decode("utf-8")
            setattr(cls, name, ref_value)

        return super().setUpClass()

    assign_ref_to_settings = PyInterfaceTests.assign_ref_to_settings
    equal_settings_to_ref = PyInterfaceTests.equal_settings_to_ref

    mock_logger = PyInterfaceTests.mock_logger

    def mock_get_energy(self, dm_ao):
        """
        this function is a mock function for the energy function
        """
        return np.sum(dm_ao)

    def mock_update_dm(self, dm_ao, fock):
        """
        this function is a mock function for the density matrix updating function
        """
        fock[:] = 2 * dm_ao

        return np.sum(dm_ao), self.mock_get_response

    def mock_get_response(self, dm_ao, response):
        """
        this function is a mock function for the response function
        """
        response[:] = 2 * dm_ao

        return

    # replace original library with mock library
    @patch("pyopentrustregion.python_interface.lib.oao_factory", lib.mock_oao_factory)
    def test_oao_factory_py_interface(self):
        """
        this function tests the OAO factory python interface (only tests if dm_ao,
        mock_get_energy and mock_update_dm_jk are passed correctly for the open-shell
        case since everything else is the same in the closed-shell case)
        """
        ao_overlap = np.full(2 * (n_ao,), 2.0, dtype=np.float64)

        # initialize test flag
        test_passed = True

        # initialize settings object
        settings = OAOSettings()
        self.assign_ref_to_settings(settings)

        # initialize logging boolean
        self.test_logger = True

        # number of particles
        n_particle = 1

        # initialize density matrix
        dm_ao = np.full(2 * (n_ao,), 1.0, dtype=np.float64)

        # call OAO factory python interface
        obj_func_oao, update_orbs_oao, project_oao = oao_factory(
            dm_ao,
            ao_overlap,
            n_particle,
            n_ao,
            self.mock_get_energy,
            self.mock_update_dm,
            settings,
        )

        # check if logger was called correctly
        if not self.test_logger:
            print(
                " test_oao_factory_py_interface failed: Called logging function wrong."
            )
            test_passed = False

        # call returned OAO objective function
        kappa = np.ones(n_param, dtype=np.float64)
        try:
            func = obj_func_oao(kappa)
        except RuntimeError:
            print(
                " test_oao_factory_py_interface failed: Returned OAO objective "
                "function raises error."
            )
            test_passed = False

        # check results
        if func != 3.0:
            print(
                " test_oao_factory_py_interface failed: Returned function value of "
                "returned OAO objective function wrong."
            )
            test_passed = False

        # call returned OAO orbital updating function
        grad = np.empty(n_param, dtype=np.float64)
        h_diag = np.empty(n_param, dtype=np.float64)
        try:
            func, hess_x_funptr = update_orbs_oao(kappa, grad, h_diag)
        except RuntimeError:
            print(
                " test_oao_factory_py_interface failed: Returned OAO orbital updating "
                "function raises error."
            )
            test_passed = False

        # check if density matrix was updated
        if not np.allclose(dm_ao, np.full(2 * (n_ao,), 2.0, dtype=np.float64)):
            print(
                " test_oao_factory_py_interface failed: Density matrix not updated "
                "correctly."
            )
            test_passed = False

        # check results
        if func != 3.0:
            print(
                " test_oao_factory_py_interface failed: Returned function value of "
                "returned OAO orbital updating function wrong."
            )
            test_passed = False
        if not np.allclose(grad, np.full(n_param, 2.0, dtype=np.float64)):
            print(
                " test_oao_factory_py_interface failed: Returned gradient of returned "
                "OAO orbital updating function wrong."
            )
            test_passed = False
        if not np.allclose(h_diag, np.full(n_param, 3.0, dtype=np.float64)):
            print(
                " test_oao_factory_py_interface failed: Returned OAO Hessian diagonal "
                "of returned OAO orbital updating function wrong."
            )
            test_passed = False

        # call returned hess_x function
        x = np.ones(n_param, dtype=np.float64)
        hess_x = np.empty(n_param, dtype=np.float64)

        try:
            hess_x_funptr(x, hess_x)
        except RuntimeError:
            print(
                " test_oao_factory_py_interface failed: Returned OAO Hessian linear "
                "transformation function of returned OAO orbital updating function "
                "raises error."
            )
            test_passed = False

        # check results
        if not np.allclose(hess_x, np.full(n_param, 4.0, dtype=np.float64)):
            print(
                " test_oao_factory_py_interface failed: Returned OAO Hessian linear "
                "transformation of returned OAO orbital updating function wrong."
            )
            test_passed = False

        # call returned OAO projection function
        vector = np.full(n_param, 1.0, dtype=np.float64)
        try:
            project_oao(vector)
        except RuntimeError:
            print(
                " test_oao_factory_py_interface failed: Returned OAO projection "
                "function raises error."
            )
            test_passed = False

        # check results
        if not np.allclose(vector, np.full(n_param, 2.0, dtype=np.float64)):
            print(
                " test_oao_factory_py_interface failed: Returned projected vector "
                "of returned OAO projection function wrong."
            )
            test_passed = False

        # number of particles
        n_particle = 2

        # initialize density matrix
        dm_ao = np.full((n_particle, n_ao, n_ao), 1.0, dtype=np.float64)

        # call OAO factory python interface
        oao_factory(
            dm_ao,
            ao_overlap,
            n_particle,
            n_ao,
            self.mock_get_energy,
            self.mock_update_dm,
            settings,
        )

        self.assertTrue(
            c_bool.in_dll(lib, "test_oao_factory_interface").value and test_passed,
            "test_oao_factory_py_interface failed",
        )
        print(" test_oao_factory_py_interface PASSED")

    # replace original library with mock library
    @patch(
        "pyopentrustregion.python_interface.lib.oao_deconstructor",
        lib.mock_oao_deconstructor,
    )
    def test_oao_deconstructor_py_interface(self):
        """
        this function tests the OAO deconstructor python interface
        """
        # initialize test flag
        test_passed = True

        # call OAO deconstructor python interface
        oao_deconstructor()

        # check if deconstructor was called correctly
        if not c_bool.in_dll(lib, "test_oao_deconstructor_interface").value:
            print(
                " test_oao_deconstructor_py_interface failed: Deconstructor called "
                "wrong."
            )
            test_passed = False

        self.assertTrue(test_passed, "test_oao_deconstructor_py_interface failed")
        print(" test_oao_deconstructor_py_interface PASSED")

    @patch.object(OAOSettings, "init_c_struct", lib.mock_init_oao_settings)
    def test_oao_settings(self):
        """
        this function ensure the OAOSettings object is properly initialized and
        synchronized with the underlying C struct
        """
        test_passed = True
        settings = OAOSettings()
        if not self.equal_settings_to_ref(settings):
            print(" test_oao_settings failed: Settings not initialized correctly.")
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
                " test_oao_settings failed: Optional callbacks are not set correctly."
            )
            test_passed = False

        self.assertTrue(test_passed, "test_oao_settings failed")
        print(" test_oao_settings PASSED")
