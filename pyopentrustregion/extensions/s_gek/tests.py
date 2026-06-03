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
from pyopentrustregion.extensions.s_gek import (
    SGEKSettings,
    update_orbs_s_gek_factory,
    update_orbs_s_gek_deconstructor,
)

if NUMPY_AVAILABLE:
    import numpy as np


# define all tests in alphabetical order
fortran_tests = {
    "s_gek_c_interface_tests": [
        "assign_s_gek_c_f",
        "assign_s_gek_f_c",
        "hess_x_s_gek_c_wrapper",
        "init_s_gek_settings_c",
        "update_orbs_orig_s_gek_f_wrapper",
        "change_reference_f_wrapper",
        "update_orbs_s_gek_c_wrapper",
        "update_orbs_s_gek_deconstructor_c_wrapper",
        "update_orbs_s_gek_factory_c_wrapper",
    ],
}


@add_tests
class SGEKCInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the S-GEK C interface
    """

    tests = fortran_tests["s_gek_c_interface_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for S-GEK C interface...")
        return super().setUpClass()


@unittest.skipUnless(NUMPY_AVAILABLE, "NumPy not available.")
class SGEKPyInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the Python interface
    """

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for S-GEK Python interface...")

        # get fields
        fields = SGEKSettings.c_struct._fields_

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
        lib.get_reference_s_gek_values.argtypes = [POINTER(RefSettingsC)]
        lib.get_reference_s_gek_values.restype = None
        lib.get_reference_s_gek_values(byref(ref_settings))

        # extract Python values
        for name, _ in ref_fields:
            ref_value = getattr(ref_settings, name)
            if isinstance(ref_value, bytes):
                ref_value = ref_value.decode("utf-8")
            setattr(cls, name, ref_value)

        return super().setUpClass()

    assign_ref_to_settings = PyInterfaceTests.assign_ref_to_settings

    mock_update_orbs = PyInterfaceTests.mock_update_orbs
    mock_hess_x = PyInterfaceTests.mock_hess_x

    def mock_change_reference(self, new_ref, kappa_list, local_grad_list, grad_list):
        """
        this function is a mock function for the change reference function
        """
        self.change_reference_new_ref = np.allclose(
            new_ref, np.ones(n_param, dtype=np.float64)
        )

        kappa_list *= 2
        local_grad_list *= 3
        grad_list *= 4

        return

    mock_logger = PyInterfaceTests.mock_logger

    # replace original library with mock library
    @patch(
        "pyopentrustregion.python_interface.lib.update_orbs_s_gek_factory",
        lib.mock_update_orbs_s_gek_factory,
    )
    def test_update_orbs_s_gek_factory_py_interface(self):
        """
        this function tests the S-GEK orbital updating factory python interface
        """
        # initialize test flag
        test_passed = True

        # initialize settings object
        settings = SGEKSettings()
        self.assign_ref_to_settings(settings)

        # initialize logging boolean
        self.test_logger = True

        # call S-GEK orbital updating factory python interface
        update_orbs_s_gek = update_orbs_s_gek_factory(
            self.mock_update_orbs, self.mock_change_reference, n_param, settings
        )

        # check if logger was called correctly
        if not self.test_logger:
            print(
                " test_update_orbs_s_gek_factory_py_interface failed: Called logging "
                "function wrong."
            )
            test_passed = False

        # call returned S-GEK orbital updating function
        kappa = np.ones(n_param, dtype=np.float64)
        grad = np.empty(n_param, dtype=np.float64)
        h_diag = np.empty(n_param, dtype=np.float64)
        try:
            func, hess_x_funptr = update_orbs_s_gek(kappa, grad, h_diag)
        except RuntimeError:
            print(
                " test_update_orbs_s_gek_factory_py_interface failed: Returned S-GEK "
                "orbital updating raises error."
            )
            test_passed = False

        # check results
        if self.change_reference_new_ref is not True:
            print(
                " test_update_orbs_s_gek_factory_py_interface failed: New "
                "reference parameters inside given change reference function wrong."
            )
            test_passed = False
        if func != 3.0:
            print(
                " test_update_orbs_s_gek_factory_py_interface failed: Returned"
                " function value of returned S-GEK orbital updating function wrong."
            )
            test_passed = False
        if not np.allclose(grad, np.full(n_param, 2.0, dtype=np.float64)):
            print(
                " test_update_orbs_s_gek_factory_py_interface failed: Returned "
                "gradient of returned S-GEK orbital updating function wrong."
            )
            test_passed = False
        if not np.allclose(h_diag, np.full(n_param, 3.0, dtype=np.float64)):
            print(
                " test_update_orbs_s_gek_factory_py_interface failed: Returned "
                "S-GEK Hessian diagonal of returned S-GEK orbital updating function "
                "wrong."
            )
            test_passed = False

        # call returned hess_x function
        x = np.ones(n_param, dtype=np.float64)
        hess_x = np.empty(n_param, dtype=np.float64)

        try:
            hess_x_funptr(x, hess_x)
        except RuntimeError:
            print(
                " test_update_orbs_s_gek_factory_py_interface failed: Returned "
                "S-GEK Hessian linear transformation function of returned S-GEK "
                "orbital updating function raises error."
            )
            test_passed = False

        # check results
        if not np.allclose(hess_x, np.full(n_param, 4.0, dtype=np.float64)):
            print(
                " test_update_orbs_s_gek_factory_py_interface failed: Returned Hessian "
                "linear transformation of returned S-GEK orbital updating function "
                "wrong."
            )
            test_passed = False

        self.assertTrue(
            c_bool.in_dll(lib, "test_update_orbs_s_gek_factory_interface").value
            and test_passed,
            "test_update_orbs_s_gek_factory_py_interface failed",
        )
        print(" test_update_orbs_s_gek_factory_py_interface PASSED")

    # replace original library with mock library
    @patch(
        "pyopentrustregion.python_interface.lib.update_orbs_s_gek_deconstructor",
        lib.mock_update_orbs_s_gek_deconstructor,
    )
    def test_update_orbs_deconstructor_py_interface(self):
        """
        this function tests the S-GEK orbital updating deconstructor python
        interface
        """
        # initialize test flag
        test_passed = True

        # call S-GEK orbital updating deconstructor python interface
        update_orbs_s_gek_deconstructor()

        # check if deconstructor was called correctly
        if not c_bool.in_dll(
            lib, "test_update_orbs_s_gek_deconstructor_interface"
        ).value:
            print(
                " test_update_orbs_s_gek_deconstructor_py_interface failed: Called "
                "logging deconstructor wrong."
            )
            test_passed = False

        self.assertTrue(
            test_passed, "test_update_orbs_s_gek_deconstructor_py_interface failed"
        )
        print(" test_update_orbs_s_gek_deconstructor_py_interface PASSED")

    @patch.object(SGEKSettings, "init_c_struct", lib.mock_init_s_gek_settings)
    def test_s_gek_settings(self):
        """
        this function ensure the SGEKSettings object is properly initialized and
        synchronized with the underlying C struct
        """
        test_passed = True
        settings = SGEKSettings()
        if not self.equal_settings_to_ref(settings):
            print(" test_s_gek_settings failed: Settings not initialized correctly.")

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
                " test_s_gek_settings failed: Optional callbacks are not set correctly."
            )
            test_passed = False

        self.assertTrue(test_passed, "test_s_gek_settings failed")
        print(" test_s_gek_settings PASSED")
