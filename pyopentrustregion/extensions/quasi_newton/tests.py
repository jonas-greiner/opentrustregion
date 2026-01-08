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
from pyopentrustregion.extensions.quasi_newton import (
    QNSettings,
    update_orbs_qn_factory,
    update_orbs_qn_deconstructor,
)

if NUMPY_AVAILABLE:
    import numpy as np


# define all tests in alphabetical order
fortran_tests = {
    "qn_c_interface_tests": [
        "assign_qn_c_f",
        "assign_qn_f_c",
        "hess_x_qn_c_wrapper",
        "init_qn_settings_c",
        "update_orbs_orig_qn_f_wrapper",
        "change_reference_qn_f_wrapper",
        "update_orbs_qn_c_wrapper",
        "update_orbs_qn_deconstructor_c_wrapper",
        "update_orbs_qn_factory_c_wrapper",
    ],
}


@add_tests
class QNCInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the quasi-Newton C interface
    """

    tests = fortran_tests["qn_c_interface_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for quasi-Newton C interface...")
        return super().setUpClass()


@unittest.skipUnless(NUMPY_AVAILABLE, "NumPy not available.")
class QNPyInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the Python interface
    """

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for quasi-Newton Python interface...")

        # get fields
        fields = QNSettings.c_struct._fields_

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
        lib.get_reference_qn_values.argtypes = [POINTER(RefSettingsC)]
        lib.get_reference_qn_values.restype = None
        lib.get_reference_qn_values(byref(ref_settings))

        # extract Python values
        for name, _ in ref_fields:
            ref_value = getattr(ref_settings, name)
            if isinstance(ref_value, bytes):
                ref_value = ref_value.decode("utf-8")
            setattr(cls, name, ref_value)

        return super().setUpClass()

    # replace original library with mock library
    @patch(
        "pyopentrustregion.python_interface.lib.update_orbs_qn_factory",
        lib.mock_update_orbs_qn_factory,
    )
    def test_update_orbs_qn_factory_py_interface(self):
        """
        this function tests the quasi-Newton orbital updating factory python interface
        """
        n_param = 3

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

        def mock_change_reference(new_ref, kappa_list, local_grad_list, grad_list):
            """
            this function is a mock function for the change reference function
            """
            nonlocal test_passed

            if not np.allclose(new_ref, np.ones(n_param, dtype=np.float64)):
                print(
                    " test_update_orbs_s_gek_factory_py_interface failed: New "
                    "reference parameters inside given change reference function wrong."
                )
                test_passed = False

            kappa_list *= 2
            local_grad_list *= 3
            grad_list *= 4

            return

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
        settings = QNSettings()
        settings.logger = mock_logger
        for field_info in settings.c_struct._fields_:
            field_name, field_type = field_info[:2]
            if field_type == c_void_p or field_name == "initialized":
                continue
            setattr(settings, field_name, getattr(self, field_name + "_ref"))

        # initialize logging boolean
        test_logger = False

        # call quasi-Newton orbital updating factory python interface
        update_orbs_qn = update_orbs_qn_factory(
            mock_update_orbs, mock_change_reference, n_param, settings
        )

        # check if logger was called correctly
        if not test_logger:
            print(
                " test_update_orbs_qn_factory_py_interface failed: Called logging "
                "function wrong."
            )
            test_passed = False

        # call returned quasi-Newton orbital updating function
        kappa = np.ones(n_param, dtype=np.float64)
        grad = np.empty(n_param, dtype=np.float64)
        h_diag = np.empty(n_param, dtype=np.float64)
        try:
            func, hess_x_funptr = update_orbs_qn(kappa, grad, h_diag)
        except RuntimeError:
            print(
                " test_update_orbs_qn_factory_py_interface failed: Returned "
                "quasi-Newton orbital updating raises error."
            )
            test_passed = False

        # check results
        if func != 3.0:
            print(
                " test_update_orbs_qn_factory_py_interface failed: Returned function "
                "value of returned quasi-Newton orbital updating function "
                "wrong."
            )
            test_passed = False
        if not np.allclose(grad, np.full(n_param, 2.0, dtype=np.float64)):
            print(
                " test_update_orbs_qn_factory_py_interface failed: Returned gradient "
                "of returned quasi-Newton orbital updating function wrong."
            )
            test_passed = False
        if not np.allclose(h_diag, np.full(n_param, 3.0, dtype=np.float64)):
            print(
                " test_update_orbs_qn_factory_py_interface failed: Returned "
                "quasi-Newton Hessian diagonal of returned quasi-Newton orbital "
                "updating function wrong."
            )
            test_passed = False

        # call returned hess_x function
        x = np.ones(n_param, dtype=np.float64)
        hess_x = np.empty(n_param, dtype=np.float64)

        try:
            hess_x_funptr(x, hess_x)
        except RuntimeError:
            print(
                " test_update_orbs_qn_factory_py_interface failed: Returned "
                "quasi-Newton Hessian linear transformation function of returned "
                "quasi-Newton orbital updating function raises error."
            )
            test_passed = False

        # check results
        if not np.allclose(hess_x, np.full(n_param, 4.0, dtype=np.float64)):
            print(
                " test_update_orbs_qn_factory_py_interface failed: Returned Hessian "
                "linear transformation of returned quasi-Newton orbital updating "
                "function wrong."
            )
            test_passed = False

        self.assertTrue(
            c_bool.in_dll(lib, "test_update_orbs_qn_factory_interface").value
            and test_passed,
            "test_update_orbs_qn_factory_py_interface failed",
        )
        print(" test_update_orbs_qn_factory_py_interface PASSED")

    # replace original library with mock library
    @patch(
        "pyopentrustregion.python_interface.lib.update_orbs_qn_deconstructor",
        lib.mock_update_orbs_qn_deconstructor,
    )
    def test_update_orbs_deconstructor_py_interface(self):
        """
        this function tests the quasi-Newton orbital updating deconstructor python
        interface
        """
        # initialize test flag
        test_passed = True

        # call quasi-Newton orbital updating deconstructor python interface
        update_orbs_qn_deconstructor()

        # check if deconstructor was called correctly
        if not c_bool.in_dll(lib, "test_update_orbs_qn_deconstructor_interface").value:
            print(
                " test_update_orbs_qn_deconstructor_py_interface failed: Called "
                "logging deconstructor wrong."
            )
            test_passed = False

        self.assertTrue(
            test_passed, "test_update_orbs_qn_deconstructor_py_interface failed"
        )
        print(" test_update_orbs_qn_deconstructor_py_interface PASSED")

    @patch.object(QNSettings, "init_c_struct", lib.mock_init_qn_settings)
    def test_qn_settings(self):
        """
        this function ensure the QNSettings object is properly initialized and
        synchronized with the underlying C struct
        """
        settings = QNSettings()
        test_passed = True
        for field_info in settings.c_struct._fields_:
            field_name, field_type = field_info[:2]
            if field_type == c_void_p:
                if (
                    getattr(settings, field_name) is not None
                    or getattr(settings.settings_c, field_name) is not None
                ):
                    print(
                        " test_qn_settings failed: Optional function pointer "
                        f"{field_name} not initialized correctly."
                    )
                    test_passed = False
            elif field_name == "initialized":
                if not getattr(settings, field_name):
                    print(
                        " test_qn_settings failed: Field initialized not "
                        "initialized correctly."
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
                        f" test_qn_settings failed: Field {field_name} not "
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
            print(" test_qn_settings failed: Optional callbacks are not set correctly.")
            test_passed = False

        self.assertTrue(test_passed, "test_qn_settings failed")
        print(" test_qn_settings PASSED")
