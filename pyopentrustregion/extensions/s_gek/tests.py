# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import unittest
from ctypes import (
    c_bool,
    c_char,
    Array,
    byref,
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
    LegacyPyInterfaceHelpers,
)
from pyopentrustregion.python_interface import (
    c_real,
    c_int,
    SolverSettings,
    SolverState,
    logger_interface_type,
)
from pyopentrustregion.extensions.s_gek import (
    SGEKSettings,
    s_gek_factory,
)
from pyopentrustregion.extensions.common.tests import (
    setup_wired_solver_settings,
    check_wired_solver_settings,
    check_extension_owner,
)

if NUMPY_AVAILABLE:
    import numpy as np


# define all tests in alphabetical order
fortran_tests = {
    "s_gek_c_system_tests": ["s_gek_settings_init"],
    "s_gek_c_interface_tests": [
        "assign_s_gek_c_f",
        "assign_s_gek_f_c",
        "change_reference_f_wrapper",
        "init_s_gek_settings_c",
        "s_gek_factory_c_wrapper",
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

    assign_ref_to_settings = LegacyPyInterfaceHelpers.assign_ref_to_settings
    equal_settings_to_ref = LegacyPyInterfaceHelpers.equal_settings_to_ref

    mock_update_orbs = LegacyPyInterfaceHelpers.mock_update_orbs
    mock_hess_x = LegacyPyInterfaceHelpers.mock_hess_x
    mock_obj_func = LegacyPyInterfaceHelpers.mock_obj_func

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

    mock_logger = LegacyPyInterfaceHelpers.mock_logger

    # replace original library with mock library
    @patch(
        "pyopentrustregion.python_interface.lib.s_gek_factory",
        lib.mock_s_gek_factory,
    )
    def test_s_gek_factory_py_interface(self):
        """
        this function tests the S-GEK orbital updating factory python interface
        """
        # initialize test flag
        test_passed = True

        # initialize settings objects
        settings = SGEKSettings()
        self.assign_ref_to_settings(settings)
        solver_settings, python_functions = setup_wired_solver_settings()

        # initialize logging boolean
        self.test_logger = True

        # call S-GEK orbital updating factory python interface
        update_orbs_s_gek, obj_func_s_gek = s_gek_factory(
            n_param,
            self.mock_update_orbs,
            self.mock_obj_func,
            self.mock_change_reference,
            solver_settings,
            SolverState(),
            settings,
        )

        # check if logger was called correctly
        if not self.test_logger:
            print(
                " test_s_gek_factory_py_interface failed: Called logging "
                "function wrong."
            )
            test_passed = False

        # check that the settings were passed
        if not c_bool.in_dll(lib, "s_gek_settings_passed").value:
            print(" test_s_gek_factory_py_interface failed: Settings not passed.")
            test_passed = False

        # call returned S-GEK orbital updating function
        kappa = np.ones(n_param, dtype=np.float64)
        grad = np.empty(n_param, dtype=np.float64)
        h_diag = np.empty(n_param, dtype=np.float64)
        try:
            func, hess_x_funptr = update_orbs_s_gek(kappa, grad, h_diag)
        except RuntimeError:
            print(
                " test_s_gek_factory_py_interface failed: Returned S-GEK "
                "orbital updating raises error."
            )
            test_passed = False

        # check results
        if self.change_reference_new_ref is not True:
            print(
                " test_s_gek_factory_py_interface failed: New "
                "reference parameters inside given change reference function wrong."
            )
            test_passed = False
        if func != 3.0:
            print(
                " test_s_gek_factory_py_interface failed: Returned"
                " function value of returned S-GEK orbital updating function wrong."
            )
            test_passed = False
        if not np.allclose(grad, np.full(n_param, 2.0, dtype=np.float64)):
            print(
                " test_s_gek_factory_py_interface failed: Returned "
                "gradient of returned S-GEK orbital updating function wrong."
            )
            test_passed = False
        if not np.allclose(h_diag, np.full(n_param, 3.0, dtype=np.float64)):
            print(
                " test_s_gek_factory_py_interface failed: Returned "
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
                " test_s_gek_factory_py_interface failed: Returned "
                "S-GEK Hessian linear transformation function of returned S-GEK "
                "orbital updating function raises error."
            )
            test_passed = False

        # check results
        if not np.allclose(hess_x, np.full(n_param, 4.0, dtype=np.float64)):
            print(
                " test_s_gek_factory_py_interface failed: Returned Hessian "
                "linear transformation of returned S-GEK orbital updating function "
                "wrong."
            )
            test_passed = False

        # call returned objective function
        try:
            func = obj_func_s_gek(kappa)
        except RuntimeError:
            print(
                " test_s_gek_factory_py_interface failed: Returned objective "
                "function raises error."
            )
            test_passed = False
        else:
            if func != 3.0:
                print(
                    " test_s_gek_factory_py_interface failed: Returned function value "
                    "of returned objective function wrong."
                )
                test_passed = False

        # check the solver settings the factory wired
        if not check_wired_solver_settings(
            solver_settings,
            python_functions,
            "s_gek_factory_py_interface",
            (update_orbs_s_gek, obj_func_s_gek),
            settings,
        ):
            test_passed = False

        # call s_gek_factory python interface with the host's functions for the solver
        # settings the factory wired, which has to raise an error before the Fortran
        # factory is called
        try:
            s_gek_factory(
                n_param,
                self.mock_update_orbs,
                self.mock_obj_func,
                self.mock_change_reference,
                solver_settings,
                SolverState(),
                settings,
            )
            print(
                " test_s_gek_factory_py_interface failed: Host's functions for wired "
                "solver settings did not raise an error."
            )
            test_passed = False
        except ValueError:
            pass

        # call s_gek_factory python interface without settings, which have to be passed as
        # none so that the library uses its default settings
        settings_flag = c_bool.in_dll(lib, "s_gek_settings_passed")
        settings_flag.value = True
        try:
            s_gek_factory(
                n_param,
                self.mock_update_orbs,
                self.mock_obj_func,
                self.mock_change_reference,
                setup_wired_solver_settings()[0],
                SolverState(),
            )
        except RuntimeError:
            print(
                " test_s_gek_factory_py_interface failed: Call without settings raises "
                "error."
            )
            test_passed = False
        if settings_flag.value:
            print(
                " test_s_gek_factory_py_interface failed: Settings passed although none were "
                "given."
            )
            test_passed = False

        # the factory hands the extension it receives to an owner, which the solver
        # settings keep alive and which keeps the state alive
        if not check_extension_owner(
            lambda solver_settings, solver_state: s_gek_factory(
                n_param,
                self.mock_update_orbs,
                self.mock_obj_func,
                self.mock_change_reference,
                solver_settings,
                solver_state,
                settings,
            ),
            "s_gek_factory_py_interface",
        ):
            test_passed = False

        # call s_gek_factory python interface for a host callback which cannot be
        # called with the arguments the C interface passes, which has to raise an error
        # before the Fortran factory is called
        try:
            s_gek_factory(
                n_param,
                self.mock_update_orbs,
                self.mock_obj_func,
                lambda new_ref, kappa_list, local_grad_list: None,
                SolverSettings(),
                SolverState(),
                settings,
            )
            print(
                " test_s_gek_factory_py_interface failed: change_reference with a "
                "wrong number of arguments did not raise an error."
            )
            test_passed = False
        except TypeError:
            pass

        self.assertTrue(
            c_bool.in_dll(lib, "test_s_gek_factory_interface").value and test_passed,
            "test_s_gek_factory_py_interface failed",
        )
        print(" test_s_gek_factory_py_interface PASSED")

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
            test_passed = False

        # assigning a logging function writes the C callback function the settings
        # built for it into the C settings, which calls it whatever context it
        # receives, and reading returns the logging function
        messages = []

        def logger(message):
            messages.append(message)

        settings.logger = logger
        logger_interface_type(settings.settings_c.logger)(b"test", 12345)
        if settings.logger is not logger or messages != ["test"]:
            print(" test_s_gek_settings failed: Logging function not set correctly.")
            test_passed = False

        self.assertTrue(test_passed, "test_s_gek_settings failed")
        print(" test_s_gek_settings PASSED")


@add_tests
class SGEKCSystemTests(unittest.TestCase):
    """
    this class contains system tests for the S-GEK C interface
    """

    tests = fortran_tests["s_gek_c_system_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running system tests for S-GEK C interface...")
        return super().setUpClass()
