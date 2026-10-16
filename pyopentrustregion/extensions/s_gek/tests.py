# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import unittest
from ctypes import POINTER, c_bool, c_size_t
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
    SolverSettings,
    SolverState,
    logger_interface_type,
)
from pyopentrustregion.extensions.s_gek.python_interface import (
    ChangeReferenceInterface,
)
from pyopentrustregion.extensions.s_gek import (
    SGEKSettings,
    s_gek_factory,
)
from pyopentrustregion.extensions.common.tests import (
    check_settings_layout,
    setup_wired_solver_settings,
    check_wired_solver_settings,
    check_extension_owner,
    release_mock_extension,
    assign_reference_settings,
    settings_hold_reference,
    ExtensionCallbackMocks,
)

# define all tests in alphabetical order
fortran_tests = {
    "s_gek_c_integration_tests": ["s_gek_settings_init", "s_gek_settings_layout"],
    "s_gek_c_interface_tests": [
        "assign_s_gek_c_f",
        "assign_s_gek_f_c",
        "change_reference_f_wrapper",
        "init_s_gek_settings_c",
        "s_gek_factory_c_wrapper",
    ],
    "s_gek_tests": [
        "add_s_gek",
        "init_s_gek",
        "init_s_gek_settings",
        "s_gek_factory",
        "s_gek_sanity_check",
        "update_orbs_s_gek",
    ],
}


@add_tests
class SGEKTests(unittest.TestCase):
    """
    this class contains unit tests for S-GEK
    """

    tests = fortran_tests["s_gek_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for S-GEK...")
        return super().setUpClass()


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


class SGEKPyInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the Python interface
    """

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for S-GEK Python interface...")

        return super().setUpClass()

    mock_update_orbs = ExtensionCallbackMocks.mock_update_orbs
    mock_hess_x = ExtensionCallbackMocks.mock_hess_x
    mock_obj_func = ExtensionCallbackMocks.mock_obj_func

    def mock_change_reference(self, new_ref, kappa_list, local_grad_list, grad_list):
        """
        this function is a mock function for the change reference function, which
        changes the stored points to the new reference in a flat parameter space,
        subtracting it from every point and taking the gradients over from the local
        gradients, and triples the local gradients, as the change of reference checkers
        expect
        """
        kappa_list -= new_ref
        grad_list[:] = local_grad_list
        local_grad_list *= 3

    def test_callback_interfaces(self):
        """
        this function tests that the interface to the change of reference function
        reports an exception raised by the function as an error and keeps it
        """
        test_passed = True
        vector = np.zeros(n_param)
        vector_ptr = vector.ctypes.data_as(POINTER(c_real))

        def raising(*args):
            raise ValueError("callback function failure")

        for name, call in [
            (
                "change of reference function",
                lambda exc: ChangeReferenceInterface(raising, n_param, exc).call(
                    vector_ptr, 1, vector_ptr, vector_ptr, vector_ptr
                ),
            ),
        ]:
            exception = {}
            if call(exception) != 1:
                print(
                    f" test_callback_interfaces failed: Exception raised by {name} "
                    "not reported as error."
                )
                test_passed = False
            if not isinstance(exception.get("exc"), ValueError):
                print(
                    f" test_callback_interfaces failed: Exception raised by {name} "
                    "not kept."
                )
                test_passed = False

        self.assertTrue(test_passed, "test_callback_interfaces failed")
        print(" test_callback_interfaces PASSED")

    # replace original library with mock library, also for the release of the
    # extensions the factory returns
    @patch(
        "pyopentrustregion.python_interface.lib.s_gek_factory",
        lib.mock_s_gek_factory,
    )
    @patch(
        "pyopentrustregion.extensions.common.python_interface.release_extension",
        release_mock_extension,
    )
    def test_s_gek_factory_py_interface(self):
        """
        this function tests the S-GEK factory python interface
        """
        # initialize test flag
        test_passed = True

        # the logging function records the messages, and the mock factory records in a
        # flag whether it received what it expects, which is cleared before and read
        # after every call of it
        messages = []
        mock_passed = c_bool.in_dll(lib, "test_s_gek_factory_interface")

        # initialize settings objects
        settings = SGEKSettings()
        assign_reference_settings(settings, lib.reference_s_gek_field)
        settings.logger = messages.append
        solver_settings, python_functions = setup_wired_solver_settings()

        # call S-GEK factory python interface
        mock_passed.value = False
        update_orbs_s_gek, obj_func_s_gek = s_gek_factory(
            n_param,
            self.mock_update_orbs,
            self.mock_obj_func,
            self.mock_change_reference,
            solver_settings,
            SolverState(),
            settings,
        )

        # check if the mock factory received the correct input
        if not mock_passed.value:
            print(
                " test_s_gek_factory_py_interface failed: Mock factory received wrong "
                "input for settings."
            )
            test_passed = False

        # check if logger was called correctly
        if "test" not in messages:
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
        mock_passed.value = False
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

        # check if the mock factory received the correct input
        if not mock_passed.value:
            print(
                " test_s_gek_factory_py_interface failed: Mock factory received wrong "
                "input without settings."
            )
            test_passed = False

        # the factory hands the extension it receives to an owner, which the solver
        # settings keep alive and which keeps the state alive
        mock_passed.value = False
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

        # check if the mock factory received the correct input
        if not mock_passed.value:
            print(
                " test_s_gek_factory_py_interface failed: Mock factory received wrong "
                "input for the check of the owner."
            )
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

        self.assertTrue(test_passed, "test_s_gek_factory_py_interface failed")
        print(" test_s_gek_factory_py_interface PASSED")

    @patch.object(SGEKSettings, "init_c_struct", lib.mock_init_s_gek_settings)
    def test_s_gek_settings(self):
        """
        this function ensure the SGEKSettings object is properly initialized and
        synchronized with the underlying C struct
        """
        test_passed = True
        settings = SGEKSettings()
        if not settings_hold_reference(
            "test_s_gek_settings", settings, lib.reference_s_gek_field
        ):
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
class SGEKCIntegrationTests(unittest.TestCase):
    """
    this class contains integration tests that drive the real library through the
    S-GEK C header
    """

    tests = fortran_tests["s_gek_c_integration_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running integration tests for S-GEK C interface...")
        return super().setUpClass()


class SGEKPyIntegrationTests(unittest.TestCase):
    """
    this class contains integration tests that drive the real library (not the mock
    library) through the S-GEK Python interface
    """

    @classmethod
    def setUpClass(cls):
        print_separator("Running integration tests for S-GEK Python interface...")
        return super().setUpClass()

    def test_s_gek_settings_layout(self):
        """
        this function checks that the Python S-GEK settings structure has the size of
        the bind(C) settings type and that every field is read back under its own name
        from settings whose fields Fortran sets one by one to the reference values
        """
        lib.s_gek_settings_size.restype = c_size_t
        self.assertTrue(
            check_settings_layout(
                "test_s_gek_settings_layout",
                SGEKSettings,
                lib.get_reference_s_gek_values,
                lib.reference_s_gek_field,
                lib.s_gek_settings_size(),
            ),
            "test_s_gek_settings_layout failed",
        )
        print(" test_s_gek_settings_layout PASSED")
