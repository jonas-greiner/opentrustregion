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
)
from pyopentrustregion.extensions.rotation_qn.python_interface import (
    TransportInterface,
    InitHessInterface,
)
from pyopentrustregion.extensions.rotation_qn import (
    RotationQNSettings,
    rotation_qn_factory,
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
    "rotation_qn_c_integration_tests": [
        "rotation_qn_settings_init",
        "rotation_qn_settings_layout",
    ],
    "rotation_qn_c_interface_tests": [
        "assign_rotation_qn_c_f",
        "assign_rotation_qn_f_c",
        "init_hess_f_wrapper",
        "init_rotation_qn_settings_c",
        "rotation_qn_factory_c_wrapper",
        "transport_f_wrapper",
    ],
    "rotation_qn_tests": [
        "add_rotation_qn",
        "hess_x_bfgs",
        "hess_x_sr1",
        "init_rotation_qn",
        "rotation_qn_factory",
        "rotation_qn_sanity_check",
        "sr1_skip_add",
        "update_orbs_rotation_qn",
    ],
}


@add_tests
class RotationQNTests(unittest.TestCase):
    """
    this class contains unit tests for rotation QN
    """

    tests = fortran_tests["rotation_qn_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for rotation QN...")
        return super().setUpClass()


@add_tests
class RotationQNCInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the rotation QN C interface
    """

    tests = fortran_tests["rotation_qn_c_interface_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for rotation QN C interface...")
        return super().setUpClass()


class RotationQNPyInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the Python interface
    """

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for rotation QN Python interface...")

        return super().setUpClass()

    mock_update_orbs = ExtensionCallbackMocks.mock_update_orbs
    mock_hess_x = ExtensionCallbackMocks.mock_hess_x
    mock_obj_func = ExtensionCallbackMocks.mock_obj_func

    def mock_transport(self, geodesic, tangent_vector):
        """
        this function is a mock function for the transport function
        """
        tangent_vector[:] = geodesic * tangent_vector

        return

    def mock_init_hess(self, vector):
        """
        this function is a mock function for the initial Hessian function
        """
        vector[:] = 2 * vector

        return

    def test_callback_interfaces(self):
        """
        this function tests that the interfaces to the transport and initial Hessian
        functions report an exception raised by the function as an error and keep it
        """
        test_passed = True
        vector = np.zeros(n_param)
        vector_ptr = vector.ctypes.data_as(POINTER(c_real))

        def raising(*args):
            raise ValueError("callback function failure")

        for name, call in [
            (
                "transport function",
                lambda exc: TransportInterface(raising, n_param, exc).call(
                    vector_ptr, vector_ptr
                ),
            ),
            (
                "initial Hessian function",
                lambda exc: InitHessInterface(raising, n_param, exc).call(vector_ptr),
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
        "pyopentrustregion.python_interface.lib.rotation_qn_factory",
        lib.mock_rotation_qn_factory,
    )
    @patch(
        "pyopentrustregion.extensions.common.python_interface.release_extension",
        release_mock_extension,
    )
    def test_rotation_qn_factory_py_interface(self):
        """
        this function tests the rotation QN factory python interface
        """
        # initialize test flag
        test_passed = True

        # the logging function records the messages, and the mock factory records in a
        # flag whether it received what it expects, which is cleared before and read
        # after every call of it
        messages = []
        mock_passed = c_bool.in_dll(lib, "test_rotation_qn_factory_interface")

        # initialize settings objects
        settings = RotationQNSettings()
        assign_reference_settings(settings, lib.reference_rotation_qn_field)
        settings.logger = messages.append
        solver_settings, python_functions = setup_wired_solver_settings()

        # call rotation QN factory python interface
        mock_passed.value = False
        update_orbs_rotation_qn, obj_func_rotation_qn = rotation_qn_factory(
            n_param,
            self.mock_update_orbs,
            self.mock_obj_func,
            self.mock_transport,
            self.mock_init_hess,
            solver_settings,
            SolverState(),
            settings,
        )

        # check if the mock factory received the correct input
        if not mock_passed.value:
            print(
                " test_rotation_qn_factory_py_interface failed: Mock factory received "
                "wrong "
                "input for settings."
            )
            test_passed = False

        # check if logger was called correctly
        if "test" not in messages:
            print(
                " test_rotation_qn_factory_py_interface failed: Called logging "
                "function wrong."
            )
            test_passed = False

        # check that the settings were passed
        if not c_bool.in_dll(lib, "rotation_qn_settings_passed").value:
            print(" test_rotation_qn_factory_py_interface failed: Settings not passed.")
            test_passed = False

        # call returned rotation QN orbital updating function
        kappa = np.ones(n_param, dtype=np.float64)
        grad = np.empty(n_param, dtype=np.float64)
        h_diag = np.empty(n_param, dtype=np.float64)
        try:
            func, hess_x_funptr = update_orbs_rotation_qn(kappa, grad, h_diag)
        except RuntimeError:
            print(
                " test_rotation_qn_factory_py_interface failed: Returned "
                "rotation QN orbital updating raises error."
            )
            test_passed = False

        # check results
        if func != 3.0:
            print(
                " test_rotation_qn_factory_py_interface failed: Returned function "
                "value of returned rotation QN orbital updating function "
                "wrong."
            )
            test_passed = False
        if not np.allclose(grad, np.full(n_param, 2.0, dtype=np.float64)):
            print(
                " test_rotation_qn_factory_py_interface failed: Returned gradient "
                "of returned rotation QN orbital updating function wrong."
            )
            test_passed = False
        if not np.allclose(h_diag, np.full(n_param, 3.0, dtype=np.float64)):
            print(
                " test_rotation_qn_factory_py_interface failed: Returned "
                "quasi-Newton Hessian diagonal of returned rotation QN orbital "
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
                " test_rotation_qn_factory_py_interface failed: Returned "
                "quasi-Newton Hessian linear transformation function of returned "
                "rotation QN orbital updating function raises error."
            )
            test_passed = False

        # check results
        if not np.allclose(hess_x, np.full(n_param, 4.0, dtype=np.float64)):
            print(
                " test_rotation_qn_factory_py_interface failed: Returned Hessian "
                "linear transformation of returned rotation QN orbital updating "
                "function wrong."
            )
            test_passed = False

        # call returned objective function
        try:
            func = obj_func_rotation_qn(kappa)
        except RuntimeError:
            print(
                " test_rotation_qn_factory_py_interface failed: Returned objective "
                "function "
                "raises error."
            )
            test_passed = False
        else:
            if func != 3.0:
                print(
                    " test_rotation_qn_factory_py_interface failed: Returned function "
                    "value of "
                    "returned objective function wrong."
                )
                test_passed = False

        # check the solver settings the factory wired
        if not check_wired_solver_settings(
            solver_settings,
            python_functions,
            "test_rotation_qn_factory_py_interface",
            (update_orbs_rotation_qn, obj_func_rotation_qn),
            settings,
        ):
            test_passed = False

        # call rotation_qn_factory python interface with the host's functions for the
        # solver settings the factory wired, which has to raise an error before the
        # Fortran factory is called
        try:
            rotation_qn_factory(
                n_param,
                self.mock_update_orbs,
                self.mock_obj_func,
                self.mock_transport,
                self.mock_init_hess,
                solver_settings,
                SolverState(),
                settings,
            )
            print(
                " test_rotation_qn_factory_py_interface failed: Host's functions for "
                "wired "
                "solver settings did not raise an error."
            )
            test_passed = False
        except ValueError:
            pass

        # call rotation_qn_factory python interface without settings, which have to be
        # passed as none so that the library uses its default settings
        settings_flag = c_bool.in_dll(lib, "rotation_qn_settings_passed")
        settings_flag.value = True
        mock_passed.value = False
        try:
            rotation_qn_factory(
                n_param,
                self.mock_update_orbs,
                self.mock_obj_func,
                self.mock_transport,
                self.mock_init_hess,
                setup_wired_solver_settings()[0],
                SolverState(),
            )
        except RuntimeError:
            print(
                " test_rotation_qn_factory_py_interface failed: Call without settings "
                "raises "
                "error."
            )
            test_passed = False
        if settings_flag.value:
            print(
                " test_rotation_qn_factory_py_interface failed: Settings passed "
                "although none were "
                "given."
            )
            test_passed = False

        # check if the mock factory received the correct input
        if not mock_passed.value:
            print(
                " test_rotation_qn_factory_py_interface failed: Mock factory received "
                "wrong "
                "input without settings."
            )
            test_passed = False

        # the factory hands the extension it receives to an owner, which the solver
        # settings keep alive and which keeps the state alive
        mock_passed.value = False
        if not check_extension_owner(
            lambda solver_settings, solver_state: rotation_qn_factory(
                n_param,
                self.mock_update_orbs,
                self.mock_obj_func,
                self.mock_transport,
                self.mock_init_hess,
                solver_settings,
                solver_state,
                settings,
            ),
            "test_rotation_qn_factory_py_interface",
        ):
            test_passed = False

        # check if the mock factory received the correct input
        if not mock_passed.value:
            print(
                " test_rotation_qn_factory_py_interface failed: Mock factory received "
                "wrong "
                "input for the check of the owner."
            )
            test_passed = False

        # call rotation_qn_factory python interface for host callbacks which cannot be
        # called with the arguments the C interface passes, which have to raise an
        # error before the Fortran factory is called
        invalid_cases = (
            (
                "transport with a wrong number of arguments",
                (lambda geodesic: None, self.mock_init_hess),
            ),
            (
                "init_hess with a wrong number of arguments",
                (self.mock_transport, lambda: None),
            ),
        )
        for case, callbacks in invalid_cases:
            try:
                rotation_qn_factory(
                    n_param,
                    self.mock_update_orbs,
                    self.mock_obj_func,
                    *callbacks,
                    SolverSettings(),
                    SolverState(),
                    settings,
                )
                print(
                    f" test_rotation_qn_factory_py_interface failed: {case} did not "
                    f"raise an "
                    "error."
                )
                test_passed = False
            except TypeError:
                pass

        self.assertTrue(test_passed, "test_rotation_qn_factory_py_interface failed")
        print(" test_rotation_qn_factory_py_interface PASSED")

    @patch.object(
        RotationQNSettings, "init_c_struct", lib.mock_init_rotation_qn_settings
    )
    def test_rotation_qn_settings(self):
        """
        this function ensures the RotationQNSettings object is properly initialized and
        synchronized with the underlying C struct
        """
        self.assertTrue(
            settings_hold_reference(
                "test_rotation_qn_settings",
                RotationQNSettings(),
                lib.reference_rotation_qn_field,
            ),
            "test_rotation_qn_settings failed",
        )
        print(" test_rotation_qn_settings PASSED")


@add_tests
class RotationQNCIntegrationTests(unittest.TestCase):
    """
    this class contains integration tests that drive the real library through the
    rotation QN C header
    """

    tests = fortran_tests["rotation_qn_c_integration_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running integration tests for rotation QN C interface...")
        return super().setUpClass()


class RotationQNPyIntegrationTests(unittest.TestCase):
    """
    this class contains integration tests that drive the real library (not the mock
    library) through the rotation QN Python interface
    """

    @classmethod
    def setUpClass(cls):
        print_separator("Running integration tests for rotation QN Python interface...")
        return super().setUpClass()

    def test_rotation_qn_settings_layout(self):
        """
        this function checks that the Python rotation QN settings structure has the
        size of the bind(C) settings type and that every field is read back under its
        own name from settings whose fields Fortran sets one by one to the reference
        values
        """
        lib.rotation_qn_settings_size.restype = c_size_t
        self.assertTrue(
            check_settings_layout(
                "test_rotation_qn_settings_layout",
                RotationQNSettings,
                lib.get_reference_rotation_qn_values,
                lib.reference_rotation_qn_field,
                lib.rotation_qn_settings_size(),
            ),
            "test_rotation_qn_settings_layout failed",
        )
        print(" test_rotation_qn_settings_layout PASSED")
