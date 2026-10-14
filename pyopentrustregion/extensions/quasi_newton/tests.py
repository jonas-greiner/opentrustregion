# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import unittest
from ctypes import c_bool
from unittest.mock import patch

import numpy as np

from pyopentrustregion.tests import (
    lib,
    n_param,
    add_tests,
    print_separator,
)
from pyopentrustregion.python_interface import (
    SolverSettings,
    SolverState,
    logger_interface_type,
)
from pyopentrustregion.extensions.quasi_newton import (
    QNSettings,
    qn_factory,
)
from pyopentrustregion.extensions.common.tests import (
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
    "qn_c_system_tests": ["qn_settings_init"],
    "qn_c_interface_tests": [
        "assign_qn_c_f",
        "assign_qn_f_c",
        "init_hess_f_wrapper",
        "init_qn_settings_c",
        "qn_factory_c_wrapper",
        "transport_f_wrapper",
    ],
    "qn_tests": [
        "add_qn",
        "hess_x_bfgs",
        "hess_x_sr1",
        "init_qn",
        "init_qn_settings",
        "qn_factory",
        "qn_sanity_check",
        "sr1_skip_add",
    ],
}


@add_tests
class QNTests(unittest.TestCase):
    """
    this class contains unit tests for quasi-Newton
    """

    tests = fortran_tests["qn_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for quasi-Newton...")
        return super().setUpClass()


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


class QNPyInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the Python interface
    """

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for quasi-Newton Python interface...")

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

    mock_logger = ExtensionCallbackMocks.mock_logger

    # replace original library with mock library, also for the release of the
    # extensions the factory returns
    @patch(
        "pyopentrustregion.python_interface.lib.qn_factory",
        lib.mock_qn_factory,
    )
    @patch(
        "pyopentrustregion.extensions.common.python_interface.release_extension",
        release_mock_extension,
    )
    def test_qn_factory_py_interface(self):
        """
        this function tests the quasi-Newton orbital updating factory python interface
        """
        # initialize test flag
        test_passed = True

        # initialize settings objects
        settings = QNSettings()
        assign_reference_settings(settings, lib.reference_qn_field)
        settings.logger = self.mock_logger
        solver_settings, python_functions = setup_wired_solver_settings()

        # initialize logging boolean
        self.test_logger = True

        # call quasi-Newton orbital updating factory python interface
        update_orbs_qn, obj_func_qn = qn_factory(
            n_param,
            self.mock_update_orbs,
            self.mock_obj_func,
            self.mock_transport,
            self.mock_init_hess,
            solver_settings,
            SolverState(),
            settings,
        )

        # check if logger was called correctly
        if not self.test_logger:
            print(
                " test_qn_factory_py_interface failed: Called logging "
                "function wrong."
            )
            test_passed = False

        # check that the settings were passed
        if not c_bool.in_dll(lib, "qn_settings_passed").value:
            print(" test_qn_factory_py_interface failed: Settings not passed.")
            test_passed = False

        # call returned quasi-Newton orbital updating function
        kappa = np.ones(n_param, dtype=np.float64)
        grad = np.empty(n_param, dtype=np.float64)
        h_diag = np.empty(n_param, dtype=np.float64)
        try:
            func, hess_x_funptr = update_orbs_qn(kappa, grad, h_diag)
        except RuntimeError:
            print(
                " test_qn_factory_py_interface failed: Returned "
                "quasi-Newton orbital updating raises error."
            )
            test_passed = False

        # check results
        if func != 3.0:
            print(
                " test_qn_factory_py_interface failed: Returned function "
                "value of returned quasi-Newton orbital updating function "
                "wrong."
            )
            test_passed = False
        if not np.allclose(grad, np.full(n_param, 2.0, dtype=np.float64)):
            print(
                " test_qn_factory_py_interface failed: Returned gradient "
                "of returned quasi-Newton orbital updating function wrong."
            )
            test_passed = False
        if not np.allclose(h_diag, np.full(n_param, 3.0, dtype=np.float64)):
            print(
                " test_qn_factory_py_interface failed: Returned "
                "quasi-Newton Hessian diagonal of returned quasi-Newton orbital "
                "updating function wrong."
            )
            test_passed = False

        # the returned Hessian linear transformation keeps the orbital updating function
        # alive, which keeps what the factory's object calls alive
        if hess_x_funptr._owner is not update_orbs_qn:
            print(
                " test_qn_factory_py_interface failed: Returned Hessian linear "
                "transformation does not keep the orbital updating function alive."
            )
            test_passed = False

        # call returned hess_x function
        x = np.ones(n_param, dtype=np.float64)
        hess_x = np.empty(n_param, dtype=np.float64)

        try:
            hess_x_funptr(x, hess_x)
        except RuntimeError:
            print(
                " test_qn_factory_py_interface failed: Returned "
                "quasi-Newton Hessian linear transformation function of returned "
                "quasi-Newton orbital updating function raises error."
            )
            test_passed = False

        # check results
        if not np.allclose(hess_x, np.full(n_param, 4.0, dtype=np.float64)):
            print(
                " test_qn_factory_py_interface failed: Returned Hessian "
                "linear transformation of returned quasi-Newton orbital updating "
                "function wrong."
            )
            test_passed = False

        # call returned objective function
        try:
            func = obj_func_qn(kappa)
        except RuntimeError:
            print(
                " test_qn_factory_py_interface failed: Returned objective function "
                "raises error."
            )
            test_passed = False
        else:
            if func != 3.0:
                print(
                    " test_qn_factory_py_interface failed: Returned function value of "
                    "returned objective function wrong."
                )
                test_passed = False

        # check the solver settings the factory wired
        if not check_wired_solver_settings(
            solver_settings,
            python_functions,
            "qn_factory_py_interface",
            (update_orbs_qn, obj_func_qn),
            settings,
        ):
            test_passed = False

        # call qn_factory python interface with the host's functions for the solver
        # settings the factory wired, which has to raise an error before the Fortran
        # factory is called
        try:
            qn_factory(
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
                " test_qn_factory_py_interface failed: Host's functions for wired "
                "solver settings did not raise an error."
            )
            test_passed = False
        except ValueError:
            pass

        # call qn_factory python interface without settings, which have to be passed as
        # none so that the library uses its default settings
        settings_flag = c_bool.in_dll(lib, "qn_settings_passed")
        settings_flag.value = True
        try:
            qn_factory(
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
                " test_qn_factory_py_interface failed: Call without settings raises "
                "error."
            )
            test_passed = False
        if settings_flag.value:
            print(
                " test_qn_factory_py_interface failed: Settings passed although none were "
                "given."
            )
            test_passed = False

        # the factory hands the extension it receives to an owner, which the solver
        # settings keep alive and which keeps the state alive
        if not check_extension_owner(
            lambda solver_settings, solver_state: qn_factory(
                n_param,
                self.mock_update_orbs,
                self.mock_obj_func,
                self.mock_transport,
                self.mock_init_hess,
                solver_settings,
                solver_state,
                settings,
            ),
            "qn_factory_py_interface",
        ):
            test_passed = False

        # call qn_factory python interface for host callbacks which cannot be called
        # with the arguments the C interface passes, which have to raise an error
        # before the Fortran factory is called
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
                qn_factory(
                    n_param,
                    self.mock_update_orbs,
                    self.mock_obj_func,
                    *callbacks,
                    SolverSettings(),
                    SolverState(),
                    settings,
                )
                print(
                    f" test_qn_factory_py_interface failed: {case} did not raise an "
                    "error."
                )
                test_passed = False
            except TypeError:
                pass

        self.assertTrue(
            c_bool.in_dll(lib, "test_qn_factory_interface").value and test_passed,
            "test_qn_factory_py_interface failed",
        )
        print(" test_qn_factory_py_interface PASSED")

    @patch.object(QNSettings, "init_c_struct", lib.mock_init_qn_settings)
    def test_qn_settings(self):
        """
        this function ensure the QNSettings object is properly initialized and
        synchronized with the underlying C struct
        """
        test_passed = True
        settings = QNSettings()
        if not settings_hold_reference(
            "test_qn_settings", settings, lib.reference_qn_field
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
            print(" test_qn_settings failed: Logging function not set correctly.")
            test_passed = False

        self.assertTrue(test_passed, "test_qn_settings failed")
        print(" test_qn_settings PASSED")


@add_tests
class QNCSystemTests(unittest.TestCase):
    """
    this class contains system tests for the quasi-Newton C interface
    """

    tests = fortran_tests["qn_c_system_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running system tests for quasi-Newton C interface...")
        return super().setUpClass()
