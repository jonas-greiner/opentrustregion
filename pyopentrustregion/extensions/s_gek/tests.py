# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import unittest
from ctypes import POINTER, c_bool, c_void_p
from unittest.mock import patch

import numpy as np

from pyopentrustregion.tests import (
    lib,
    n_param,
    add_tests,
    print_separator,
    check_callback_interfaces,
    reference_state,
    Hartmann6DProblem,
)
from pyopentrustregion.python_interface import (
    c_int,
    c_real,
    hess_x_interface_type,
    SolverSettings,
    SolverState,
)
from pyopentrustregion.extensions.s_gek.python_interface import (
    ChangeReferenceInterface,
)
from pyopentrustregion.extensions.common.python_interface import (
    obj_func_grad_interface_type,
)
from pyopentrustregion.extensions.s_gek import (
    SGEKSettings,
    s_gek_solver,
)
from pyopentrustregion.extensions.common.tests import (
    check_settings_layout,
    assign_reference_settings,
    settings_hold_reference,
    ExtensionCallbackMocks,
)

# define all tests in alphabetical order
fortran_tests = {
    "s_gek_c_integration_tests": [
        "s_gek_settings_init",
        "s_gek_settings_layout",
        "s_gek_solver_c",
    ],
    "s_gek_c_interface_tests": [
        "assign_s_gek_c_f",
        "assign_s_gek_f_c",
        "change_reference_f_wrapper",
        "init_s_gek_settings_c",
        "s_gek_solver_c_wrapper",
    ],
    "s_gek_tests": [
        "conv_check_surrogate",
        "dispersion_gek",
        "distances_gek",
        "energy_gek",
        "floored_h_diag",
        "gradient_gek",
        "hess_x_surrogate",
        "hessian_gek",
        "inner_converged",
        "inverse_bfgs_step",
        "invoke_change_reference",
        "matern_kernel",
        "minimize_surrogate_otr",
        "minimize_surrogate_rs_rfo",
        "obj_func_surrogate",
        "rs_rfo",
        "s_gek_sanity_check",
        "s_gek_solver",
        "s_gek_step",
        "setup_gek",
        "surrogate_of",
        "update_orbs_surrogate",
        "updated_prior",
        "updated_threshold_factor",
        "variance_threshold_of",
    ],
}

# define the signature of the Fortran function providing the reference settings
lib.get_reference_s_gek_values.argtypes = [POINTER(SGEKSettings.c_struct)]
lib.get_reference_s_gek_values.restype = None


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
        this function is a mock function for the chart transformation, which changes
        the points to the new reference in a flat parameter space, subtracting it from
        every point and taking the gradients over from the local gradients, and triples
        the local gradients, as the chart transformation checkers expect
        """
        kappa_list -= new_ref
        grad_list[:] = local_grad_list
        local_grad_list *= 3

    def mock_obj_func_grad(self, kappa, grad):
        """
        this function is a mock function for the objective function of the settings,
        which returns the sum of the step as the energy and twice the step as the
        gradient
        """
        grad[:] = 2 * kappa

        return np.sum(kappa)

    def test_callback_interfaces(self):
        """
        this function tests that the interface to the chart transformation reports an
        exception or an interruption raised by the function as an error and keeps it
        with the call whose library called it, and that the function receives the new
        reference read-only
        """
        vector = np.zeros(n_param)
        vector_ptr = vector.ctypes.data_as(POINTER(c_real))
        output = np.zeros(n_param)
        output_ptr = output.ctypes.data_as(POINTER(c_real))
        interfaces = {
            "chart transformation": lambda func: ChangeReferenceInterface(func)(
                vector_ptr, 1, output_ptr, output_ptr, output_ptr, None
            ),
        }
        test_passed = check_callback_interfaces(
            "test_callback_interfaces", interfaces, list(interfaces), vector, n_param
        )

        self.assertTrue(test_passed, "test_callback_interfaces failed")
        print(" test_callback_interfaces PASSED")

    # replace original library with mock library
    @patch("pyopentrustregion.python_interface.lib.s_gek_solver", lib.mock_s_gek_solver)
    def test_s_gek_solver_py_interface(self):
        """
        this function tests the S-GEK solver python interface
        """
        # initialize test flag
        test_passed = True

        # the logging function records the messages, and the mock S-GEK solver records
        # in a flag whether it received what it expects, which is cleared before and
        # read after every call of it
        messages = []
        mock_passed = c_bool.in_dll(lib, "test_s_gek_solver_interface")
        settings_flag = c_bool.in_dll(lib, "s_gek_solver_settings_passed")
        s_gek_settings_flag = c_bool.in_dll(lib, "s_gek_settings_passed")

        # initialize S-GEK settings with the reference values and every callback
        # function
        s_gek_settings = SGEKSettings()
        assign_reference_settings(s_gek_settings, lib.reference_s_gek_field)
        s_gek_settings.obj_func_grad = self.mock_obj_func_grad
        s_gek_settings.init_hess = self.mock_hess_x
        s_gek_settings.logger = messages.append

        # call S-GEK solver python interface with the optional arguments
        mock_passed.value = False
        state = SolverState()
        s_gek_solver(
            self.mock_update_orbs,
            self.mock_obj_func,
            self.mock_change_reference,
            n_param,
            SolverSettings(),
            state,
            s_gek_settings,
        )
        if not mock_passed.value:
            print(
                " test_s_gek_solver_py_interface failed: Mock S-GEK solver received "
                "wrong input."
            )
            test_passed = False
        if not (settings_flag.value and s_gek_settings_flag.value):
            print(" test_s_gek_solver_py_interface failed: Settings not passed.")
            test_passed = False
        if "test" not in messages:
            print(
                " test_s_gek_solver_py_interface failed: Called logging function "
                "wrong."
            )
            test_passed = False

        # check that the state is passed, which the mock fills with the reference
        # values
        if state.n_update_orbs != reference_state("n_update_orbs", c_int):
            print(" test_s_gek_solver_py_interface failed: State not passed.")
            test_passed = False

        # a callback function that raises makes the S-GEK solver fail with the
        # exception as the cause of the reported error
        def raising_obj_func(kappa):
            raise ValueError("objective function failure")

        try:
            s_gek_solver(
                self.mock_update_orbs,
                raising_obj_func,
                self.mock_change_reference,
                n_param,
            )
            print(
                " test_s_gek_solver_py_interface failed: Exception raised by "
                "objective function was not reported."
            )
            test_passed = False
        except RuntimeError as e:
            if not isinstance(e.__cause__, ValueError):
                print(
                    " test_s_gek_solver_py_interface failed: Exception raised by "
                    "objective function was not the cause of the reported error."
                )
                test_passed = False

        # call S-GEK solver python interface without settings, which have to be passed
        # as none so that the library uses its default settings, the results of the
        # mock start opposite to the expected ones
        settings_flag.value = True
        s_gek_settings_flag.value = True
        mock_passed.value = False
        s_gek_solver(
            self.mock_update_orbs,
            self.mock_obj_func,
            self.mock_change_reference,
            n_param,
        )
        if not mock_passed.value:
            print(
                " test_s_gek_solver_py_interface failed: Mock S-GEK solver received "
                "wrong input without settings."
            )
            test_passed = False
        if settings_flag.value or s_gek_settings_flag.value:
            print(
                " test_s_gek_solver_py_interface failed: Settings passed although "
                "none were given."
            )
            test_passed = False

        self.assertTrue(test_passed, "test_s_gek_solver_py_interface failed")
        print(" test_s_gek_solver_py_interface PASSED")

    @patch.object(SGEKSettings, "init_c_struct", lib.mock_init_s_gek_settings)
    def test_s_gek_settings(self):
        """
        this function ensures the SGEKSettings object is properly initialized and
        synchronized with the underlying C struct, and that the C functions of every
        callback field in the C settings are read as callables around them
        """
        test_passed = settings_hold_reference(
            "test_s_gek_settings", SGEKSettings(), lib.reference_s_gek_field
        )

        # the C objective function and initial Hessian in the C settings are read as
        # callables around them, which hand them their inputs and the context of the
        # settings, write their outputs back, also into arrays the library cannot write
        # in place, and raise their errors
        received_contexts = []
        error_code = 0

        def c_obj_func_grad(kappa_ptr, grad_ptr, func_ptr, context):
            kappa = np.ctypeslib.as_array(kappa_ptr, shape=(n_param,))
            np.ctypeslib.as_array(grad_ptr, shape=(n_param,))[:] = 2.0 * kappa
            func_ptr[0] = np.sum(kappa)
            received_contexts.append(context)
            return error_code

        def c_init_hess(x_ptr, hess_x_ptr, context):
            x = np.ctypeslib.as_array(x_ptr, shape=(n_param,))
            np.ctypeslib.as_array(hess_x_ptr, shape=(n_param,))[:] = 3.0 * x
            received_contexts.append(context)
            return error_code

        # outputs, as strided views, and the values the C functions have to write into
        # them for the input x
        x = np.random.random(n_param)

        def call_obj_func_grad(function):
            output = np.zeros(2 * n_param)[::2]
            return [np.array(function(x, output)), output]

        def call_init_hess(function):
            output = np.zeros(2 * n_param)[::2]
            function(x, output)
            return [output]

        for name, c_function, call, expected in (
            (
                "obj_func_grad",
                obj_func_grad_interface_type(c_obj_func_grad),
                call_obj_func_grad,
                [np.sum(x), 2.0 * x],
            ),
            (
                "init_hess",
                hess_x_interface_type(c_init_hess),
                call_init_hess,
                [3.0 * x],
            ),
        ):
            settings = SGEKSettings()
            settings.settings_c.context = 6789
            setattr(settings.settings_c, name, c_void_p.from_buffer(c_function).value)
            received_contexts = []
            error_code = 0
            outputs = call(getattr(settings, name))
            if not all(
                np.allclose(output, values) for output, values in zip(outputs, expected)
            ) or received_contexts != [6789]:
                print(
                    f" test_s_gek_settings failed: C function of {name} not called "
                    "with the inputs and the context of the settings or its outputs "
                    "not returned."
                )
                test_passed = False
            error_code = 3
            try:
                call(getattr(settings, name))
                print(
                    f" test_s_gek_settings failed: Error of the C function of {name} "
                    "not raised."
                )
                test_passed = False
            except RuntimeError:
                pass

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


class SGEKPyIntegrationTests(Hartmann6DProblem, unittest.TestCase):
    """
    this class contains integration tests that drive the real library (not the mock
    library) through the S-GEK Python interface
    """

    @classmethod
    def setUpClass(cls):
        print_separator("Running integration tests for S-GEK Python interface...")
        cls.load_hartmann6d()
        return super().setUpClass()

    def test_s_gek_solver_py(self):
        """
        this function solves the Hartmann 6D problem with the real S-GEK solver through
        the Python interface in its flat chart, with every callback function of the
        S-GEK settings, the exact Hessian as the initial Hessian, and checks that the
        solve converges to the minimum, that every callback function is called and
        finds the state current and that the logging function of the S-GEK settings
        receives the error messages of S-GEK
        """
        test_passed = True
        curr_vars = self.near_minimum.copy()
        curr_hess = np.empty((self.n_param, self.n_param))
        state = SolverState()
        counts = dict.fromkeys(
            ("update_orbs", "change_reference", "obj_func_grad", "init_hess"), 0
        )
        messages = []
        state_not_current = False

        def check_state_current():
            nonlocal state_not_current
            if state.n_update_orbs != counts["update_orbs"]:
                state_not_current = True

        def exact_hess_x(x, hess_x):
            hess_x[:] = curr_hess @ x

        def update_orbs(kappa, grad, h_diag):
            check_state_current()
            counts["update_orbs"] += 1
            curr_vars[:] += kappa
            grad[:] = self._grad(curr_vars)
            curr_hess[:] = self._hess(curr_vars)
            h_diag[:] = np.diag(curr_hess)
            return self._func(curr_vars), exact_hess_x

        def obj_func(kappa):
            return self._func(curr_vars + kappa)

        def change_reference(new_ref, kappa_list, local_grad_list, grad_list):
            check_state_current()
            counts["change_reference"] += 1
            kappa_list -= new_ref
            grad_list[:] = local_grad_list

        def obj_func_grad(kappa, grad):
            nonlocal state_not_current
            check_state_current()
            if state.n_obj_func != counts["obj_func_grad"]:
                state_not_current = True
            counts["obj_func_grad"] += 1
            grad[:] = self._grad(curr_vars + kappa)
            return self._func(curr_vars + kappa)

        def init_hess(x, hess_x):
            check_state_current()
            counts["init_hess"] += 1
            exact_hess_x(x, hess_x)

        s_gek_settings = SGEKSettings()
        s_gek_settings.obj_func_grad = obj_func_grad
        s_gek_settings.init_hess = init_hess
        s_gek_settings.logger = messages.append
        s_gek_solver(
            update_orbs,
            obj_func,
            change_reference,
            self.n_param,
            SolverSettings(),
            state,
            s_gek_settings,
        )
        if not np.allclose(curr_vars, self.minimum1, rtol=0.0, atol=1e-4):
            print(
                " test_s_gek_solver_py failed: Solve did not converge to the minimum."
            )
            test_passed = False
        if any(count == 0 for count in counts.values()):
            print(" test_s_gek_solver_py failed: Callback functions not called.")
            test_passed = False
        if (
            state_not_current
            or state.n_update_orbs != counts["update_orbs"]
            or state.n_obj_func != counts["obj_func_grad"]
        ):
            print(
                " test_s_gek_solver_py failed: State was not current in the callbacks "
                "or not returned."
            )
            test_passed = False

        # the logging function of the S-GEK settings receives the error messages of
        # S-GEK, here for S-GEK settings with an option S-GEK does not know
        s_gek_settings.inner_solver = "unknown"
        s_gek_settings.verbose = 1
        messages.clear()
        try:
            s_gek_solver(
                update_orbs,
                obj_func,
                change_reference,
                self.n_param,
                SolverSettings(),
                s_gek_settings=s_gek_settings,
            )
            print(
                " test_s_gek_solver_py failed: No error produced with an unknown "
                "option."
            )
            test_passed = False
        except RuntimeError:
            pass
        if not messages:
            print(
                " test_s_gek_solver_py failed: Error message not logged with an "
                "unknown option."
            )
            test_passed = False

        self.assertTrue(test_passed, "test_s_gek_solver_py failed")
        print(" test_s_gek_solver_py PASSED")

    def test_s_gek_settings_layout(self):
        """
        this function checks that every field of the Python S-GEK settings structure is
        read back under its own name from settings whose fields Fortran sets one by one
        to the reference values
        """
        test_passed = check_settings_layout(
            "test_s_gek_settings_layout",
            SGEKSettings,
            lib.get_reference_s_gek_values,
            lib.reference_s_gek_field,
        )

        self.assertTrue(test_passed, "test_s_gek_settings_layout failed")
        print(" test_s_gek_settings_layout PASSED")
