# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import unittest
from ctypes import POINTER, c_bool, c_void_p, pointer
from unittest.mock import patch

import numpy as np

from pyopentrustregion.tests import (
    lib,
    n_param,
    add_tests,
    print_separator,
    check_callback_interfaces,
    Hartmann6DProblem,
)
from pyopentrustregion.python_interface import (
    c_real,
    hess_x_interface_type,
    precond_interface_type,
    SolverSettings,
    SolverState,
    solver,
)
from pyopentrustregion.extensions.rotation_qn.python_interface import (
    TransportInterface,
    ObjFuncGradInterface,
    obj_func_grad_interface_type,
    transport_interface_type,
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
        "rotation_qn_factory_c",
        "rotation_qn_settings_init",
        "rotation_qn_settings_layout",
    ],
    "rotation_qn_c_interface_tests": [
        "assign_rotation_qn_c_f",
        "assign_rotation_qn_f_c",
        "init_rotation_qn_settings_c",
        "obj_func_grad_f_wrapper",
        "rotation_qn_factory_c_wrapper",
        "transport_f_wrapper",
    ],
    "rotation_qn_tests": [
        "add_pair_rotation_qn",
        "add_rotation_qn",
        "build_model_bfgs",
        "build_model_sr1",
        "hess_x_rotation_qn",
        "init_precond_rotation_qn",
        "init_rotation_qn",
        "invoke_init_hess",
        "invoke_init_precond",
        "invoke_obj_func_grad",
        "invoke_transport",
        "keep_pair",
        "model_x_rotation_qn",
        "obj_func_rotation_qn",
        "precond_rotation_qn",
        "rotation_qn_factory",
        "rotation_qn_sanity_check",
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

    def mock_init_precond(self, residual, mu, precond_residual):
        """
        this function is a mock function for the preconditioner of the settings
        """
        precond_residual[:] = mu * residual

    def mock_obj_func_grad(self, kappa, grad):
        """
        this function is a mock function for the objective function of the settings,
        which returns the sum of the step as the energy and twice the step as the
        gradient
        """
        grad[:] = 2 * kappa

        return np.sum(kappa)

    def mock_transport(self, geodesic, tangent_vector):
        """
        this function is a mock function for the transport function
        """
        tangent_vector[:] = geodesic * tangent_vector

        return

    def test_callback_interfaces(self):
        """
        this function tests that the interfaces to the objective function and the
        transport function of the settings report an exception or an interruption
        raised by the function as an error and keep it with the call whose library
        called them, and that they receive the step and the geodesic read-only
        """
        vector = np.zeros(n_param)
        vector_ptr = vector.ctypes.data_as(POINTER(c_real))
        output = np.zeros(n_param)
        output_ptr = output.ctypes.data_as(POINTER(c_real))
        func_ptr = pointer(c_real())
        interfaces = {
            "objective function of the settings": lambda func: ObjFuncGradInterface(
                {"obj_func_grad": func}, "obj_func_grad"
            )(vector_ptr, output_ptr, func_ptr, None),
            "transport function of the settings": lambda func: TransportInterface(
                {"transport": func}, "transport"
            )(vector_ptr, output_ptr, None),
        }
        test_passed = check_callback_interfaces(
            "test_callback_interfaces",
            interfaces,
            [
                "objective function of the settings",
                "transport function of the settings",
            ],
            vector,
            n_param,
        )

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
        settings.init_hess = self.mock_hess_x
        settings.init_precond = self.mock_init_precond
        settings.obj_func_grad = self.mock_obj_func_grad
        settings.transport = self.mock_transport
        solver_settings, python_functions = setup_wired_solver_settings()

        # call rotation QN factory python interface
        mock_passed.value = False
        update_orbs_rotation_qn, obj_func_rotation_qn = rotation_qn_factory(
            n_param,
            self.mock_update_orbs,
            self.mock_obj_func,
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

        self.assertTrue(test_passed, "test_rotation_qn_factory_py_interface failed")
        print(" test_rotation_qn_factory_py_interface PASSED")

    @patch.object(
        RotationQNSettings, "init_c_struct", lib.mock_init_rotation_qn_settings
    )
    def test_rotation_qn_settings(self):
        """
        this function ensures the RotationQNSettings object is properly initialized and
        synchronized with the underlying C struct, and that the C functions of every
        callback field in the C settings are read as callables around them
        """
        test_passed = settings_hold_reference(
            "test_rotation_qn_settings",
            RotationQNSettings(),
            lib.reference_rotation_qn_field,
        )

        # a C objective function in the C settings is read as a callable around it,
        # which hands it the step, returns its energy and gradient, calls it with the
        # context of the settings and raises its error
        settings = RotationQNSettings()
        settings.settings_c.context = 6789
        received_contexts = []
        error_code = 0

        def c_obj_func_grad(kappa_ptr, grad_ptr, func_ptr, context):
            kappa = np.ctypeslib.as_array(kappa_ptr, shape=(n_param,))
            np.ctypeslib.as_array(grad_ptr, shape=(n_param,))[:] = 2.0 * kappa
            func_ptr[0] = np.sum(kappa)
            received_contexts.append(context)
            return error_code

        c_function = obj_func_grad_interface_type(c_obj_func_grad)
        settings.settings_c.obj_func_grad = c_void_p.from_buffer(c_function).value
        kappa = np.random.random(n_param)
        # a strided gradient, which the library cannot write in place, has to receive
        # the values written into its copy
        grad = np.zeros(2 * n_param)[::2]
        func = settings.obj_func_grad(kappa, grad)
        if (
            not np.isclose(func, np.sum(kappa))
            or not np.allclose(grad, 2.0 * kappa)
            or received_contexts != [6789]
        ):
            print(
                " test_rotation_qn_settings failed: C objective function not called "
                "with the step and the context of the settings or its energy and "
                "gradient not returned."
            )
            test_passed = False
        error_code = 3
        try:
            settings.obj_func_grad(kappa, grad)
            print(
                " test_rotation_qn_settings failed: Error of the C objective function "
                "not raised."
            )
            test_passed = False
        except RuntimeError:
            pass

        # the C initial Hessian, its preconditioner and transport function in the C
        # settings are read as callables around them too, which hand them their inputs
        # and the context of the settings, write their outputs back, also into an array
        # the library cannot write in place, and raise their errors
        def c_init_hess(x_ptr, hess_x_ptr, context):
            x = np.ctypeslib.as_array(x_ptr, shape=(n_param,))
            np.ctypeslib.as_array(hess_x_ptr, shape=(n_param,))[:] = 3.0 * x
            received_contexts.append(context)
            return error_code

        def c_init_precond(residual_ptr, mu, precond_residual_ptr, context):
            residual = np.ctypeslib.as_array(residual_ptr, shape=(n_param,))
            precond_residual = np.ctypeslib.as_array(
                precond_residual_ptr, shape=(n_param,)
            )
            precond_residual[:] = residual / (3.0 - mu[0])
            received_contexts.append(context)
            return error_code

        def c_transport(geodesic_ptr, tangent_vector_ptr, context):
            geodesic = np.ctypeslib.as_array(geodesic_ptr, shape=(n_param,))
            np.ctypeslib.as_array(tangent_vector_ptr, shape=(n_param,))[:] += geodesic
            received_contexts.append(context)
            return error_code

        x, y = np.random.random(n_param), np.random.random(n_param)
        for name, c_function, call, expected in (
            (
                "init_hess",
                hess_x_interface_type(c_init_hess),
                lambda function, output: function(x, output),
                3.0 * x,
            ),
            (
                "init_precond",
                precond_interface_type(c_init_precond),
                lambda function, output: function(x, 0.5, output),
                x / 2.5,
            ),
            (
                "transport",
                transport_interface_type(c_transport),
                lambda function, output: function(y, output),
                x + y,
            ),
        ):
            settings = RotationQNSettings()
            settings.settings_c.context = 6789
            setattr(settings.settings_c, name, c_void_p.from_buffer(c_function).value)
            received_contexts = []
            error_code = 0
            output = np.zeros(2 * n_param)[::2]
            output[:] = x
            call(getattr(settings, name), output)
            if not np.allclose(output, expected) or received_contexts != [6789]:
                print(
                    f" test_rotation_qn_settings failed: C function of {name} not "
                    "called with the inputs and the context of the settings or its "
                    "output not returned."
                )
                test_passed = False
            error_code = 3
            try:
                call(getattr(settings, name), output)
                print(
                    f" test_rotation_qn_settings failed: Error of the C function of "
                    f"{name} not raised."
                )
                test_passed = False
            except RuntimeError:
                pass

        self.assertTrue(test_passed, "test_rotation_qn_settings failed")
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


class RotationQNPyIntegrationTests(Hartmann6DProblem, unittest.TestCase):
    """
    this class contains integration tests that drive the real library (not the mock
    library) through the rotation QN Python interface
    """

    @classmethod
    def setUpClass(cls):
        print_separator("Running integration tests for rotation QN Python interface...")
        cls.load_hartmann6d()
        return super().setUpClass()

    def test_rotation_qn_factory_py(self):
        """
        this function solves the Hartmann 6D problem with the real factory and solver
        through the Python interface, with either updating scheme, SR1 with every
        callback function of the settings, the Hessian diagonal as the initial Hessian,
        its level-shifted inverse as its preconditioner, the objective function of the
        settings, which has to replace the host's, and an identity as the transport
        function, and BFGS with none, whose defaults are the same apart from the
        objective function, so that the host's is called, and checks that the solve
        converges to the minimum, that every callback function is called, while the
        exact Hessian the orbital update returns is not, and that the factory lets the
        solver refresh the model after a rejected step only with the objective
        function of the settings
        """
        test_passed = True
        for qn_method in ("sr1", "bfgs"):
            callbacks = qn_method == "sr1"
            curr_vars = self.near_minimum.copy()
            curr_hess = np.empty((self.n_param, self.n_param))
            counts = dict.fromkeys(
                (
                    "update_orbs",
                    "obj_func",
                    "exact_hess_x",
                    "init_hess",
                    "init_precond",
                    "obj_func_grad",
                    "transport",
                ),
                0,
            )

            def exact_hess_x(x, hess_x):
                counts["exact_hess_x"] += 1
                hess_x[:] = curr_hess @ x

            def init_hess(x, hess_x):
                counts["init_hess"] += 1
                hess_x[:] = np.diag(curr_hess) * x

            def init_precond(residual, mu, precond_residual):
                counts["init_precond"] += 1
                divisors = np.diag(curr_hess) - mu
                divisors[np.abs(divisors) < 1e-10] = 1e-10
                precond_residual[:] = residual / divisors

            def update_orbs(kappa, grad, h_diag):
                counts["update_orbs"] += 1
                curr_vars[:] += kappa
                grad[:] = self._grad(curr_vars)
                curr_hess[:] = self._hess(curr_vars)
                h_diag[:] = np.diag(curr_hess)
                return self._func(curr_vars), exact_hess_x

            def obj_func(kappa):
                counts["obj_func"] += 1
                return self._func(curr_vars + kappa)

            def obj_func_grad(kappa, grad):
                counts["obj_func_grad"] += 1
                grad[:] = self._grad(curr_vars + kappa)
                return self._func(curr_vars + kappa)

            def transport(geodesic, tangent_vector):
                counts["transport"] += 1

            settings = SolverSettings()
            state = SolverState()
            rotation_qn_settings = RotationQNSettings()
            rotation_qn_settings.qn_method = qn_method
            if callbacks:
                rotation_qn_settings.init_hess = init_hess
                rotation_qn_settings.init_precond = init_precond
                rotation_qn_settings.obj_func_grad = obj_func_grad
                rotation_qn_settings.transport = transport
            update_orbs_rotation_qn, obj_func_rotation_qn = rotation_qn_factory(
                self.n_param,
                update_orbs,
                obj_func,
                settings,
                state,
                rotation_qn_settings,
            )
            if settings.refresh_hess != callbacks:
                print(
                    " test_rotation_qn_factory_py failed: Refresh of the model after a "
                    "rejected step not set with exactly the objective function of the "
                    f"settings for {qn_method}."
                )
                test_passed = False

            solver(
                update_orbs_rotation_qn,
                obj_func_rotation_qn,
                self.n_param,
                settings,
                state,
            )
            if not np.allclose(curr_vars, self.minimum1, rtol=0.0, atol=1e-4):
                print(
                    " test_rotation_qn_factory_py failed: Solve did not converge to "
                    f"the minimum for {qn_method}."
                )
                test_passed = False
            if (
                counts["update_orbs"] == 0
                or counts["exact_hess_x"] != 0
                or (counts["obj_func"] > 0) == callbacks
                or any(
                    (counts[name] > 0) != callbacks
                    for name in (
                        "init_hess",
                        "init_precond",
                        "obj_func_grad",
                        "transport",
                    )
                )
            ):
                print(
                    " test_rotation_qn_factory_py failed: Callback functions not "
                    "called, host's objective function called together with the one "
                    f"of the settings or exact Hessian called for {qn_method}."
                )
                test_passed = False

        self.assertTrue(test_passed, "test_rotation_qn_factory_py failed")
        print(" test_rotation_qn_factory_py PASSED")

    def test_rotation_qn_settings_layout(self):
        """
        this function checks that every field of the Python rotation QN settings
        structure is read back under its own name from settings whose fields Fortran
        sets one by one to the reference values
        """
        self.assertTrue(
            check_settings_layout(
                "test_rotation_qn_settings_layout",
                RotationQNSettings,
                lib.get_reference_rotation_qn_values,
                lib.reference_rotation_qn_field,
            ),
            "test_rotation_qn_settings_layout failed",
        )
        print(" test_rotation_qn_settings_layout PASSED")
