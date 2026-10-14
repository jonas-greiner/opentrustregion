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
    c_real,
    SolverSettings,
    SolverState,
    logger_interface_type,
)
from pyopentrustregion.extensions.common.python_interface import OrbitalOwner
from pyopentrustregion.extensions.common.tests import (
    n_ao,
    setup_wired_solver_settings,
    check_wired_solver_settings,
    check_extension_owner,
    release_mock_extension,
    assign_reference_settings,
    settings_hold_reference,
    ExtensionCallbackMocks,
)
from pyopentrustregion.extensions.oao import OAOSettings, oao_factory

# define all tests in alphabetical order
fortran_tests = {
    "oao_c_system_tests": ["oao_settings_init"],
    "oao_tests": [
        "calculate_grad_h_diag_oao",
        "get_extra_trial_vectors_oao",
        "get_hess_eigval_pairs_oao",
        "hess_x_oao",
        "hess_x_static_oao",
        "init_oao_settings",
        "oao_factory_common",
        "oao_factory_cs",
        "oao_factory_os",
        "oao_sanity_check",
        "obj_func_oao",
        "pack_asymm",
        "project_asymm",
        "project_oao",
        "project_symm",
        "purify",
        "refresh_hess_eigen_oao",
        "rotate_dm_ao",
        "rotate_from_hess_eigenbasis_oao",
        "rotate_orbitals_oao",
        "rotate_to_hess_eigenbasis_oao",
        "symmetric_transformation",
        "unpack_asymm",
        "update_orbs_oao",
    ],
    "oao_c_interface_tests": [
        "assign_oao_c_f",
        "assign_oao_f_c",
        "init_oao_settings_c",
        "oao_factory_c_wrapper",
    ],
}

# multiples of the density matrix the mock density matrix evaluating function returns
# for the Fock matrix and the response
evaluate_dm_factors = list((c_real * 2).in_dll(lib, "test_evaluate_dm_factors"))


@add_tests
class OAOTests(unittest.TestCase):
    """
    this class contains unit tests for OAO
    """

    tests = fortran_tests["oao_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for OAO...")
        return super().setUpClass()


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


class OAOPyInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the Python interface
    """

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for OAO Python interface...")

        return super().setUpClass()

    mock_logger = ExtensionCallbackMocks.mock_logger

    def mock_evaluate_dm(self, dm_ao, fock, get_response):
        """
        this function is a mock function for the density matrix evaluating function,
        which builds only what it is asked for
        """
        if fock is not None:
            fock[:] = evaluate_dm_factors[0] * dm_ao

        return np.sum(dm_ao), (self.mock_get_response if get_response else None)

    def mock_get_response(self, dm_ao, response):
        """
        this function is a mock function for the response function
        """
        response[:] = evaluate_dm_factors[1] * dm_ao

        return

    # replace original library with mock library, also for the release of the
    # extensions the factory returns
    @patch("pyopentrustregion.python_interface.lib.oao_factory", lib.mock_oao_factory)
    @patch(
        "pyopentrustregion.extensions.common.python_interface.release_extension",
        release_mock_extension,
    )
    def test_oao_factory_py_interface(self):
        """
        this function tests the OAO factory python interface
        """
        ao_overlap = np.full(2 * (n_ao,), 2.0, dtype=np.float64)

        # initialize test flag
        test_passed = True

        # initialize settings object
        settings = OAOSettings()
        assign_reference_settings(settings, lib.reference_oao_field)
        settings.logger = self.mock_logger

        # initialize logging boolean
        self.test_logger = True

        # initialize density matrix
        dm_ao = np.full(2 * (n_ao,), 1.0, dtype=np.float64)

        # initialize solver settings object
        solver_settings, python_functions = setup_wired_solver_settings()

        # call OAO factory python interface
        update_orbs_oao, obj_func_oao = oao_factory(
            dm_ao,
            ao_overlap,
            1,
            n_ao,
            self.mock_evaluate_dm,
            solver_settings,
            SolverState(),
            settings,
        )

        # check that the settings and the density matrix evaluating function were
        # passed
        if not (
            c_bool.in_dll(lib, "oao_settings_passed").value
            and c_bool.in_dll(lib, "oao_evaluate_dm_passed").value
        ):
            print(
                " test_oao_factory_py_interface failed: Settings or density matrix "
                "evaluating function not passed."
            )
            test_passed = False

        # check the solver settings the factory wired
        if not check_wired_solver_settings(
            solver_settings,
            python_functions,
            "oao_factory_py_interface",
            (update_orbs_oao, obj_func_oao),
            settings,
        ):
            test_passed = False

        # check that the factory received the density matrix in place
        if not np.allclose(dm_ao, np.full(2 * (n_ao,), 2.0, dtype=np.float64)):
            print(
                " test_oao_factory_py_interface failed: Density matrix not passed in "
                "place."
            )
            test_passed = False

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
            self.mock_evaluate_dm,
            setup_wired_solver_settings()[0],
            SolverState(),
            settings,
        )

        # call OAO factory python interface without settings and without a density
        # matrix evaluating function, which have to be passed as none, and determine if
        # the owner of the object holds the dimensions of the orbitals for the ARH
        # factory and keeps the density matrix alive
        settings_flag = c_bool.in_dll(lib, "oao_settings_passed")
        evaluate_dm_flag = c_bool.in_dll(lib, "oao_evaluate_dm_passed")
        settings_flag.value = True
        evaluate_dm_flag.value = True
        dm_ao = np.full((n_particle, n_ao, n_ao), 1.0, dtype=np.float64)
        solver_settings = setup_wired_solver_settings()[0]
        try:
            oao_factory(
                dm_ao,
                ao_overlap,
                n_particle,
                n_ao,
                None,
                solver_settings,
                SolverState(),
            )
        except RuntimeError:
            print(
                " test_oao_factory_py_interface failed: Call without settings and "
                "density matrix evaluating function raises error."
            )
            test_passed = False
        if settings_flag.value or evaluate_dm_flag.value:
            print(
                " test_oao_factory_py_interface failed: Settings or density matrix "
                "evaluating function passed although none were given."
            )
            test_passed = False
        owner = getattr(solver_settings, "_extension_owner", None)
        if (
            not isinstance(owner, OrbitalOwner)
            or owner.n_ao != n_ao
            or owner.n_particle != n_particle
            or owner.saved_objects.get("dm_ao") is not dm_ao
        ):
            print(
                " test_oao_factory_py_interface failed: Owner of the object does not "
                "hold the dimensions of the orbitals or the density matrix."
            )
            test_passed = False

        # the factory hands the extension it receives to an owner, which the solver
        # settings keep alive and which keeps the state alive
        if not check_extension_owner(
            lambda solver_settings, solver_state: oao_factory(
                np.full(2 * (n_ao,), 1.0, dtype=np.float64),
                ao_overlap,
                1,
                n_ao,
                self.mock_evaluate_dm,
                solver_settings,
                solver_state,
                settings,
            ),
            "oao_factory_py_interface",
        ):
            test_passed = False

        # call oao_factory python interface for invalid density and overlap matrices,
        # which have to raise an error before the Fortran factory is called
        read_only_dm_ao = np.full(2 * (n_ao,), 1.0, dtype=np.float64)
        read_only_dm_ao.flags.writeable = False
        invalid_cases = (
            (
                "density matrix not matching the number of particles",
                (np.full(2 * (n_ao,), 1.0), ao_overlap, 2, n_ao, self.mock_evaluate_dm),
                ValueError,
            ),
            (
                "density matrix not matching the number of AOs",
                (
                    np.full(2 * (n_ao,), 1.0),
                    np.full(2 * (n_ao + 1,), 2.0),
                    1,
                    n_ao + 1,
                    self.mock_evaluate_dm,
                ),
                ValueError,
            ),
            (
                "non-contiguous density matrix",
                (
                    np.full((n_ao, 2 * n_ao), 1.0)[:, ::2],
                    ao_overlap,
                    1,
                    n_ao,
                    self.mock_evaluate_dm,
                ),
                ValueError,
            ),
            (
                "single-precision density matrix",
                (
                    np.full(2 * (n_ao,), 1.0, dtype=np.float32),
                    ao_overlap,
                    1,
                    n_ao,
                    self.mock_evaluate_dm,
                ),
                ValueError,
            ),
            (
                "read-only density matrix",
                (read_only_dm_ao, ao_overlap, 1, n_ao, self.mock_evaluate_dm),
                ValueError,
            ),
            (
                "wrongly shaped AO overlap matrix",
                (
                    np.full(2 * (n_ao,), 1.0),
                    ao_overlap[:-1, :-1],
                    1,
                    n_ao,
                    self.mock_evaluate_dm,
                ),
                ValueError,
            ),
            (
                "evaluate_dm with a wrong number of arguments",
                (
                    np.full(2 * (n_ao,), 1.0),
                    ao_overlap,
                    1,
                    n_ao,
                    lambda dm, fock: (0.0, None),
                ),
                TypeError,
            ),
        )
        for case, args, error_type in invalid_cases:
            try:
                oao_factory(*args, SolverSettings(), SolverState(), settings)
                print(
                    f" test_oao_factory_py_interface failed: {case} did not raise an "
                    "error."
                )
                test_passed = False
            except error_type:
                pass

        self.assertTrue(
            c_bool.in_dll(lib, "test_oao_factory_interface").value and test_passed,
            "test_oao_factory_py_interface failed",
        )
        print(" test_oao_factory_py_interface PASSED")

    @patch.object(OAOSettings, "init_c_struct", lib.mock_init_oao_settings)
    def test_oao_settings(self):
        """
        this function ensure the OAOSettings object is properly initialized and
        synchronized with the underlying C struct
        """
        test_passed = True
        settings = OAOSettings()
        if not settings_hold_reference(
            "test_oao_settings", settings, lib.reference_oao_field
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
            print(" test_oao_settings failed: Logging function not set correctly.")
            test_passed = False

        self.assertTrue(test_passed, "test_oao_settings failed")
        print(" test_oao_settings PASSED")


@add_tests
class OAOCSystemTests(unittest.TestCase):
    """
    this class contains system tests for the OAO C interface
    """

    tests = fortran_tests["oao_c_system_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running system tests for OAO C interface...")
        return super().setUpClass()
