# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import unittest
from ctypes import c_bool, c_size_t
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
from pyopentrustregion.extensions.common.python_interface import OrbitalOwner
from pyopentrustregion.extensions.common.tests import (
    OrbitalCallbackMocks,
    ToyRun,
    check_released,
    evaluate_dm_missing_warning,
    toy_model,
    check_settings_layout,
    n_ao,
    setup_wired_solver_settings,
    check_wired_solver_settings,
    check_extension_owner,
    release_mock_extension,
    assign_reference_settings,
    settings_hold_reference,
    ao_overlap_pattern,
)
from pyopentrustregion.extensions.oao import OAOSettings, oao_factory

# define all tests in alphabetical order
fortran_tests = {
    "oao_c_integration_tests": [
        "oao_factory_c",
        "oao_settings_init",
        "oao_settings_layout",
    ],
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


def density_pattern(n_particle, offset=0.0):
    """
    this function returns density matrices of the given number of particle channels
    whose values encode their channel and the pair of their AO indices regardless of
    their order together with an offset, as density_pattern of the Fortran test
    reference does, so that they are symmetric while their values reveal matrices of
    the wrong shape or channel
    """
    return offset + np.stack(
        [100 * k + ao_overlap_pattern() for k in range(1, n_particle + 1)]
    )


def setup_oao_toy_run(n_particle):
    """
    this function prepares a solve of the toy model in the OAO basis and the given shell
    through the OAO factory and returns the run and the functions the factory returned,
    also for the integration tests of ARH
    """
    run = ToyRun(False, n_particle)
    update_orbs, obj_func = oao_factory(
        run.orbitals,
        toy_model().s,
        n_particle,
        n_ao,
        run.evaluate_dm,
        run.settings,
        run.state,
    )
    return run, update_orbs, obj_func


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

    mock_evaluate_dm = OrbitalCallbackMocks.mock_evaluate_dm
    mock_get_response = OrbitalCallbackMocks.mock_get_response

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
        ao_overlap = ao_overlap_pattern()

        # initialize test flag
        test_passed = True

        # the logging function records the messages, and the mock factory records in a
        # flag whether it received what it expects, which is cleared before and read
        # after every call of it
        messages = []
        mock_passed = c_bool.in_dll(lib, "test_oao_factory_interface")

        # initialize settings object
        settings = OAOSettings()
        assign_reference_settings(settings, lib.reference_oao_field)
        settings.logger = messages.append

        # initialize density matrix
        dm_ao = density_pattern(1)[0]

        # initialize solver settings object
        solver_settings, python_functions = setup_wired_solver_settings()

        # call OAO factory python interface
        mock_passed.value = False
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

        # check if the mock factory received the correct input
        if not mock_passed.value:
            print(
                " test_oao_factory_py_interface failed: Mock factory received wrong "
                "input for settings."
            )
            test_passed = False

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
            "test_oao_factory_py_interface",
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
        if "test" not in messages:
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
        dm_ao = density_pattern(n_particle)

        # call OAO factory python interface with an AO overlap matrix of single
        # precision, which the factory has to convert
        mock_passed.value = False
        oao_factory(
            dm_ao,
            ao_overlap.astype(np.float32),
            n_particle,
            n_ao,
            self.mock_evaluate_dm,
            setup_wired_solver_settings()[0],
            SolverState(),
            settings,
        )

        # check if the mock factory received the correct input
        if not mock_passed.value:
            print(
                " test_oao_factory_py_interface failed: Mock factory received wrong "
                "input for an AO overlap matrix of single precision."
            )
            test_passed = False

        # call OAO factory python interface without settings and without a density
        # matrix evaluating function, which have to be passed as none, and determine if
        # the owner of the object holds the dimensions of the orbitals for the ARH
        # factory and keeps the density matrix alive
        settings_flag = c_bool.in_dll(lib, "oao_settings_passed")
        evaluate_dm_flag = c_bool.in_dll(lib, "oao_evaluate_dm_passed")
        settings_flag.value = True
        evaluate_dm_flag.value = True
        mock_passed.value = False
        dm_ao = density_pattern(n_particle)
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

        # check if the mock factory received the correct input
        if not mock_passed.value:
            print(
                " test_oao_factory_py_interface failed: Mock factory received wrong "
                "input without settings."
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
        mock_passed.value = False
        if not check_extension_owner(
            lambda solver_settings, solver_state: oao_factory(
                density_pattern(1)[0],
                ao_overlap,
                1,
                n_ao,
                self.mock_evaluate_dm,
                solver_settings,
                solver_state,
                settings,
            ),
            "test_oao_factory_py_interface",
        ):
            test_passed = False

        # check if the mock factory received the correct input
        if not mock_passed.value:
            print(
                " test_oao_factory_py_interface failed: Mock factory received wrong "
                "input for the check of the owner."
            )
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

        self.assertTrue(test_passed, "test_oao_factory_py_interface failed")
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
class OAOCIntegrationTests(unittest.TestCase):
    """
    this class contains integration tests that drive the real library through the
    OAO C header
    """

    tests = fortran_tests["oao_c_integration_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running integration tests for OAO C interface...")
        return super().setUpClass()


class OAOPyIntegrationTests(unittest.TestCase):
    """
    this class contains integration tests that drive the real library (not the mock
    library) through the OAO Python interface
    """

    @classmethod
    def setUpClass(cls):
        print_separator("Running integration tests for OAO Python interface...")
        return super().setUpClass()

    def test_oao_factory_py(self):
        """
        this function solves the toy model through the OAO factory in the closed and
        the open shell, checks that the OAO object is released once nothing keeps it
        alive and that the factory warns through the logging function of the OAO
        settings about a missing density matrix evaluating function
        """
        test_passed = True

        # solve the toy model in the closed and the open shell
        for n_particle, case_name in ((1, "closed-shell"), (2, "open-shell")):
            run, update_orbs, obj_func = setup_oao_toy_run(n_particle)
            releases = [run.settings._extension_owner.release]
            test_passed &= run.solve(
                update_orbs, obj_func, "test_oao_factory_py", case_name
            )
            del run, update_orbs, obj_func
            test_passed &= check_released(releases, "test_oao_factory_py", case_name)

        # check that the factory warns about a missing density matrix evaluating
        # function
        run = ToyRun(False, 1)
        messages = []
        settings = OAOSettings()
        settings.logger = messages.append
        settings.verbose = 3
        oao_factory(
            run.orbitals,
            toy_model().s,
            1,
            n_ao,
            None,
            run.settings,
            run.state,
            settings,
        )
        if evaluate_dm_missing_warning not in "".join(messages):
            print(
                " test_oao_factory_py failed: No warning about a missing density matrix "
                "evaluating function."
            )
            test_passed = False

        self.assertTrue(test_passed, "test_oao_factory_py failed")
        print(" test_oao_factory_py PASSED")

    def test_oao_settings_layout(self):
        """
        this function checks that the Python OAO settings structure has the size of
        the bind(C) settings type and that every field is read back under its own name
        from settings whose fields Fortran sets one by one to the reference values
        """
        lib.oao_settings_size.restype = c_size_t
        self.assertTrue(
            check_settings_layout(
                "test_oao_settings_layout",
                OAOSettings,
                lib.get_reference_oao_values,
                lib.reference_oao_field,
                lib.oao_settings_size(),
            ),
            "test_oao_settings_layout failed",
        )
        print(" test_oao_settings_layout PASSED")
