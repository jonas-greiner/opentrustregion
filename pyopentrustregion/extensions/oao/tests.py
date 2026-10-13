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
    logger_interface_type,
)
from pyopentrustregion.extensions.common.tests import (
    n_ao,
    setup_wired_solver_settings,
    check_wired_solver_settings,
)
from pyopentrustregion.extensions.oao import OAOSettings, oao_factory, oao_deconstructor

if NUMPY_AVAILABLE:
    import numpy as np


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
        "oao_deconstructor",
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
        "evaluate_dm_oao_f_wrapper",
        "get_response_oao_f_wrapper",
        "init_oao_settings_c",
        "oao_deconstructor_c_wrapper",
        "oao_factory_c_wrapper",
        "update_orbs_oao_c_wrapper",
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


@unittest.skipUnless(NUMPY_AVAILABLE, "NumPy not available.")
class OAOPyInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the Python interface
    """

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for OAO Python interface...")

        # get fields
        fields = OAOSettings.c_struct._fields_

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
        lib.get_reference_oao_values.argtypes = [POINTER(RefSettingsC)]
        lib.get_reference_oao_values.restype = None
        lib.get_reference_oao_values(byref(ref_settings))

        # extract Python values
        for name, _ in ref_fields:
            ref_value = getattr(ref_settings, name)
            if isinstance(ref_value, bytes):
                ref_value = ref_value.decode("utf-8")
            setattr(cls, name, ref_value)

        return super().setUpClass()

    assign_ref_to_settings = LegacyPyInterfaceHelpers.assign_ref_to_settings
    equal_settings_to_ref = LegacyPyInterfaceHelpers.equal_settings_to_ref

    mock_logger = LegacyPyInterfaceHelpers.mock_logger

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

    # replace original library with mock library
    @patch("pyopentrustregion.python_interface.lib.oao_factory", lib.mock_oao_factory)
    def test_oao_factory_py_interface(self):
        """
        this function tests the OAO factory python interface
        """
        ao_overlap = np.full(2 * (n_ao,), 2.0, dtype=np.float64)

        # initialize test flag
        test_passed = True

        # initialize settings object
        settings = OAOSettings()
        self.assign_ref_to_settings(settings)

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
            settings,
        )

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
            settings,
        )

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
                oao_factory(*args, SolverSettings(), settings)
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

    # replace original library with mock library
    @patch(
        "pyopentrustregion.python_interface.lib.oao_deconstructor",
        lib.mock_oao_deconstructor,
    )
    def test_oao_deconstructor_py_interface(self):
        """
        this function tests the OAO deconstructor python interface
        """
        # initialize test flag
        test_passed = True

        # call OAO deconstructor python interface
        oao_deconstructor()

        # check if deconstructor was called correctly
        if not c_bool.in_dll(lib, "test_oao_deconstructor_interface").value:
            print(
                " test_oao_deconstructor_py_interface failed: Deconstructor called "
                "wrong."
            )
            test_passed = False

        self.assertTrue(test_passed, "test_oao_deconstructor_py_interface failed")
        print(" test_oao_deconstructor_py_interface PASSED")

    @patch.object(OAOSettings, "init_c_struct", lib.mock_init_oao_settings)
    def test_oao_settings(self):
        """
        this function ensure the OAOSettings object is properly initialized and
        synchronized with the underlying C struct
        """
        test_passed = True
        settings = OAOSettings()
        if not self.equal_settings_to_ref(settings):
            print(" test_oao_settings failed: Settings not initialized correctly.")
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
