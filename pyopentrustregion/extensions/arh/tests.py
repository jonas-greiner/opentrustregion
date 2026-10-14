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
from pyopentrustregion.extensions.arh import ARHSettings, arh_factory
from pyopentrustregion.extensions.common.python_interface import (
    ExtensionOwner,
    OrbitalOwner,
)
from pyopentrustregion.extensions.common.tests import (
    n_ao,
    setup_wired_solver_settings,
    check_wired_solver_settings,
    check_extension_owner,
    release_mock_extension,
)

if NUMPY_AVAILABLE:
    import numpy as np


# define all tests in alphabetical order
fortran_tests = {
    "arh_c_system_tests": ["arh_settings_init"],
    "arh_tests": [
        "apply_ms_sr1_skip",
        "arh_factory_common",
        "arh_factory_cs",
        "arh_factory_os",
        "arh_sanity_check",
        "arh_set_solver_settings",
        "build_a_block_linear_os",
        "build_a_block_nonlinear_os",
        "build_a_part",
        "build_a_transformed",
        "build_hess_model_cs",
        "build_hess_model_os",
        "cache_channel_split_dirs",
        "cache_combined_channel_dirs",
        "combine_channels",
        "congruence_transform",
        "construct_arh_mo",
        "construct_arh_oao",
        "cross_symmetrize",
        "density_in_history",
        "factorize_history",
        "get_extra_trial_vectors_arh",
        "get_low_rank_hess_factors",
        "get_ms_a_inv",
        "get_ms_a_inv_jk_cs",
        "get_ms_a_inv_jk_os",
        "get_ms_jk_inv",
        "hess_x_arh",
        "hess_x_static_arh_mo",
        "hess_x_static_arh_oao",
        "history_channel_rows_arh_mo",
        "history_channel_rows_arh_oao",
        "history_columns_arh_mo",
        "history_columns_arh_oao",
        "history_dm_arh_mo",
        "history_dm_arh_oao",
        "history_step_mask",
        "init_arh_settings",
        "inv_hess_x_arh",
        "median",
        "obj_func_arh",
        "packed_channel_rows_arh_mo",
        "packed_channel_rows_arh_oao",
        "precond_arh",
        "precond_pd_arh",
        "prepend",
        "project_arh_oao",
        "rebase_dirs",
        "rebuild_stale_hess_model",
        "resolvable_residual",
        "response_gram",
        "rotate_trial_arh_mo",
        "rotate_trial_arh_oao",
        "spectral_to_dense",
        "to_history_basis_arh_mo",
        "to_history_basis_arh_oao",
        "truncated_eigval_inv",
        "update_orbs_arh",
    ],
    "arh_c_interface_tests": [
        "arh_factory_c_wrapper",
        "assign_arh_c_f",
        "assign_arh_f_c",
        "evaluate_dm_cs_f_wrapper",
        "evaluate_dm_os_f_wrapper",
        "init_arh_settings_c",
    ],
}

# multiples of the density matrix the mock density matrix evaluating function returns
# for each optional output, in the order of its argument list
arh_evaluate_dm_factors = list((c_real * 4).in_dll(lib, "test_arh_evaluate_dm_factors"))


def setup_orbital_owner(solver_settings, n_particle, copy_back=None):
    """
    this function lets solver settings keep an owner of an orbital basis object with
    the dimensions of its orbitals and the function copying them back, as an orbital
    basis factory sets one up, without the owner's own setup
    """
    owner = OrbitalOwner.__new__(OrbitalOwner)
    owner.n_ao = n_ao
    owner.n_particle = n_particle
    owner.copy_back = copy_back
    solver_settings._extension_owner = owner
    return owner


@add_tests
class ARHTests(unittest.TestCase):
    """
    this class contains unit tests for ARH
    """

    tests = fortran_tests["arh_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for ARH...")
        return super().setUpClass()


@add_tests
class ARHCInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the ARH C interface
    """

    tests = fortran_tests["arh_c_interface_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for ARH C interface...")
        return super().setUpClass()


@unittest.skipUnless(NUMPY_AVAILABLE, "NumPy not available.")
class ARHPyInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the Python interface
    """

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for ARH Python interface...")

        # get fields
        fields = ARHSettings.c_struct._fields_

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
        lib.get_reference_arh_values.argtypes = [POINTER(RefSettingsC)]
        lib.get_reference_arh_values.restype = None
        lib.get_reference_arh_values(byref(ref_settings))

        # extract Python values
        for name, _ in ref_fields:
            ref_value = getattr(ref_settings, name)
            if isinstance(ref_value, bytes):
                ref_value = ref_value.decode("utf-8")
            setattr(cls, name, ref_value)

        return super().setUpClass()

    assign_ref_to_settings = LegacyPyInterfaceHelpers.assign_ref_to_settings
    equal_settings_to_ref = LegacyPyInterfaceHelpers.equal_settings_to_ref

    def mock_arh_evaluate_dm(self, dm_ao, fock, v_coulomb, v_exchange, v_nonlinear):
        """
        this function is a mock function for the density matrix evaluating function
        with separate Coulomb, exact-exchange and non-linear potential contributions
        for both shells
        """
        if fock is not None:
            fock[:] = arh_evaluate_dm_factors[0] * dm_ao
        if v_coulomb is not None:
            v_coulomb[:] = arh_evaluate_dm_factors[1] * dm_ao
        if v_exchange is not None:
            v_exchange[:] = arh_evaluate_dm_factors[2] * dm_ao
        if v_nonlinear is not None:
            v_nonlinear[:] = arh_evaluate_dm_factors[3] * dm_ao

        return np.sum(dm_ao)

    mock_logger = LegacyPyInterfaceHelpers.mock_logger

    # replace original library with mock library, also for the release of the
    # extensions the factory returns
    @patch("pyopentrustregion.python_interface.lib.arh_factory", lib.mock_arh_factory)
    @patch(
        "pyopentrustregion.extensions.common.python_interface.release_extension",
        release_mock_extension,
    )
    def test_arh_factory_py_interface(self):
        """
        this function tests the ARH factory python interface
        """
        # initialize test flag
        test_passed = True

        # initialize settings object
        settings = ARHSettings()
        self.assign_ref_to_settings(settings)

        # call the ARH factory python interface for solver settings an orbital basis
        # factory wired last for either shell, whose owner holds the dimensions of the
        # orbitals and the function copying them back
        n_particle_flag = c_int.in_dll(lib, "arh_n_particle")
        for n_particle in (1, 2):
            n_particle_flag.value = n_particle
            solver_settings, python_functions = setup_wired_solver_settings()
            copy_back_calls = []
            orbitals = setup_orbital_owner(
                solver_settings, n_particle, lambda: copy_back_calls.append(True)
            )
            self.test_logger = True
            update_orbs_arh, obj_func_arh = arh_factory(
                self.mock_arh_evaluate_dm, solver_settings, SolverState(), settings
            )

            # check that the settings were passed and the logging function called
            if not c_bool.in_dll(lib, "arh_settings_passed").value:
                print(" test_arh_factory_py_interface failed: Settings not passed.")
                test_passed = False
            if not self.test_logger:
                print(
                    " test_arh_factory_py_interface failed: Called logging function "
                    "wrong."
                )
                test_passed = False

            # check the solver settings the factory wired and adjusted to the
            # approximate Hessian and that the owner of the ARH object keeps the owner
            # of the orbital basis object
            if not check_wired_solver_settings(
                solver_settings,
                python_functions,
                "arh_factory_py_interface",
                (update_orbs_arh, obj_func_arh),
                settings,
            ):
                test_passed = False
            if not (
                solver_settings.refresh_hess
                and solver_settings.n_micro == 300
                and not solver_settings.hess_symm
            ):
                print(
                    " test_arh_factory_py_interface failed: Solver settings not "
                    "adjusted to the approximate Hessian."
                )
                test_passed = False
            if (
                solver_settings._extension_owner.saved_objects["wired_before"]
                is not orbitals
            ):
                print(
                    " test_arh_factory_py_interface failed: Owner of the orbital basis "
                    "object not kept."
                )
                test_passed = False

            # call returned ARH objective function
            kappa = np.ones(n_param, dtype=np.float64)
            if obj_func_arh(kappa) != 3.0:
                print(
                    " test_arh_factory_py_interface failed: Returned function value of "
                    "returned ARH objective function wrong."
                )
                test_passed = False

            # call returned ARH orbital updating function, which has to copy the
            # orbitals back like that of the orbital basis object
            grad = np.empty(n_param, dtype=np.float64)
            h_diag = np.empty(n_param, dtype=np.float64)
            func, hess_x_funptr = update_orbs_arh(kappa, grad, h_diag)
            if (
                func != 3.0
                or not np.allclose(grad, 2.0)
                or not np.allclose(h_diag, 3.0)
            ):
                print(
                    " test_arh_factory_py_interface failed: Returned values of "
                    "returned ARH orbital updating function wrong."
                )
                test_passed = False
            if copy_back_calls != [True]:
                print(
                    " test_arh_factory_py_interface failed: Orbitals not copied back "
                    "by the returned ARH orbital updating function."
                )
                test_passed = False

            # call returned hess_x function
            hess_x = np.empty(n_param, dtype=np.float64)
            hess_x_funptr(np.ones(n_param, dtype=np.float64), hess_x)
            if not np.allclose(hess_x, 4.0):
                print(
                    " test_arh_factory_py_interface failed: Returned Hessian linear "
                    "transformation of returned ARH orbital updating function wrong."
                )
                test_passed = False

        # call the ARH factory python interface without settings, which have to be
        # passed as none so that the library uses its default settings
        settings_flag = c_bool.in_dll(lib, "arh_settings_passed")
        settings_flag.value = True
        solver_settings = setup_wired_solver_settings()[0]
        setup_orbital_owner(solver_settings, 2)
        arh_factory(self.mock_arh_evaluate_dm, solver_settings, SolverState())
        if settings_flag.value:
            print(
                " test_arh_factory_py_interface failed: Settings passed although none "
                "were given."
            )
            test_passed = False

        # the factory hands the extension it receives to an owner, which the solver
        # settings keep alive and which keeps the state alive
        def call_factory(solver_settings, solver_state):
            setup_orbital_owner(solver_settings, 2)
            arh_factory(
                self.mock_arh_evaluate_dm, solver_settings, solver_state, settings
            )

        if not check_extension_owner(call_factory, "arh_factory_py_interface"):
            test_passed = False

        # call the ARH factory python interface for solver settings no orbital basis
        # factory wired last and for an evaluate_dm that cannot be called with the
        # arguments the C interface passes, which have to raise an error before the
        # library is called
        other_settings = SolverSettings()
        other_settings._extension_owner = ExtensionOwner.__new__(ExtensionOwner)
        orbital_settings = SolverSettings()
        setup_orbital_owner(orbital_settings, 2)
        invalid_cases = (
            (
                "solver settings no factory wired",
                SolverSettings(),
                self.mock_arh_evaluate_dm,
                ValueError,
            ),
            (
                "solver settings another factory wired last",
                other_settings,
                self.mock_arh_evaluate_dm,
                ValueError,
            ),
            (
                "evaluate_dm with a wrong number of arguments",
                orbital_settings,
                lambda dm, fock: 0.0,
                TypeError,
            ),
        )
        for case, case_settings, evaluate_dm, error_type in invalid_cases:
            try:
                arh_factory(evaluate_dm, case_settings, SolverState(), settings)
                print(
                    f" test_arh_factory_py_interface failed: {case} did not raise an "
                    "error."
                )
                test_passed = False
            except error_type:
                pass

        self.assertTrue(
            c_bool.in_dll(lib, "test_arh_factory_interface").value and test_passed,
            "test_arh_factory_py_interface failed",
        )
        print(" test_arh_factory_py_interface PASSED")

    @patch.object(ARHSettings, "init_c_struct", lib.mock_init_arh_settings)
    def test_arh_settings(self):
        """
        this function ensure the ARHSettings object is properly initialized and
        synchronized with the underlying C struct
        """
        test_passed = True
        settings = ARHSettings()
        if not self.equal_settings_to_ref(settings):
            print(" test_arh_settings failed: Settings not initialized correctly.")
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
            print(" test_arh_settings failed: Logging function not set correctly.")
            test_passed = False

        self.assertTrue(test_passed, "test_arh_settings failed")
        print(" test_arh_settings PASSED")


@add_tests
class ARHCSystemTests(unittest.TestCase):
    """
    this class contains system tests for the ARH C interface
    """

    tests = fortran_tests["arh_c_system_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running system tests for ARH C interface...")
        return super().setUpClass()
