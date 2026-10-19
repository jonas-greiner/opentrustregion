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
    c_real,
    c_int,
    SolverSettings,
    SolverState,
)
from pyopentrustregion.extensions.density_qn import (
    DensityQNSettings,
    density_qn_factory,
)
from pyopentrustregion.extensions.mo.tests import setup_mo_toy_run
from pyopentrustregion.extensions.oao.tests import setup_oao_toy_run
from pyopentrustregion.extensions.common.python_interface import (
    ExtensionOwner,
    OrbitalOwner,
)
from pyopentrustregion.extensions.common.tests import (
    check_released,
    check_settings_layout,
    n_ao,
    setup_wired_solver_settings,
    check_wired_solver_settings,
    check_extension_owner,
    release_mock_extension,
    assign_reference_settings,
    settings_hold_reference,
)

# define all tests in alphabetical order
fortran_tests = {
    "density_qn_c_integration_tests": [
        "density_qn_factory_c",
        "density_qn_settings_init",
        "density_qn_settings_layout",
    ],
    "density_qn_tests": [
        "apply_ms_sr1_skip",
        "density_qn_factory_common",
        "density_qn_factory_cs",
        "density_qn_factory_os",
        "density_qn_sanity_check",
        "density_qn_set_solver_settings",
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
        "construct_density_qn_mo",
        "construct_density_qn_oao",
        "cross_symmetrize",
        "density_in_history",
        "factorize_history",
        "get_extra_trial_vectors_density_qn",
        "get_low_rank_hess_factors",
        "get_ms_a_inv",
        "get_ms_a_inv_jk_cs",
        "get_ms_a_inv_jk_os",
        "get_ms_jk_inv",
        "hess_x_density_qn",
        "hess_x_static_density_qn_mo",
        "hess_x_static_density_qn_oao",
        "history_channel_rows_density_qn_mo",
        "history_channel_rows_density_qn_oao",
        "history_columns_density_qn_mo",
        "history_columns_density_qn_oao",
        "history_dm_density_qn_mo",
        "history_dm_density_qn_oao",
        "history_step_mask",
        "history_window",
        "inv_hess_x_density_qn",
        "obj_func_density_qn",
        "packed_channel_rows_density_qn_mo",
        "packed_channel_rows_density_qn_oao",
        "precond_density_qn",
        "precond_pd_density_qn",
        "prepend",
        "rebase_dirs",
        "rebuild_stale_hess_model",
        "response_gram",
        "rotate_trial_density_qn_mo",
        "rotate_trial_density_qn_oao",
        "spectral_to_dense",
        "to_history_basis_density_qn_mo",
        "to_history_basis_density_qn_oao",
        "truncated_eigval_inv",
        "update_orbs_density_qn",
    ],
    "density_qn_c_interface_tests": [
        "density_qn_evaluate_dm_cs_f_wrapper",
        "density_qn_evaluate_dm_os_f_wrapper",
        "density_qn_factory_c_wrapper",
        "assign_density_qn_c_f",
        "assign_density_qn_f_c",
        "init_density_qn_settings_c",
    ],
}

# multiples of the density matrix the mock density matrix evaluating function returns
# for each optional output, in the order of its argument list
density_qn_evaluate_dm_factors = list(
    (c_real * 4).in_dll(lib, "test_density_qn_evaluate_dm_factors")
)


def setup_orbital_owner(solver_settings, n_particle, copy_back=None):
    """
    this function lets solver settings keep an owner of an orbital basis object with
    the dimensions of its orbitals and the function copying them back, as an orbital
    basis factory sets one up, without the owner's own setup
    """
    owner = OrbitalOwner.__new__(OrbitalOwner)
    owner.n_ao = n_ao
    owner.n_particle = n_particle
    owner.exception = {}
    owner.copy_back = copy_back
    solver_settings._extension_owner = owner
    return owner


@add_tests
class DensityQNTests(unittest.TestCase):
    """
    this class contains unit tests for density QN
    """

    tests = fortran_tests["density_qn_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for density QN...")
        return super().setUpClass()


@add_tests
class DensityQNCInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the density QN C interface
    """

    tests = fortran_tests["density_qn_c_interface_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for density QN C interface...")
        return super().setUpClass()


class DensityQNPyInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the Python interface
    """

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for density QN Python interface...")

        return super().setUpClass()

    def mock_density_qn_evaluate_dm(
        self, dm_ao, fock, v_coulomb, v_exchange, v_nonlinear
    ):
        """
        this function is a mock function for the density matrix evaluating function
        with separate Coulomb, exact-exchange and non-linear potential contributions
        for both shells
        """
        if fock is not None:
            fock[:] = density_qn_evaluate_dm_factors[0] * dm_ao
        if v_coulomb is not None:
            v_coulomb[:] = density_qn_evaluate_dm_factors[1] * dm_ao
        if v_exchange is not None:
            v_exchange[:] = density_qn_evaluate_dm_factors[2] * dm_ao
        if v_nonlinear is not None:
            v_nonlinear[:] = density_qn_evaluate_dm_factors[3] * dm_ao

        return np.sum(dm_ao)

    # replace original library with mock library, also for the release of the
    # extensions the factory returns
    @patch(
        "pyopentrustregion.python_interface.lib.density_qn_factory",
        lib.mock_density_qn_factory,
    )
    @patch(
        "pyopentrustregion.extensions.common.python_interface.release_extension",
        release_mock_extension,
    )
    def test_density_qn_factory_py_interface(self):
        """
        this function tests the density QN factory python interface
        """
        # initialize test flag
        test_passed = True

        # the logging function records the messages, and the mock factory records in a
        # flag whether it received what it expects, which is cleared before and read
        # after every call of it
        messages = []
        mock_passed = c_bool.in_dll(lib, "test_density_qn_factory_interface")

        # initialize settings object
        settings = DensityQNSettings()
        assign_reference_settings(settings, lib.reference_density_qn_field)
        settings.logger = messages.append

        # call the density QN factory python interface for solver settings an orbital
        # basis factory wired last for either shell, whose owner holds the dimensions
        # of the orbitals and the function copying them back
        n_particle_flag = c_int.in_dll(lib, "density_qn_n_particle")
        for n_particle in (1, 2):
            n_particle_flag.value = n_particle
            solver_settings, python_functions = setup_wired_solver_settings()
            copy_back_calls = []
            orbitals = setup_orbital_owner(
                solver_settings, n_particle, lambda: copy_back_calls.append(True)
            )
            messages.clear()
            mock_passed.value = False
            update_orbs_density_qn, obj_func_density_qn = density_qn_factory(
                self.mock_density_qn_evaluate_dm,
                solver_settings,
                SolverState(),
                settings,
            )

            # check if the mock factory received the correct input
            if not mock_passed.value:
                print(
                    " test_density_qn_factory_py_interface failed: Mock factory "
                    "received "
                    f"wrong input for {n_particle} particle channel(s)."
                )
                test_passed = False

            # check that the settings were passed and the logging function called
            if not c_bool.in_dll(lib, "density_qn_settings_passed").value:
                print(
                    " test_density_qn_factory_py_interface failed: Settings not passed."
                )
                test_passed = False
            if "test" not in messages:
                print(
                    " test_density_qn_factory_py_interface failed: Called logging "
                    "function "
                    "wrong."
                )
                test_passed = False

            # the returned functions share the collector of the orbital basis factory,
            # whose functions the object calls
            if not (
                update_orbs_density_qn._otr_exception is orbitals.exception
                and obj_func_density_qn._otr_exception is orbitals.exception
            ):
                print(
                    " test_density_qn_factory_py_interface failed: Collector of the "
                    "orbital "
                    "basis factory not shared."
                )
                test_passed = False

            # check the solver settings the factory wired and adjusted to the
            # approximate Hessian and that the owner of the density QN object keeps the
            # owner of the orbital basis object
            if not check_wired_solver_settings(
                solver_settings,
                python_functions,
                "test_density_qn_factory_py_interface",
                (update_orbs_density_qn, obj_func_density_qn),
                settings,
            ):
                test_passed = False
            if not (
                solver_settings.refresh_hess
                and solver_settings.n_micro == 300
                and not solver_settings.hess_symm
            ):
                print(
                    " test_density_qn_factory_py_interface failed: Solver settings not "
                    "adjusted to the approximate Hessian."
                )
                test_passed = False
            if (
                solver_settings._extension_owner.saved_objects["wired_before"]
                is not orbitals
            ):
                print(
                    " test_density_qn_factory_py_interface failed: Owner of the "
                    "orbital basis "
                    "object not kept."
                )
                test_passed = False

            # call returned density QN objective function
            kappa = np.ones(n_param, dtype=np.float64)
            if obj_func_density_qn(kappa) != 3.0:
                print(
                    " test_density_qn_factory_py_interface failed: Returned function "
                    "value of "
                    "returned density QN objective function wrong."
                )
                test_passed = False

            # call returned density QN orbital updating function, which has to copy the
            # orbitals back like that of the orbital basis object
            grad = np.empty(n_param, dtype=np.float64)
            h_diag = np.empty(n_param, dtype=np.float64)
            func, hess_x_funptr = update_orbs_density_qn(kappa, grad, h_diag)
            if (
                func != 3.0
                or not np.allclose(grad, 2.0)
                or not np.allclose(h_diag, 3.0)
            ):
                print(
                    " test_density_qn_factory_py_interface failed: Returned values of "
                    "returned density QN orbital updating function wrong."
                )
                test_passed = False
            if copy_back_calls != [True]:
                print(
                    " test_density_qn_factory_py_interface failed: Orbitals not copied "
                    "back "
                    "by the returned density QN orbital updating function."
                )
                test_passed = False

            # call returned hess_x function
            hess_x = np.empty(n_param, dtype=np.float64)
            hess_x_funptr(np.ones(n_param, dtype=np.float64), hess_x)
            if not np.allclose(hess_x, 4.0):
                print(
                    " test_density_qn_factory_py_interface failed: Returned Hessian "
                    "linear "
                    "transformation of returned density QN orbital updating function "
                    "wrong."
                )
                test_passed = False

        # call the density QN factory python interface without settings, which have to
        # be passed as none so that the library uses its default settings
        settings_flag = c_bool.in_dll(lib, "density_qn_settings_passed")
        settings_flag.value = True
        solver_settings = setup_wired_solver_settings()[0]
        setup_orbital_owner(solver_settings, 2)
        mock_passed.value = False
        density_qn_factory(
            self.mock_density_qn_evaluate_dm, solver_settings, SolverState()
        )
        if settings_flag.value:
            print(
                " test_density_qn_factory_py_interface failed: Settings passed "
                "although none "
                "were given."
            )
            test_passed = False

        # check if the mock factory received the correct input
        if not mock_passed.value:
            print(
                " test_density_qn_factory_py_interface failed: Mock factory received "
                "wrong "
                "input without settings."
            )
            test_passed = False

        # the factory hands the extension it receives to an owner, which the solver
        # settings keep alive and which keeps the state alive
        def call_factory(solver_settings, solver_state):
            setup_orbital_owner(solver_settings, 2)
            density_qn_factory(
                self.mock_density_qn_evaluate_dm,
                solver_settings,
                solver_state,
                settings,
            )

        mock_passed.value = False
        if not check_extension_owner(
            call_factory, "test_density_qn_factory_py_interface"
        ):
            test_passed = False

        # check if the mock factory received the correct input
        if not mock_passed.value:
            print(
                " test_density_qn_factory_py_interface failed: Mock factory received "
                "wrong "
                "input for the check of the owner."
            )
            test_passed = False

        # call the density QN factory python interface for solver settings no orbital
        # basis factory wired last and for an evaluate_dm that cannot be called with
        # the arguments the C interface passes, which have to raise an error before the
        # library is called
        other_settings = SolverSettings()
        other_settings._extension_owner = ExtensionOwner.__new__(ExtensionOwner)
        orbital_settings = SolverSettings()
        setup_orbital_owner(orbital_settings, 2)
        invalid_cases = (
            (
                "solver settings no factory wired",
                SolverSettings(),
                self.mock_density_qn_evaluate_dm,
                ValueError,
            ),
            (
                "solver settings another factory wired last",
                other_settings,
                self.mock_density_qn_evaluate_dm,
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
                density_qn_factory(evaluate_dm, case_settings, SolverState(), settings)
                print(
                    f" test_density_qn_factory_py_interface failed: {case} did not "
                    f"raise an "
                    "error."
                )
                test_passed = False
            except error_type:
                pass

        self.assertTrue(test_passed, "test_density_qn_factory_py_interface failed")
        print(" test_density_qn_factory_py_interface PASSED")

    @patch.object(DensityQNSettings, "init_c_struct", lib.mock_init_density_qn_settings)
    def test_density_qn_settings(self):
        """
        this function ensures the DensityQNSettings object is properly initialized and
        synchronized with the underlying C struct
        """
        self.assertTrue(
            settings_hold_reference(
                "test_density_qn_settings",
                DensityQNSettings(),
                lib.reference_density_qn_field,
            ),
            "test_density_qn_settings failed",
        )
        print(" test_density_qn_settings PASSED")


@add_tests
class DensityQNCIntegrationTests(unittest.TestCase):
    """
    this class contains integration tests that drive the real library through the
    density QN C header
    """

    tests = fortran_tests["density_qn_c_integration_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running integration tests for density QN C interface...")
        return super().setUpClass()


class DensityQNPyIntegrationTests(unittest.TestCase):
    """
    this class contains integration tests that drive the real library (not the mock
    library) through the density QN Python interface
    """

    @classmethod
    def setUpClass(cls):
        print_separator("Running integration tests for density QN Python interface...")
        return super().setUpClass()

    def test_density_qn_factory_py(self):
        """
        this function solves the toy model with density QN, built through the density
        QN factory on the objects of both orbital bases in the closed and the open
        shell, compares with the minimum the exact Hessian of the orbital basis reaches
        and checks that the density QN and the orbital basis object are released once
        nothing keeps them alive
        """
        test_passed = True

        for setup_basis_run, basis_name in (
            (setup_mo_toy_run, "MO-basis"),
            (setup_oao_toy_run, "OAO-basis"),
        ):
            for n_particle, shell_name in ((1, "closed-shell"), (2, "open-shell")):
                case_name = f"{basis_name} {shell_name}"

                # solve with the exact Hessian of the orbital basis
                run, update_orbs, obj_func = setup_basis_run(n_particle)
                if not run.solve(
                    update_orbs, obj_func, "test_density_qn_factory_py", case_name
                ):
                    test_passed = False
                    continue
                exact_energy = run.state.func

                # wire the settings with a new orbital basis object and set the exact
                # Hessian linear transformation its orbital update returns at the
                # starting orbitals for the stability check, before the density QN
                # factory wires them with density QN
                run, update_orbs_basis, _ = setup_basis_run(n_particle)
                releases = [run.settings._extension_owner.release]
                _, hess_x = update_orbs_basis(
                    np.zeros(run.n_param), np.empty(run.n_param), np.empty(run.n_param)
                )
                run.settings.stability_hess_x = hess_x
                update_orbs, obj_func = density_qn_factory(
                    run.density_qn_evaluate_dm, run.settings, run.state
                )
                releases.append(run.settings._extension_owner.release)

                # solve with density QN, which reaches the same minimum
                test_passed &= run.solve(
                    update_orbs,
                    obj_func,
                    "test_density_qn_factory_py",
                    case_name,
                    density_qn=True,
                )
                if abs(run.state.func - exact_energy) > 1e-8:
                    print(
                        " test_density_qn_factory_py failed: Minimum differs from the "
                        "one of "
                        f"the exact Hessian for the {case_name} case."
                    )
                    test_passed = False
                del run, update_orbs, obj_func, update_orbs_basis, hess_x
                test_passed &= check_released(
                    releases, "test_density_qn_factory_py", case_name
                )

        self.assertTrue(test_passed, "test_density_qn_factory_py failed")
        print(" test_density_qn_factory_py PASSED")

    def test_density_qn_settings_layout(self):
        """
        this function checks that the Python density QN settings structure has the size
        of the bind(C) settings type and that every field is read back under its own
        name from settings whose fields Fortran sets one by one to the reference values
        """
        lib.density_qn_settings_size.restype = c_size_t
        self.assertTrue(
            check_settings_layout(
                "test_density_qn_settings_layout",
                DensityQNSettings,
                lib.get_reference_density_qn_values,
                lib.reference_density_qn_field,
                lib.density_qn_settings_size(),
            ),
            "test_density_qn_settings_layout failed",
        )
        print(" test_density_qn_settings_layout PASSED")
