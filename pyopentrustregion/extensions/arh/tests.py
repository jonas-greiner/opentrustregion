# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import unittest
from ctypes import (
    c_bool,
    c_char,
    c_void_p,
    Array,
    byref,
    CFUNCTYPE,
    Structure,
    POINTER,
)
from unittest.mock import patch

from pyopentrustregion.tests import (
    lib,
    n_param,
    n_extra_trial_vectors,
    NUMPY_AVAILABLE,
    add_tests,
    print_separator,
    PyInterfaceTests,
)
from pyopentrustregion.python_interface import (
    c_real,
    c_int,
    SolverSettings,
)
from pyopentrustregion.extensions.arh import (
    ARHSettings,
    arh_factory,
    arh_factory_mo,
    arh_factory_oao,
    arh_deconstructor,
)
from pyopentrustregion.extensions.oao.tests import n_ao

if NUMPY_AVAILABLE:
    import numpy as np


# define all tests in alphabetical order
fortran_tests = {
    "arh_c_system_tests": ["arh_settings_init"],
    "arh_tests": [
        "apply_ms_sr1_skip",
        "arh_deconstructor",
        "arh_factory_mo_cs",
        "arh_factory_mo_os",
        "arh_factory_oao_cs",
        "arh_factory_oao_os",
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
        "get_extra_trial_vectors_arh_callback",
        "get_low_rank_hess_factors",
        "get_ms_a_inv",
        "get_ms_a_inv_os_linear",
        "hess_x_arh_callback",
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
        "obj_func_arh_cs_callback",
        "obj_func_arh_os_callback",
        "packed_channel_rows_arh_mo",
        "packed_channel_rows_arh_oao",
        "precond_arh_callback",
        "precond_pd_arh_callback",
        "prepend",
        "rebase_dirs",
        "rebuild_stale_hess_model",
        "resolvable_residual",
        "response_gram",
        "response_gram_os_linear",
        "rotate_trial_arh_mo",
        "rotate_trial_arh_oao",
        "spectral_to_dense",
        "to_history_basis_arh_mo",
        "to_history_basis_arh_oao",
        "truncated_eigval_inv",
        "update_orbs_arh_cs_callback",
        "update_orbs_arh_os_callback",
    ],
    "arh_c_interface_tests": [
        "arh_deconstructor_c_wrapper",
        "arh_factory_mo_c_wrapper",
        "arh_factory_oao_c_wrapper",
        "assign_arh_c_f",
        "assign_arh_f_c",
        "evaluate_dm_cs_f_wrapper",
        "evaluate_dm_os_f_wrapper",
        "hess_x_arh_c_wrapper",
        "init_arh_settings_c",
        "update_orbs_arh_c_wrapper",
    ],
}

# multiples of the density matrix the mock density matrix evaluating functions return
# for each optional output, in the order of their argument lists
evaluate_dm_cs_factors = list((c_real * 2).in_dll(lib, "test_evaluate_dm_cs_factors"))
evaluate_dm_os_factors = list((c_real * 4).in_dll(lib, "test_evaluate_dm_os_factors"))

# number of MOs and occupied orbitals
n_mo = c_int.in_dll(lib, "test_n_mo").value
n_occ = list((c_int * 2).in_dll(lib, "test_n_occ"))


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

    assign_ref_to_settings = PyInterfaceTests.assign_ref_to_settings
    equal_settings_to_ref = PyInterfaceTests.equal_settings_to_ref

    def mock_evaluate_dm_cs(self, dm_ao, fock, v_nonlinear):
        """
        this function is a mock function for the density matrix evaluating function
        with a separate non-linear potential contribution for the closed-shell case
        """
        if fock is not None:
            fock[:] = evaluate_dm_cs_factors[0] * dm_ao
        if v_nonlinear is not None:
            v_nonlinear[:] = evaluate_dm_cs_factors[1] * dm_ao

        return np.sum(dm_ao)

    def mock_evaluate_dm_os(
        self, dm_ao, fock, v_same_spin, v_opposite_spin, v_nonlinear
    ):
        """
        this function is a mock function for the density matrix evaluating function
        with separate same- and opposite spin and non-linear contributions for the
        open-shell case
        """
        if fock is not None:
            fock[:] = evaluate_dm_os_factors[0] * dm_ao
        if v_same_spin is not None:
            v_same_spin[:] = evaluate_dm_os_factors[1] * dm_ao
        if v_opposite_spin is not None:
            v_opposite_spin[:] = evaluate_dm_os_factors[2] * dm_ao
        if v_nonlinear is not None:
            v_nonlinear[:] = evaluate_dm_os_factors[3] * dm_ao

        return np.sum(dm_ao)

    mock_logger = PyInterfaceTests.mock_logger

    # replace original library with mock library
    @patch(
        "pyopentrustregion.python_interface.lib.arh_factory_mo",
        lib.mock_arh_factory_mo,
    )
    def test_arh_factory_mo_py_interface(self):
        """
        this function tests the ARH factory python interface for orbitals
        parameterized in the MO basis
        """
        ao_overlap = np.full(2 * (n_ao,), 2.0, dtype=np.float64)

        # initialize test flag
        test_passed = True

        # initialize settings object
        settings = ARHSettings()
        self.assign_ref_to_settings(settings)

        def mo_coeff_pattern(n_particle, offset):
            # MO coefficients encoding their particle channel, AO and MO index
            k, i, j = np.meshgrid(
                np.arange(1, n_particle + 1),
                np.arange(1, n_ao + 1),
                np.arange(1, n_mo + 1),
                indexing="ij",
            )
            return (offset + 100 * k + 10 * i + j).astype(np.float64)

        # MO coefficients stored row-major are copied back, those stored column-major
        # are rotated in place and those of a different precision are converted
        cases = (
            ("closed-shell row-major", 1, False, np.float64, self.mock_evaluate_dm_cs),
            (
                "closed-shell column-major",
                1,
                True,
                np.float64,
                self.mock_evaluate_dm_cs,
            ),
            (
                "closed-shell single-precision",
                1,
                True,
                np.float32,
                self.mock_evaluate_dm_cs,
            ),
            ("open-shell row-major", 2, False, np.float64, self.mock_evaluate_dm_os),
            ("open-shell column-major", 2, True, np.float64, self.mock_evaluate_dm_os),
        )
        mock_passed = c_bool.in_dll(lib, "test_arh_factory_mo_interface")
        for case, n_particle, column_major, dtype, evaluate_dm in cases:
            # initialize MO coefficients
            mo_coeff = mo_coeff_pattern(n_particle, 0.0).astype(dtype)
            if n_particle == 1:
                mo_coeff = mo_coeff[0]
            if column_major:
                mo_coeff = np.ascontiguousarray(mo_coeff.swapaxes(-1, -2)).swapaxes(
                    -1, -2
                )

            # initialize logging boolean
            self.test_logger = True

            # initialize solver settings object, without a projection for the
            # closed-shell cases and with a projection supplied by the caller, which
            # has to be kept, for the open-shell cases, already passed on as a C
            # callback by a previous solver call
            solver_settings = SolverSettings()
            caller_project = None if n_particle == 1 else lambda vector: None
            solver_settings.project = caller_project
            solver_settings.stability_settings.project = caller_project
            solver_settings.set_optional_callbacks(n_param, {})

            # call ARH MO factory python interface
            obj_func_arh, update_orbs_arh = arh_factory_mo(
                mo_coeff,
                ao_overlap,
                n_occ[0] if n_particle == 1 else n_occ[:n_particle],
                n_particle,
                n_ao,
                n_mo,
                evaluate_dm,
                solver_settings,
                settings,
            )

            # determine if the ARH routines are wired into the solver settings without
            # a projection
            if any(
                callback is None
                for callback in (
                    solver_settings.precond,
                    solver_settings.precond_pd,
                    solver_settings.get_extra_trial_vectors,
                    solver_settings.stability_settings.precond,
                    solver_settings.stability_settings.get_extra_trial_vectors,
                )
            ):
                print(
                    " test_arh_factory_mo_py_interface failed: ARH routines not wired "
                    f"into solver settings for the {case} case."
                )
                test_passed = False
            if (
                solver_settings.project is not caller_project
                or solver_settings.stability_settings.project is not caller_project
            ):
                print(
                    " test_arh_factory_mo_py_interface failed: Projection wired or "
                    f"projection supplied by the caller replaced for the {case} case."
                )
                test_passed = False
            if not solver_settings.refresh_hess or solver_settings.hess_symm:
                print(
                    " test_arh_factory_mo_py_interface failed: Hessian refresh or "
                    f"symmetry not passed on for the {case} case."
                )
                test_passed = False
            if solver_settings.n_micro != 300:
                print(
                    " test_arh_factory_mo_py_interface failed: Micro iteration limit "
                    f"not passed on for the {case} case."
                )
                test_passed = False

            # check if the mock factory received the correct input and reset its flag
            if not mock_passed.value:
                print(
                    " test_arh_factory_mo_py_interface failed: Mock factory received "
                    f"wrong input for the {case} case."
                )
                test_passed = False
            mock_passed.value = True

            # check if logger was called correctly
            if not self.test_logger:
                print(
                    " test_arh_factory_mo_py_interface failed: Called logging "
                    f"function wrong for the {case} case."
                )
                test_passed = False

            # call returned ARH objective function
            kappa = np.ones(n_param, dtype=np.float64)
            func = obj_func_arh(kappa)
            if func != 3.0:
                print(
                    " test_arh_factory_mo_py_interface failed: Returned function "
                    f"value of returned ARH objective function wrong for the {case} "
                    "case."
                )
                test_passed = False

            # call returned ARH orbital updating function
            grad = np.empty(n_param, dtype=np.float64)
            h_diag = np.empty(n_param, dtype=np.float64)
            func, hess_x_funptr = update_orbs_arh(kappa, grad, h_diag)
            if (
                func != 3.0
                or not np.allclose(grad, 2.0)
                or not np.allclose(h_diag, 3.0)
            ):
                print(
                    " test_arh_factory_mo_py_interface failed: Returned values of "
                    f"returned ARH orbital updating function wrong for the {case} "
                    "case."
                )
                test_passed = False

            # call returned hess_x function
            hess_x = np.empty(n_param, dtype=np.float64)
            hess_x_funptr(np.ones(n_param, dtype=np.float64), hess_x)
            if not np.allclose(hess_x, 4.0):
                print(
                    " test_arh_factory_mo_py_interface failed: Returned Hessian "
                    "linear transformation of returned ARH orbital updating function "
                    f"wrong for the {case} case."
                )
                test_passed = False

            # check if the rotated MO coefficients reached the caller
            expected = mo_coeff_pattern(n_particle, 1000.0)
            if n_particle == 1:
                expected = expected[0]
            if not np.allclose(mo_coeff, expected):
                print(
                    " test_arh_factory_mo_py_interface failed: MO coefficients not "
                    f"updated correctly for the {case} case."
                )
                test_passed = False

            # call wired ARH level-shifted preconditioner function
            residual = np.full(n_param, 1.0, dtype=np.float64)
            precond_residual = np.empty(n_param, dtype=np.float64)
            mu = 5.0
            solver_settings.precond(residual, mu, precond_residual)
            if not np.allclose(precond_residual, mu):
                print(
                    " test_arh_factory_mo_py_interface failed: Returned "
                    "preconditioned residual of wired ARH level-shifted "
                    f"preconditioner function wrong for the {case} case."
                )
                test_passed = False

            # call wired ARH positive-definite preconditioner function
            precond_pd_residual = np.empty(n_param, dtype=np.float64)
            solver_settings.precond_pd(residual, precond_pd_residual)
            if not np.allclose(precond_pd_residual, 3.0):
                print(
                    " test_arh_factory_mo_py_interface failed: Returned "
                    "preconditioned residual of wired ARH positive-definite "
                    f"preconditioner function wrong for the {case} case."
                )
                test_passed = False

            # call wired ARH extra trial vector function
            trial_vectors = np.empty((n_extra_trial_vectors, n_param), dtype=np.float64)
            solver_settings.get_extra_trial_vectors(trial_vectors)
            if not np.allclose(
                trial_vectors,
                np.arange(1, n_extra_trial_vectors + 1, dtype=np.float64)[:, None],
            ):
                print(
                    " test_arh_factory_mo_py_interface failed: Returned trial "
                    "vectors of wired ARH extra trial vector function wrong for "
                    f"the {case} case."
                )
                test_passed = False

        # invalid inputs have to raise an error before the library is called
        read_only_mo_coeff = mo_coeff_pattern(1, 0.0)[0]
        read_only_mo_coeff.flags.writeable = False
        cs_mo_coeff = mo_coeff_pattern(1, 0.0)[0]
        cs_evaluate_dm = self.mock_evaluate_dm_cs
        invalid_cases = (
            (
                "read-only MO coefficients",
                (
                    read_only_mo_coeff,
                    ao_overlap,
                    n_occ[0],
                    1,
                    n_ao,
                    n_mo,
                    cs_evaluate_dm,
                ),
                ValueError,
            ),
            (
                "complex MO coefficients",
                (
                    cs_mo_coeff.astype(np.complex128),
                    ao_overlap,
                    n_occ[0],
                    1,
                    n_ao,
                    n_mo,
                    cs_evaluate_dm,
                ),
                ValueError,
            ),
            (
                "integer MO coefficients",
                (
                    cs_mo_coeff.astype(np.int64),
                    ao_overlap,
                    n_occ[0],
                    1,
                    n_ao,
                    n_mo,
                    cs_evaluate_dm,
                ),
                ValueError,
            ),
            (
                "missing occupations",
                (
                    mo_coeff_pattern(2, 0.0),
                    ao_overlap,
                    n_occ[0],
                    2,
                    n_ao,
                    n_mo,
                    self.mock_evaluate_dm_os,
                ),
                ValueError,
            ),
            (
                "MO coefficients not matching the number of particles",
                (
                    cs_mo_coeff,
                    ao_overlap,
                    n_occ,
                    2,
                    n_ao,
                    n_mo,
                    self.mock_evaluate_dm_os,
                ),
                ValueError,
            ),
            (
                "MO coefficients not matching the number of AOs",
                (
                    cs_mo_coeff,
                    np.full(2 * (n_ao + 1,), 2.0),
                    n_occ[0],
                    1,
                    n_ao + 1,
                    n_mo,
                    cs_evaluate_dm,
                ),
                ValueError,
            ),
            (
                "MO coefficients not matching the number of MOs",
                (cs_mo_coeff, ao_overlap, n_occ[0], 1, n_ao, n_mo + 1, cs_evaluate_dm),
                ValueError,
            ),
            (
                "wrongly shaped AO overlap matrix",
                (
                    cs_mo_coeff,
                    ao_overlap[:-1, :-1],
                    n_occ[0],
                    1,
                    n_ao,
                    n_mo,
                    cs_evaluate_dm,
                ),
                ValueError,
            ),
            (
                "non-integer occupations",
                (
                    cs_mo_coeff,
                    ao_overlap,
                    n_occ[0] + 0.5,
                    1,
                    n_ao,
                    n_mo,
                    cs_evaluate_dm,
                ),
                TypeError,
            ),
            (
                "open-shell evaluate_dm for closed-shell MO coefficients",
                (
                    cs_mo_coeff,
                    ao_overlap,
                    n_occ[0],
                    1,
                    n_ao,
                    n_mo,
                    self.mock_evaluate_dm_os,
                ),
                TypeError,
            ),
        )
        for case, args, error_type in invalid_cases:
            try:
                arh_factory_mo(*args, SolverSettings(), settings)
                print(
                    f" test_arh_factory_mo_py_interface failed: {case} did not raise "
                    "an error."
                )
                test_passed = False
            except error_type:
                pass

        self.assertTrue(test_passed, "test_arh_factory_mo_py_interface failed")
        print(" test_arh_factory_mo_py_interface PASSED")

    # replace original library with mock library
    @patch(
        "pyopentrustregion.python_interface.lib.arh_factory_oao",
        lib.mock_arh_factory_oao,
    )
    def test_arh_factory_oao_py_interface(self):
        """
        this function tests the ARH factory python interface (only tests whether dm_ao
        and mock_evaluate_dm_os are passed correctly for the open-shell case since
        everything else is the same in the closed-shell case)
        """
        ao_overlap = np.full(2 * (n_ao,), 2.0, dtype=np.float64)

        # initialize test flag
        test_passed = True

        # initialize settings object
        settings = ARHSettings()
        self.assign_ref_to_settings(settings)

        # initialize logging boolean
        self.test_logger = True

        # initialize density matrix
        dm_ao = np.full(2 * (n_ao,), 1.0, dtype=np.float64)

        # initialize solver settings object
        solver_settings = SolverSettings()

        # call ARH factory python interface
        obj_func_arh, update_orbs_arh = arh_factory_oao(
            dm_ao,
            ao_overlap,
            1,
            n_ao,
            self.mock_evaluate_dm_cs,
            solver_settings,
            settings,
        )

        # determine if the ARH routines are wired into the solver settings
        if any(
            callback is None
            for callback in (
                solver_settings.precond,
                solver_settings.precond_pd,
                solver_settings.project,
                solver_settings.get_extra_trial_vectors,
                solver_settings.stability_settings.precond,
                solver_settings.stability_settings.project,
                solver_settings.stability_settings.get_extra_trial_vectors,
            )
        ):
            print(
                " test_arh_factory_oao_py_interface failed: ARH routines not wired "
                "into solver settings."
            )
            test_passed = False
        if not solver_settings.refresh_hess:
            print(
                " test_arh_factory_oao_py_interface failed: Hessian refresh not "
                "requested."
            )
            test_passed = False
        if solver_settings.hess_symm:
            print(
                " test_arh_factory_oao_py_interface failed: Symmetry of the "
                "approximate Hessian not passed on."
            )
            test_passed = False
        if solver_settings.n_micro != 300:
            print(
                " test_arh_factory_oao_py_interface failed: Micro iteration limit not "
                "passed on."
            )
            test_passed = False

        # check if logger was called correctly
        if not self.test_logger:
            print(
                " test_arh_factory_oao_py_interface failed: Called logging function "
                "wrong."
            )
            test_passed = False

        # call returned ARH objective function
        kappa = np.ones(n_param, dtype=np.float64)
        try:
            func = obj_func_arh(kappa)
        except RuntimeError:
            print(
                " test_arh_factory_oao_py_interface failed: Returned ARH objective "
                "function raises error."
            )
            test_passed = False

        # check results
        if func != 3.0:
            print(
                " test_arh_factory_oao_py_interface failed: Returned function value of "
                "returned ARH objective function wrong."
            )
            test_passed = False

        # call returned ARH orbital updating function
        grad = np.empty(n_param, dtype=np.float64)
        h_diag = np.empty(n_param, dtype=np.float64)
        try:
            func, hess_x_funptr = update_orbs_arh(kappa, grad, h_diag)
        except RuntimeError:
            print(
                " test_arh_factory_oao_py_interface failed: Returned ARH orbital "
                "updating function raises error."
            )
            test_passed = False

        # check if density matrix was updated
        if not np.allclose(dm_ao, np.full(2 * (n_ao,), 2.0, dtype=np.float64)):
            print(
                " test_arh_factory_oao_py_interface failed: Density matrix not updated "
                "correctly."
            )
            test_passed = False

        # check results
        if func != 3.0:
            print(
                " test_arh_factory_oao_py_interface failed: Returned function value of "
                "returned ARH orbital updating function wrong."
            )
            test_passed = False
        if not np.allclose(grad, np.full(n_param, 2.0, dtype=np.float64)):
            print(
                " test_arh_factory_oao_py_interface failed: Returned gradient of "
                "returned ARH orbital updating function wrong."
            )
            test_passed = False
        if not np.allclose(h_diag, np.full(n_param, 3.0, dtype=np.float64)):
            print(
                " test_arh_factory_oao_py_interface failed: Returned ARH Hessian "
                "diagonal of returned ARH orbital updating function wrong."
            )
            test_passed = False

        # call returned hess_x function
        x = np.ones(n_param, dtype=np.float64)
        hess_x = np.empty(n_param, dtype=np.float64)

        try:
            hess_x_funptr(x, hess_x)
        except RuntimeError:
            print(
                " test_arh_factory_oao_py_interface failed: Returned ARH Hessian "
                "linear transformation function of returned ARH orbital updating "
                "function raises error."
            )
            test_passed = False

        # check results
        if not np.allclose(hess_x, np.full(n_param, 4.0, dtype=np.float64)):
            print(
                " test_arh_factory_oao_py_interface failed: Returned Hessian linear "
                "transformation of returned ARH orbital updating function wrong."
            )
            test_passed = False

        # call wired ARH projection function
        vector = np.full(n_param, 1.0, dtype=np.float64)
        try:
            solver_settings.project(vector)
        except RuntimeError:
            print(
                " test_arh_factory_oao_py_interface failed: Wired ARH projection "
                "function raises error."
            )
            test_passed = False

        # check results
        if not np.allclose(vector, np.full(n_param, 2.0, dtype=np.float64)):
            print(
                " test_arh_factory_oao_py_interface failed: Returned projected vector "
                "of wired ARH projection function wrong."
            )
            test_passed = False

        # call wired ARH level-shifted preconditioner function
        residual = np.full(n_param, 1.0, dtype=np.float64)
        precond_residual = np.empty(n_param, dtype=np.float64)
        mu = 5.0
        try:
            solver_settings.precond(residual, mu, precond_residual)
        except RuntimeError:
            print(
                " test_arh_factory_oao_py_interface failed: Wired ARH level-shifted "
                "preconditioner function raises error."
            )
            test_passed = False

        # check results
        if not np.allclose(precond_residual, np.full(n_param, mu, dtype=np.float64)):
            print(
                " test_arh_factory_oao_py_interface failed: Returned preconditioned "
                "residual of wired ARH level-shifted preconditioner function "
                "wrong."
            )
            test_passed = False

        # call wired ARH positive-definite preconditioner function
        precond_pd_residual = np.empty(n_param, dtype=np.float64)
        try:
            solver_settings.precond_pd(residual, precond_pd_residual)
        except RuntimeError:
            print(
                " test_arh_factory_oao_py_interface failed: Wired ARH "
                "positive-definite preconditioner function raises error."
            )
            test_passed = False

        # check results
        if not np.allclose(
            precond_pd_residual, np.full(n_param, 3.0, dtype=np.float64)
        ):
            print(
                " test_arh_factory_oao_py_interface failed: Returned preconditioned "
                "residual of wired ARH positive-definite preconditioner function "
                "wrong."
            )
            test_passed = False

        # number of particles
        n_particle = 2

        # initialize density matrix
        dm_ao = np.full((n_particle, n_ao, n_ao), 1.0, dtype=np.float64)

        # call ARH factory python interface
        arh_factory_oao(
            dm_ao,
            ao_overlap,
            n_particle,
            n_ao,
            self.mock_evaluate_dm_os,
            SolverSettings(),
            settings,
        )

        # call arh_factory_oao python interface for invalid density and overlap
        # matrices, which have to raise an error before the Fortran factory is called
        read_only_dm_ao = np.full(2 * (n_ao,), 1.0, dtype=np.float64)
        read_only_dm_ao.flags.writeable = False
        invalid_cases = (
            (
                "density matrix not matching the number of particles",
                (
                    np.full(2 * (n_ao,), 1.0),
                    ao_overlap,
                    2,
                    n_ao,
                    self.mock_evaluate_dm_os,
                ),
                ValueError,
            ),
            (
                "density matrix not matching the number of AOs",
                (
                    np.full(2 * (n_ao,), 1.0),
                    np.full(2 * (n_ao + 1,), 2.0),
                    1,
                    n_ao + 1,
                    self.mock_evaluate_dm_cs,
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
                    self.mock_evaluate_dm_cs,
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
                    self.mock_evaluate_dm_cs,
                ),
                ValueError,
            ),
            (
                "read-only density matrix",
                (read_only_dm_ao, ao_overlap, 1, n_ao, self.mock_evaluate_dm_cs),
                ValueError,
            ),
            (
                "wrongly shaped AO overlap matrix",
                (
                    np.full(2 * (n_ao,), 1.0),
                    ao_overlap[:-1, :-1],
                    1,
                    n_ao,
                    self.mock_evaluate_dm_cs,
                ),
                ValueError,
            ),
            (
                "open-shell evaluate_dm for a closed-shell density matrix",
                (
                    np.full(2 * (n_ao,), 1.0),
                    ao_overlap,
                    1,
                    n_ao,
                    self.mock_evaluate_dm_os,
                ),
                TypeError,
            ),
        )
        for case, args, error_type in invalid_cases:
            try:
                arh_factory_oao(*args, SolverSettings(), settings)
                print(
                    f" test_arh_factory_oao_py_interface failed: {case} did not raise "
                    "an error."
                )
                test_passed = False
            except error_type:
                pass

        self.assertTrue(
            c_bool.in_dll(lib, "test_arh_factory_oao_interface").value and test_passed,
            "test_arh_factory_oao_py_interface failed",
        )
        print(" test_arh_factory_oao_py_interface PASSED")

    def test_arh_factory_py_interface(self):
        """
        this function tests the ARH factory python interface which, like the generic
        Fortran factory, selects the orbital basis from the arguments, by determining
        if it forwards arguments matching those of one basis factory to that factory,
        positionally or by keyword, and rejects every other argument list
        """
        module = "pyopentrustregion.extensions.arh.python_interface"
        ao_overlap = np.full(2 * (n_ao,), 2.0, dtype=np.float64)
        dm_ao = np.full(2 * (n_ao,), 1.0, dtype=np.float64)
        mo_coeff = np.full((n_ao, n_mo), 1.0, dtype=np.float64)
        evaluate_dm = self.mock_evaluate_dm_cs
        solver_settings = SolverSettings()
        settings = ARHSettings()
        oao_args = (dm_ao, ao_overlap, 1, n_ao, evaluate_dm, solver_settings, settings)
        mo_args = (
            mo_coeff,
            ao_overlap,
            n_occ[0],
            1,
            n_ao,
            n_mo,
            evaluate_dm,
            solver_settings,
            settings,
        )
        oao_names = ("dm_ao", "ao_overlap", "n_particle", "n_ao", "evaluate_dm")
        oao_names += ("solver_settings", "settings")

        # initialize test flag
        test_passed = True

        # every valid argument list, the factory it selects and the positional and
        # keyword arguments this factory has to receive
        valid_cases = (
            ("OAO arguments", oao_args, {}, "arh_factory_oao"),
            ("MO arguments", mo_args, {}, "arh_factory_mo"),
            (
                "OAO keyword arguments",
                (),
                dict(zip(oao_names, oao_args)),
                "arh_factory_oao",
            ),
            (
                "MO arguments with keyword settings",
                mo_args[:-2],
                {"solver_settings": solver_settings, "settings": settings},
                "arh_factory_mo",
            ),
        )
        for case, args, kwargs, selected in valid_cases:
            with patch(f"{module}.arh_factory_mo", return_value="mo") as mock_mo, patch(
                f"{module}.arh_factory_oao", return_value="oao"
            ) as mock_oao:
                result = arh_factory(*args, **kwargs)
                mocks = {"arh_factory_mo": mock_mo, "arh_factory_oao": mock_oao}
            if result != mocks[selected].return_value:
                print(
                    f" test_arh_factory_py_interface failed: Result of {selected} not "
                    f"returned for {case}."
                )
                test_passed = False
            if mocks[selected].call_count != 1 or any(
                mock.called for name, mock in mocks.items() if name != selected
            ):
                print(
                    f" test_arh_factory_py_interface failed: {selected} not selected "
                    f"alone for {case}."
                )
                test_passed = False
            else:
                call = mocks[selected].call_args
                if (
                    len(call.args) != len(args)
                    or any(a is not b for a, b in zip(call.args, args))
                    or call.kwargs.keys() != kwargs.keys()
                    or any(call.kwargs[key] is not kwargs[key] for key in kwargs)
                ):
                    print(
                        f" test_arh_factory_py_interface failed: Arguments not "
                        f"forwarded to {selected} for {case}."
                    )
                    test_passed = False

        # every argument list matching neither factory, which has to raise an error
        # without calling a factory
        invalid_cases = (
            ("too few arguments", oao_args[:-1], {}),
            ("an argument count between the two factories", mo_args[:-1], {}),
            ("too many arguments", mo_args + (None,), {}),
            (
                "keyword arguments of both factories",
                (),
                {**dict(zip(oao_names, oao_args)), "n_mo": n_mo},
            ),
        )
        for case, args, kwargs in invalid_cases:
            with patch(f"{module}.arh_factory_mo") as mock_mo, patch(
                f"{module}.arh_factory_oao"
            ) as mock_oao:
                try:
                    arh_factory(*args, **kwargs)
                    print(
                        f" test_arh_factory_py_interface failed: {case} did not raise "
                        "an error."
                    )
                    test_passed = False
                except TypeError:
                    pass
            if mock_mo.called or mock_oao.called:
                print(
                    f" test_arh_factory_py_interface failed: Factory called for {case}."
                )
                test_passed = False

        self.assertTrue(test_passed, "test_arh_factory_py_interface failed")
        print(" test_arh_factory_py_interface PASSED")

    # replace original library with mock library
    @patch(
        "pyopentrustregion.python_interface.lib.arh_deconstructor",
        lib.mock_arh_deconstructor,
    )
    def test_arh_deconstructor_py_interface(self):
        """
        this function tests the ARH deconstructor python interface
        """
        # initialize test flag
        test_passed = True

        # call ARH deconstructor python interface
        arh_deconstructor()

        # check if deconstructor was called correctly
        if not c_bool.in_dll(lib, "test_arh_deconstructor_interface").value:
            print(
                " test_arh_deconstructor_py_interface failed: Deconstructor called "
                "wrong."
            )
            test_passed = False

        self.assertTrue(test_passed, "test_arh_deconstructor_py_interface failed")
        print(" test_arh_deconstructor_py_interface PASSED")

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

        dummy_error_code = 42

        def dummy_logger():
            return dummy_error_code

        settings.set_optional_callback(
            "logger", dummy_logger, lambda x: x, CFUNCTYPE(c_int)
        )

        c_ptr = getattr(settings.settings_c, "logger")
        c_interface = getattr(settings.settings_c, "logger_interface", None)

        if (
            c_ptr is None
            or (isinstance(c_ptr, c_void_p) and c_ptr.value is None)
            or not callable(c_interface)
            or c_interface() != dummy_error_code
        ):
            print(
                " test_arh_settings failed: Optional callbacks are not set correctly."
            )
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
