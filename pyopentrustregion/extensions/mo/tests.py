# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import gc
import unittest
import weakref
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
    c_int,
    update_orbs_interface_type,
    SolverSettings,
    SolverState,
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
from pyopentrustregion.extensions.mo import MOSettings, mo_factory

# define all tests in alphabetical order
fortran_tests = {
    "mo_c_integration_tests": [
        "mo_factory_c",
        "mo_settings_init",
        "mo_settings_layout",
    ],
    "mo_tests": [
        "calculate_grad_h_diag_mo",
        "count_mo_params",
        "diagonalize_per_irrep",
        "get_extra_trial_vectors_mo",
        "get_hess_eigval_pairs_mo",
        "hess_x_mo",
        "hess_x_static_mo",
        "mo_factory_common",
        "mo_factory_cs",
        "mo_factory_os",
        "mo_param_rows",
        "mo_sanity_check",
        "mo_transform",
        "obj_func_mo",
        "pack_ov",
        "refresh_hess_eigen_mo",
        "rotate_from_hess_eigenbasis_mo",
        "rotate_mo_coeff",
        "rotate_orbitals_mo",
        "rotate_to_hess_eigenbasis_mo",
        "unpack_ov",
        "update_orbs_mo",
    ],
    "mo_c_interface_tests": [
        "assign_mo_c_f",
        "assign_mo_f_c",
        "init_mo_settings_c",
        "mo_factory_c_wrapper",
    ],
}

# number of MOs and of occupied orbitals of every particle channel
n_mo = c_int.in_dll(lib, "test_n_mo").value
n_occ = list((c_int * 2).in_dll(lib, "test_n_occ"))

# irreps of the MOs of every particle channel, those of the first channel serving the
# closed-shell case
orbsym = [list(irreps) for irreps in ((c_int * n_mo) * 2).in_dll(lib, "test_orbsym")]


def setup_mo_toy_run(n_particle):
    """
    this function prepares a solve of the toy model in the MO basis and the given shell
    through the MO factory and returns the run and the functions the factory returned,
    also for the integration tests of density QN
    """
    run = ToyRun(True, n_particle)
    update_orbs, obj_func = mo_factory(
        run.orbitals,
        toy_model().s,
        run.n_occ,
        n_particle,
        n_ao,
        run.n_mo,
        run.evaluate_dm,
        run.settings,
        run.state,
    )
    return run, update_orbs, obj_func


@add_tests
class MOTests(unittest.TestCase):
    """
    this class contains unit tests for MO
    """

    tests = fortran_tests["mo_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for MO...")
        return super().setUpClass()


@add_tests
class MOCInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the MO C interface
    """

    tests = fortran_tests["mo_c_interface_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for MO C interface...")
        return super().setUpClass()


class MOPyInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the Python interface
    """

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for MO Python interface...")

        return super().setUpClass()

    mock_evaluate_dm = OrbitalCallbackMocks.mock_evaluate_dm
    mock_get_response = OrbitalCallbackMocks.mock_get_response

    # replace original library with mock library, also for the release of the
    # extensions the factory returns
    @patch("pyopentrustregion.python_interface.lib.mo_factory", lib.mock_mo_factory)
    @patch(
        "pyopentrustregion.extensions.common.python_interface.release_extension",
        release_mock_extension,
    )
    def test_mo_factory_py_interface(self):
        """
        this function tests the MO factory python interface
        """
        ao_overlap = ao_overlap_pattern()

        # initialize test flag
        test_passed = True

        # the logging function records the messages, and the mock factory records in a
        # flag whether it received what it expects, which is cleared before and read
        # after every call of it
        messages = []
        mock_passed = c_bool.in_dll(lib, "test_mo_factory_interface")

        # initialize settings object
        settings = MOSettings()
        assign_reference_settings(settings, lib.reference_mo_field)
        settings.logger = messages.append

        def mo_coeff_pattern(n_particle, offset):
            # MO coefficients encoding their particle channel, AO and MO index
            k, i, j = np.meshgrid(
                np.arange(1, n_particle + 1),
                np.arange(1, n_ao + 1),
                np.arange(1, n_mo + 1),
                indexing="ij",
            )
            return (offset + 100 * k + 10 * i + j).astype(np.float64)

        def initial_mo_coeff(n_particle, column_major, dtype):
            # MO coefficients of the given shell, storage order and precision
            mo_coeff = mo_coeff_pattern(n_particle, 0.0).astype(dtype)
            if n_particle == 1:
                mo_coeff = mo_coeff[0]
            if column_major:
                mo_coeff = np.ascontiguousarray(mo_coeff.swapaxes(-1, -2)).swapaxes(
                    -1, -2
                )
            return mo_coeff

        # MO coefficients stored row-major are copied back, those stored column-major
        # are rotated in place and those of a different precision are converted, like
        # the AO overlap matrix of the same precision, with the irreps of their MOs
        # passed on only if they are given
        cases = (
            ("closed-shell row-major", 1, False, np.float64, False),
            ("closed-shell column-major", 1, True, np.float64, False),
            ("closed-shell single-precision", 1, True, np.float32, False),
            ("open-shell row-major", 2, False, np.float64, False),
            ("open-shell column-major", 2, True, np.float64, False),
            ("closed-shell with symmetry", 1, False, np.float64, True),
            ("open-shell with symmetry", 2, True, np.float64, True),
        )
        orbsym_passed = c_bool.in_dll(lib, "test_orbsym_passed")
        for case, n_particle, column_major, dtype, with_irreps in cases:
            # initialize MO coefficients
            mo_coeff = initial_mo_coeff(n_particle, column_major, dtype)

            # clear the recorded messages
            messages.clear()

            # initialize solver settings object
            solver_settings, python_functions = setup_wired_solver_settings()

            # call MO factory python interface
            mock_passed.value = False
            update_orbs_mo, obj_func_mo = mo_factory(
                mo_coeff,
                ao_overlap.astype(dtype),
                n_occ[0] if n_particle == 1 else n_occ[:n_particle],
                n_particle,
                n_ao,
                n_mo,
                self.mock_evaluate_dm,
                solver_settings,
                SolverState(),
                settings,
                orbsym=(
                    (orbsym[0] if n_particle == 1 else orbsym) if with_irreps else None
                ),
            )

            # determine if the irreps of the MOs were passed on only if they are given
            if orbsym_passed.value != with_irreps:
                print(
                    " test_mo_factory_py_interface failed: Irreps of the MOs passed "
                    f"on wrongly for the {case} case."
                )
                test_passed = False

            # check the solver settings the factory wired
            if not check_wired_solver_settings(
                solver_settings,
                python_functions,
                "test_mo_factory_py_interface",
                (update_orbs_mo, obj_func_mo),
                settings,
            ):
                test_passed = False

            # check if the mock factory received the correct input
            if not mock_passed.value:
                print(
                    " test_mo_factory_py_interface failed: Mock factory received wrong "
                    f"input for the {case} case."
                )
                test_passed = False

            # check that the settings and the density matrix evaluating function were
            # passed
            if not (
                c_bool.in_dll(lib, "mo_settings_passed").value
                and c_bool.in_dll(lib, "mo_evaluate_dm_passed").value
            ):
                print(
                    " test_mo_factory_py_interface failed: Settings or density matrix "
                    f"evaluating function not passed for the {case} case."
                )
                test_passed = False

            # check if logger was called correctly
            if "test" not in messages:
                print(
                    " test_mo_factory_py_interface failed: Called logging function "
                    f"wrong for the {case} case."
                )
                test_passed = False

            # call returned MO objective function
            kappa = np.ones(n_param, dtype=np.float64)
            func = obj_func_mo(kappa)
            if func != 3.0:
                print(
                    " test_mo_factory_py_interface failed: Returned function value "
                    f"of returned MO objective function wrong for the {case} case."
                )
                test_passed = False

            # call returned MO orbital updating function
            grad = np.empty(n_param, dtype=np.float64)
            h_diag = np.empty(n_param, dtype=np.float64)
            func, hess_x_funptr = update_orbs_mo(kappa, grad, h_diag)
            if (
                func != 3.0
                or not np.allclose(grad, 2.0)
                or not np.allclose(h_diag, 3.0)
            ):
                print(
                    " test_mo_factory_py_interface failed: Returned values of "
                    f"returned MO orbital updating function wrong for the {case} case."
                )
                test_passed = False

            # call returned hess_x function
            hess_x = np.empty(n_param, dtype=np.float64)
            hess_x_funptr(np.ones(n_param, dtype=np.float64), hess_x)
            if not np.allclose(hess_x, 4.0):
                print(
                    " test_mo_factory_py_interface failed: Returned Hessian linear "
                    "transformation of returned MO orbital updating function wrong "
                    f"for the {case} case."
                )
                test_passed = False

            # check if MO coefficients stored column-major in double precision are
            # rotated in place and all others copied back, by the returned orbital
            # updating function and for the density QN factory by the owner of the
            # object
            in_place = column_major and dtype == np.float64
            owner = solver_settings._extension_owner
            if (
                update_orbs_mo.copy_back is None
            ) != in_place or owner.copy_back is not update_orbs_mo.copy_back:
                print(
                    " test_mo_factory_py_interface failed: MO coefficients "
                    f"{'not ' if in_place else ''}rotated in place for the {case} case."
                )
                test_passed = False

            # check if the rotated MO coefficients reached the caller
            expected = mo_coeff_pattern(n_particle, 1000.0)
            if n_particle == 1:
                expected = expected[0]
            if not np.allclose(mo_coeff, expected):
                print(
                    " test_mo_factory_py_interface failed: MO coefficients not "
                    f"updated correctly for the {case} case."
                )
                test_passed = False

        # call MO factory python interface without settings and without a density
        # matrix evaluating function, which have to be passed as none, and determine if
        # the owner of the object holds the dimensions of the orbitals for the density
        # QN factory
        settings_flag = c_bool.in_dll(lib, "mo_settings_passed")
        evaluate_dm_flag = c_bool.in_dll(lib, "mo_evaluate_dm_passed")
        settings_flag.value = True
        evaluate_dm_flag.value = True
        mock_passed.value = False
        solver_settings = setup_wired_solver_settings()[0]
        try:
            mo_factory(
                initial_mo_coeff(2, True, np.float64),
                ao_overlap,
                n_occ,
                2,
                n_ao,
                n_mo,
                None,
                solver_settings,
                SolverState(),
            )
        except RuntimeError:
            print(
                " test_mo_factory_py_interface failed: Call without settings and "
                "density matrix evaluating function raises error."
            )
            test_passed = False
        if settings_flag.value or evaluate_dm_flag.value:
            print(
                " test_mo_factory_py_interface failed: Settings or density matrix "
                "evaluating function passed although none were given."
            )
            test_passed = False

        # check if the mock factory received the correct input
        if not mock_passed.value:
            print(
                " test_mo_factory_py_interface failed: Mock factory received wrong "
                "input without settings."
            )
            test_passed = False
        owner = getattr(solver_settings, "_extension_owner", None)
        if not isinstance(owner, OrbitalOwner) or (owner.n_ao, owner.n_particle) != (
            n_ao,
            2,
        ):
            print(
                " test_mo_factory_py_interface failed: Owner of the object does not "
                "hold the dimensions of the orbitals."
            )
            test_passed = False

        # the factory hands the extension it receives to an owner, which the solver
        # settings keep alive and which keeps the state alive
        mock_passed.value = False
        if not check_extension_owner(
            lambda solver_settings, solver_state: mo_factory(
                initial_mo_coeff(1, True, np.float64),
                ao_overlap,
                n_occ[0],
                1,
                n_ao,
                n_mo,
                self.mock_evaluate_dm,
                solver_settings,
                solver_state,
                settings,
            ),
            "test_mo_factory_py_interface",
        ):
            test_passed = False

        # check if the mock factory received the correct input
        if not mock_passed.value:
            print(
                " test_mo_factory_py_interface failed: Mock factory received wrong "
                "input for the check of the owner."
            )
            test_passed = False

        # call MO factory python interface for MO coefficients stored row-major and
        # column-major, release them together with the returned functions and the
        # solver settings, and determine if neither the array owning their memory nor
        # the buffer the library rotates is kept alive
        for column_major in (False, True):
            mo_coeff = initial_mo_coeff(1, column_major, np.float64)
            mo_coeff_ref = weakref.ref(
                mo_coeff if mo_coeff.base is None else mo_coeff.base
            )
            args = (
                mo_coeff,
                ao_overlap,
                n_occ[0],
                1,
                n_ao,
                n_mo,
                self.mock_evaluate_dm,
            )
            solver_settings = setup_wired_solver_settings()[0]
            mock_passed.value = False
            returned = mo_factory(*args, solver_settings, SolverState(), settings)

            # check if the mock factory received the correct input
            if not mock_passed.value:
                print(
                    " test_mo_factory_py_interface failed: Mock factory received wrong "
                    "input for the release of the MO coefficients."
                )
                test_passed = False
            mo_coeff_buffer_ref = weakref.ref(
                solver_settings._extension_owner.saved_objects["mo_coeff_buffer"]
            )
            del mo_coeff, args, returned, solver_settings
            gc.collect()
            if mo_coeff_ref() is not None or mo_coeff_buffer_ref() is not None:
                print(
                    " test_mo_factory_py_interface failed: MO coefficients stored "
                    f"{'column' if column_major else 'row'}-major or their buffer "
                    "kept alive after being released together with the returned "
                    "functions."
                )
                test_passed = False

        # call the MO orbital updating function the factory returned for row-major MO
        # coefficients with a library routine which rotates their buffer and fails,
        # which has to raise an error and still copy the rotated MO coefficients back
        # to the caller
        mo_coeff = initial_mo_coeff(1, False, np.float64)
        solver_settings = setup_wired_solver_settings()[0]
        mock_passed.value = False
        update_orbs_mo = mo_factory(
            mo_coeff,
            ao_overlap,
            n_occ[0],
            1,
            n_ao,
            n_mo,
            self.mock_evaluate_dm,
            solver_settings,
            SolverState(),
            settings,
        )[0]
        mo_coeff_buffer = solver_settings._extension_owner.saved_objects[
            "mo_coeff_buffer"
        ]

        # check if the mock factory received the correct input
        if not mock_passed.value:
            print(
                " test_mo_factory_py_interface failed: Mock factory received wrong "
                "input for a failing orbital update."
            )
            test_passed = False

        def failing_update_orbs(
            kappa_ptr, func_ptr, grad_ptr, h_diag_ptr, hess_x_ptr, context
        ):
            mo_coeff_buffer[:] = mo_coeff_pattern(1, 1000.0)[0].T
            return 1

        update_orbs_mo.update_orbs_funptr = update_orbs_interface_type(
            failing_update_orbs
        )
        try:
            update_orbs_mo(kappa, grad, h_diag)
            print(
                " test_mo_factory_py_interface failed: Failing MO orbital updating "
                "function did not raise an error."
            )
            test_passed = False
        except RuntimeError:
            pass
        if not np.allclose(mo_coeff, mo_coeff_pattern(1, 1000.0)[0]):
            print(
                " test_mo_factory_py_interface failed: MO coefficients rotated by a "
                "failing MO orbital updating function not copied back."
            )
            test_passed = False

        # invalid inputs have to raise an error before the library is called
        read_only_mo_coeff = mo_coeff_pattern(1, 0.0)[0]
        read_only_mo_coeff.flags.writeable = False
        cs_mo_coeff = mo_coeff_pattern(1, 0.0)[0]
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
                    self.mock_evaluate_dm,
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
                    self.mock_evaluate_dm,
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
                    self.mock_evaluate_dm,
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
                    self.mock_evaluate_dm,
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
                    self.mock_evaluate_dm,
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
                    self.mock_evaluate_dm,
                ),
                ValueError,
            ),
            (
                "MO coefficients not matching the number of MOs",
                (
                    cs_mo_coeff,
                    ao_overlap,
                    n_occ[0],
                    1,
                    n_ao,
                    n_mo + 1,
                    self.mock_evaluate_dm,
                ),
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
                    self.mock_evaluate_dm,
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
                    self.mock_evaluate_dm,
                ),
                TypeError,
            ),
            (
                "evaluate_dm with a wrong number of arguments",
                (
                    cs_mo_coeff,
                    ao_overlap,
                    n_occ[0],
                    1,
                    n_ao,
                    n_mo,
                    lambda dm, fock: (0.0, None),
                ),
                TypeError,
            ),
        )
        for case, args, error_type in invalid_cases:
            try:
                mo_factory(*args, SolverSettings(), SolverState(), settings)
                print(
                    f" test_mo_factory_py_interface failed: {case} did not raise an "
                    "error."
                )
                test_passed = False
            except error_type:
                pass

        # invalid irreps of the MOs have to raise an error before the library is called
        invalid_orbsym_cases = (
            ("irreps of two particle channels for one", 1, orbsym),
            ("non-integer irreps", 1, np.array(orbsym[0], dtype=np.float64)),
            ("transposed open-shell irreps", 2, np.transpose(orbsym)),
        )
        for case, n_particle, invalid_orbsym in invalid_orbsym_cases:
            try:
                mo_factory(
                    cs_mo_coeff if n_particle == 1 else mo_coeff_pattern(2, 0.0),
                    ao_overlap,
                    n_occ[0] if n_particle == 1 else n_occ,
                    n_particle,
                    n_ao,
                    n_mo,
                    self.mock_evaluate_dm,
                    SolverSettings(),
                    SolverState(),
                    settings,
                    invalid_orbsym,
                )
                print(
                    f" test_mo_factory_py_interface failed: {case} did not raise an "
                    "error."
                )
                test_passed = False
            except ValueError:
                pass

        self.assertTrue(test_passed, "test_mo_factory_py_interface failed")
        print(" test_mo_factory_py_interface PASSED")

    @patch.object(MOSettings, "init_c_struct", lib.mock_init_mo_settings)
    def test_mo_settings(self):
        """
        this function ensures the MOSettings object is properly initialized and
        synchronized with the underlying C struct
        """
        self.assertTrue(
            settings_hold_reference(
                "test_mo_settings", MOSettings(), lib.reference_mo_field
            ),
            "test_mo_settings failed",
        )
        print(" test_mo_settings PASSED")


@add_tests
class MOCIntegrationTests(unittest.TestCase):
    """
    this class contains integration tests that drive the real library through the
    MO C header
    """

    tests = fortran_tests["mo_c_integration_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running integration tests for MO C interface...")
        return super().setUpClass()


class MOPyIntegrationTests(unittest.TestCase):
    """
    this class contains integration tests that drive the real library (not the mock
    library) through the MO Python interface
    """

    @classmethod
    def setUpClass(cls):
        print_separator("Running integration tests for MO Python interface...")
        return super().setUpClass()

    def test_mo_factory_py(self):
        """
        this function solves the toy model through the MO factory in the closed and the
        open shell, the closed-shell MO coefficients row-major, so that the library
        rotates a copy of them, which it copies back, and the open-shell ones
        column-major, which it rotates in place, checks that the MO object is released
        once nothing keeps it alive and that the factory warns through the logging
        function of the MO settings about a missing density matrix evaluating function
        """
        test_passed = True

        # solve the toy model in the closed and the open shell
        for n_particle, case_name in ((1, "closed-shell"), (2, "open-shell")):
            run, update_orbs, obj_func = setup_mo_toy_run(n_particle)
            releases = [run.settings._extension_owner.release]
            test_passed &= run.solve(
                update_orbs, obj_func, "test_mo_factory_py", case_name
            )
            del run, update_orbs, obj_func
            test_passed &= check_released(releases, "test_mo_factory_py", case_name)

        # check that the factory warns about a missing density matrix evaluating
        # function
        run = ToyRun(True, 1)
        messages = []
        settings = MOSettings()
        settings.logger = messages.append
        settings.verbose = 3
        mo_factory(
            run.orbitals,
            toy_model().s,
            run.n_occ,
            1,
            n_ao,
            run.n_mo,
            None,
            run.settings,
            run.state,
            settings,
        )
        if evaluate_dm_missing_warning not in "".join(messages):
            print(
                " test_mo_factory_py failed: No warning about a missing density matrix "
                "evaluating function."
            )
            test_passed = False

        self.assertTrue(test_passed, "test_mo_factory_py failed")
        print(" test_mo_factory_py PASSED")

    def test_mo_settings_layout(self):
        """
        this function checks that every field of the Python MO settings structure is
        read back under its own name from settings whose fields Fortran sets one by one
        to the reference values
        """
        self.assertTrue(
            check_settings_layout(
                "test_mo_settings_layout",
                MOSettings,
                lib.get_reference_mo_values,
                lib.reference_mo_field,
            ),
            "test_mo_settings_layout failed",
        )
        print(" test_mo_settings_layout PASSED")
