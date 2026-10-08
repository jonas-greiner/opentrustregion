# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import gc
import unittest
import weakref
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
    update_orbs_interface_type,
    SolverSettings,
)
from pyopentrustregion.extensions.common.tests import n_ao
from pyopentrustregion.extensions.mo import MOSettings, mo_factory, mo_deconstructor
from pyopentrustregion.extensions.mo.python_interface import UpdateOrbsMOPyInterface

if NUMPY_AVAILABLE:
    import numpy as np


# define all tests in alphabetical order
fortran_tests = {
    "mo_c_system_tests": ["mo_settings_init"],
    "mo_tests": [
        "calculate_grad_h_diag_mo",
        "count_mo_params",
        "diagonalize_per_irrep",
        "finalize_mo",
        "get_extra_trial_vectors_mo",
        "get_extra_trial_vectors_mo_callback",
        "get_hess_eigval_pairs_mo",
        "hess_x_mo_callback",
        "hess_x_static_mo",
        "init_mo_settings",
        "mo_deconstructor",
        "mo_factory_common",
        "mo_factory_cs",
        "mo_factory_os",
        "mo_param_rows",
        "mo_sanity_check",
        "mo_set_solver_settings",
        "mo_transform",
        "obj_func_mo_callback",
        "pack_ov",
        "precond_mo_callback",
        "precond_pd_mo_callback",
        "refresh_hess_eigen_mo",
        "rotate_from_hess_eigenbasis_mo",
        "rotate_mo_coeff",
        "rotate_orbitals_mo",
        "rotate_to_hess_eigenbasis_mo",
        "unpack_ov",
        "update_orbs_mo_callback",
    ],
    "mo_c_interface_tests": [
        "assign_mo_c_f",
        "assign_mo_f_c",
        "evaluate_dm_mo_f_wrapper",
        "get_extra_trial_vectors_mo_c_wrapper",
        "get_response_mo_f_wrapper",
        "hess_x_mo_c_wrapper",
        "init_mo_settings_c",
        "mo_deconstructor_c_wrapper",
        "mo_factory_c_wrapper",
        "obj_func_mo_c_wrapper",
        "precond_mo_c_wrapper",
        "precond_pd_mo_c_wrapper",
        "update_orbs_mo_c_wrapper",
    ],
}

# number of MOs and of occupied orbitals of every particle channel
n_mo = c_int.in_dll(lib, "test_n_mo").value
n_occ = list((c_int * 2).in_dll(lib, "test_n_occ"))

# irreps of the MOs of every particle channel, those of the first channel serving the
# closed-shell case
orbsym = [list(irreps) for irreps in ((c_int * n_mo) * 2).in_dll(lib, "test_orbsym")]

# multiples of the density matrix the mock density matrix evaluating function returns
# for the Fock matrix and the response
evaluate_dm_factors = list((c_real * 2).in_dll(lib, "test_evaluate_dm_factors"))


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


@unittest.skipUnless(NUMPY_AVAILABLE, "NumPy not available.")
class MOPyInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the Python interface
    """

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for MO Python interface...")

        # get fields
        fields = MOSettings.c_struct._fields_

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
        lib.get_reference_mo_values.argtypes = [POINTER(RefSettingsC)]
        lib.get_reference_mo_values.restype = None
        lib.get_reference_mo_values(byref(ref_settings))

        # extract Python values
        for name, _ in ref_fields:
            ref_value = getattr(ref_settings, name)
            if isinstance(ref_value, bytes):
                ref_value = ref_value.decode("utf-8")
            setattr(cls, name, ref_value)

        return super().setUpClass()

    assign_ref_to_settings = PyInterfaceTests.assign_ref_to_settings
    equal_settings_to_ref = PyInterfaceTests.equal_settings_to_ref

    mock_logger = PyInterfaceTests.mock_logger

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
    @patch("pyopentrustregion.python_interface.lib.mo_factory", lib.mock_mo_factory)
    def test_mo_factory_py_interface(self):
        """
        this function tests the MO factory python interface
        """
        ao_overlap = np.full(2 * (n_ao,), 2.0, dtype=np.float64)

        # initialize test flag
        test_passed = True

        # initialize settings object
        settings = MOSettings()
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
        # are rotated in place and those of a different precision are converted, with
        # the irreps of their MOs passed on only if they are given
        cases = (
            ("closed-shell row-major", 1, False, np.float64, False),
            ("closed-shell column-major", 1, True, np.float64, False),
            ("closed-shell single-precision", 1, True, np.float32, False),
            ("open-shell row-major", 2, False, np.float64, False),
            ("open-shell column-major", 2, True, np.float64, False),
            ("closed-shell with symmetry", 1, False, np.float64, True),
            ("open-shell with symmetry", 2, True, np.float64, True),
        )
        mock_passed = c_bool.in_dll(lib, "test_mo_factory_interface")
        orbsym_passed = c_bool.in_dll(lib, "test_orbsym_passed")
        for case, n_particle, column_major, dtype, with_irreps in cases:
            # initialize MO coefficients
            mo_coeff = initial_mo_coeff(n_particle, column_major, dtype)

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

            # call MO factory python interface
            obj_func_mo, update_orbs_mo = mo_factory(
                mo_coeff,
                ao_overlap,
                n_occ[0] if n_particle == 1 else n_occ[:n_particle],
                n_particle,
                n_ao,
                n_mo,
                self.mock_evaluate_dm,
                solver_settings,
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

            # determine if the MO routines are wired into the solver settings without a
            # projection
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
                    " test_mo_factory_py_interface failed: MO routines not wired "
                    f"into solver settings for the {case} case."
                )
                test_passed = False
            if (
                solver_settings.project is not caller_project
                or solver_settings.stability_settings.project is not caller_project
            ):
                print(
                    " test_mo_factory_py_interface failed: Projection wired or "
                    f"projection supplied by the caller replaced for the {case} case."
                )
                test_passed = False

            # check if the mock factory received the correct input and reset its flag
            if not mock_passed.value:
                print(
                    " test_mo_factory_py_interface failed: Mock factory received "
                    f"wrong input for the {case} case."
                )
                test_passed = False
            mock_passed.value = True

            # check if logger was called correctly
            if not self.test_logger:
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
            # rotated in place and all others copied back
            in_place = column_major and dtype == np.float64
            if (update_orbs_mo.mo_coeff is None) != in_place:
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

            # call wired MO level-shifted preconditioner function
            residual = np.full(n_param, 1.0, dtype=np.float64)
            precond_residual = np.empty(n_param, dtype=np.float64)
            mu = 5.0
            solver_settings.precond(residual, mu, precond_residual)
            if not np.allclose(precond_residual, mu):
                print(
                    " test_mo_factory_py_interface failed: Returned preconditioned "
                    "residual of wired MO level-shifted preconditioner function "
                    f"wrong for the {case} case."
                )
                test_passed = False

            # call wired MO positive-definite preconditioner function
            precond_pd_residual = np.empty(n_param, dtype=np.float64)
            solver_settings.precond_pd(residual, precond_pd_residual)
            if not np.allclose(precond_pd_residual, 3.0):
                print(
                    " test_mo_factory_py_interface failed: Returned preconditioned "
                    "residual of wired MO positive-definite preconditioner function "
                    f"wrong for the {case} case."
                )
                test_passed = False

            # call wired MO extra trial vector function
            trial_vectors = np.empty((n_extra_trial_vectors, n_param), dtype=np.float64)
            solver_settings.get_extra_trial_vectors(trial_vectors)
            if not np.allclose(
                trial_vectors,
                np.arange(1, n_extra_trial_vectors + 1, dtype=np.float64)[:, None],
            ):
                print(
                    " test_mo_factory_py_interface failed: Returned trial vectors of "
                    f"wired MO extra trial vector function wrong for the {case} case."
                )
                test_passed = False

        # call MO factory python interface for row-major MO coefficients and its
        # orbital updating function, set new starting orbitals in the same array and
        # call the factory again, which has to take them over and share the buffer the
        # library rotates with the first call, since the MO object only points to the
        # buffer of the latest call
        mo_coeff = initial_mo_coeff(1, False, np.float64)
        kappa = np.ones(n_param, dtype=np.float64)
        grad = np.empty(n_param, dtype=np.float64)
        h_diag = np.empty(n_param, dtype=np.float64)
        args = (mo_coeff, ao_overlap, n_occ[0], 1, n_ao, n_mo, self.mock_evaluate_dm)
        returned = [mo_factory(*args, SolverSettings(), settings)]
        returned[0][1](kappa, grad, h_diag)
        mo_coeff[:] = mo_coeff_pattern(1, 0.0)[0]
        returned.append(mo_factory(*args, SolverSettings(), settings))
        if not mock_passed.value:
            print(
                " test_mo_factory_py_interface failed: New starting orbitals not "
                "taken over by a repeated factory call for the same array."
            )
            test_passed = False
        mock_passed.value = True
        if returned[0][0].mo_coeff_buffer is not returned[1][0].mo_coeff_buffer:
            print(
                " test_mo_factory_py_interface failed: Buffer of the MO coefficients "
                "not shared between factory calls for the same array."
            )
            test_passed = False

        # call MO factory python interface twice for new row-major MO coefficients and
        # the orbital updating function the first call returned, whose rotated MO
        # coefficients have to reach the caller, which they would not if it copied back
        # a buffer of its own, which the library no longer rotates
        mo_coeff = initial_mo_coeff(1, False, np.float64)
        args = (mo_coeff, ao_overlap, n_occ[0], 1, n_ao, n_mo, self.mock_evaluate_dm)
        returned = [mo_factory(*args, SolverSettings(), settings) for _ in range(2)]
        returned[0][1](kappa, grad, h_diag)
        if not np.allclose(mo_coeff, mo_coeff_pattern(1, 1000.0)[0]):
            print(
                " test_mo_factory_py_interface failed: MO coefficients rotated "
                "through a later factory call not updated by the orbital updating "
                "function of an earlier one."
            )
            test_passed = False

        # call MO factory python interface for MO coefficients stored row-major and
        # column-major, release them together with the returned functions, and
        # determine if neither the array owning their memory nor the buffer the library
        # rotates is kept alive by the buffers shared between factories
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
            returned = [mo_factory(*args, SolverSettings(), settings)]
            mo_coeff_buffer_ref = weakref.ref(returned[0][0].mo_coeff_buffer)
            del mo_coeff, args, returned
            gc.collect()
            if mo_coeff_ref() is not None or mo_coeff_buffer_ref() is not None:
                print(
                    " test_mo_factory_py_interface failed: MO coefficients stored "
                    f"{'column' if column_major else 'row'}-major or their buffer "
                    "kept alive after being released together with the returned "
                    "functions."
                )
                test_passed = False

        # call an MO orbital updating function whose library routine rotates the buffer
        # of row-major MO coefficients and fails, which has to raise an error and still
        # copy the rotated MO coefficients back to the caller
        mo_coeff = initial_mo_coeff(1, False, np.float64)
        mo_coeff_buffer = np.ascontiguousarray(mo_coeff.swapaxes(-1, -2))

        def failing_update_orbs(kappa_ptr, func_ptr, grad_ptr, h_diag_ptr, hess_x_ptr):
            mo_coeff_buffer[:] = mo_coeff_pattern(1, 1000.0)[0].T
            return 1

        failing_update_orbs_funptr = update_orbs_interface_type(failing_update_orbs)
        update_orbs_mo = UpdateOrbsMOPyInterface(
            update_orbs_funptr=failing_update_orbs_funptr,
            saved_objects={},
            mo_coeff=mo_coeff,
            mo_coeff_buffer=mo_coeff_buffer,
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
                mo_factory(*args, SolverSettings(), settings)
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

    # replace original library with mock library
    @patch(
        "pyopentrustregion.python_interface.lib.mo_deconstructor",
        lib.mock_mo_deconstructor,
    )
    def test_mo_deconstructor_py_interface(self):
        """
        this function tests the MO deconstructor python interface
        """
        # initialize test flag
        test_passed = True

        # call MO deconstructor python interface
        mo_deconstructor()

        # check if deconstructor was called correctly
        if not c_bool.in_dll(lib, "test_mo_deconstructor_interface").value:
            print(
                " test_mo_deconstructor_py_interface failed: Deconstructor called "
                "wrong."
            )
            test_passed = False

        self.assertTrue(test_passed, "test_mo_deconstructor_py_interface failed")
        print(" test_mo_deconstructor_py_interface PASSED")

    @patch.object(MOSettings, "init_c_struct", lib.mock_init_mo_settings)
    def test_mo_settings(self):
        """
        this function ensure the MOSettings object is properly initialized and
        synchronized with the underlying C struct
        """
        test_passed = True
        settings = MOSettings()
        if not self.equal_settings_to_ref(settings):
            print(" test_mo_settings failed: Settings not initialized correctly.")
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
            print(" test_mo_settings failed: Optional callbacks are not set correctly.")
            test_passed = False

        self.assertTrue(test_passed, "test_mo_settings failed")
        print(" test_mo_settings PASSED")


@add_tests
class MOCSystemTests(unittest.TestCase):
    """
    this class contains system tests for the MO C interface
    """

    tests = fortran_tests["mo_c_system_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running system tests for MO C interface...")
        return super().setUpClass()
