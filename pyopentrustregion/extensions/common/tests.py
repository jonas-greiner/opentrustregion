# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import functools
import gc
import unittest
from ctypes import (
    POINTER,
    Array,
    addressof,
    byref,
    c_bool,
    c_char,
    c_void_p,
    cast,
    sizeof,
    string_at,
)
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from pyopentrustregion.tests import (
    lib,
    n_param,
    add_tests,
    print_separator,
)
from pyopentrustregion.python_interface import (
    logger_interface_type,
    c_real,
    c_int,
    hess_x_interface_type,
    update_orbs_interface_type,
    obj_func_interface_type,
    SolverSettings,
    SolverState,
    solver,
)

from pyopentrustregion.extensions.common.python_interface import (
    EvaluateDMInterface,
    ExtensionOwner,
    ObjFuncPyInterface,
    UpdateOrbsPyInterface,
    check_wrapped_functions,
    extension_functions,
)

# define all tests in alphabetical order
fortran_tests = {
    "common_tests": [
        "approx_hess_x_callback",
        "channel_rows",
        "check_evaluate_dm_orbital_basis",
        "compute_sqrt_and_inv_sqrt",
        "conv_check_callback",
        "extension_level",
        "fill_extra_trial_vectors_orbital_basis",
        "get_extra_trial_vectors_callback",
        "get_extra_trial_vectors_extension",
        "hess_x_callback",
        "hess_x_extension",
        "level_shifted_divisors",
        "logger_callback",
        "matrix_exponential",
        "obj_func_callback",
        "obj_func_extension",
        "positive_definite_divisors",
        "precond_callback",
        "precond_extension",
        "precond_orbital_basis",
        "precond_pd_callback",
        "precond_pd_extension",
        "precond_pd_orbital_basis",
        "project_callback",
        "project_extension",
        "refresh_response_orbital_basis",
        "stability_hess_x_callback",
        "update_orbs_callback",
        "update_orbs_extension",
        "wire_extension_solver_settings",
        "wire_solver_settings_extension",
    ],
    "common_c_interface_tests": [
        "approx_hess_x_c_wrapper",
        "assign_solver_state_f_c",
        "assign_stability_state_f_c",
        "bundle_settings_c",
        "conv_check_c_wrapper",
        "evaluate_dm_cs_f_wrapper",
        "evaluate_dm_os_f_wrapper",
        "finish_extension_holder",
        "get_extra_trial_vectors_c_wrapper",
        "get_response_cs_f_wrapper",
        "get_response_os_f_wrapper",
        "hess_x_c_wrapper",
        "holder_level",
        "logger_c_wrapper",
        "obj_func_c_wrapper",
        "precond_c_wrapper",
        "precond_pd_c_wrapper",
        "project_c_wrapper",
        "refresh_state",
        "release_extension",
        "setup_extension_holder",
        "stability_hess_x_c_wrapper",
        "update_orbs_c_wrapper",
        "wire_solver_settings_c",
    ],
}

# the shared extension routines are only built together with an extension which uses
# them, so raise an AttributeError on import otherwise
getattr(lib, "test_n_ao")

# number of AOs
n_ao = c_int.in_dll(lib, "test_n_ao").value

# warning the factories of the orbital bases print without a density matrix evaluating
# function
evaluate_dm_missing_warning = string_at(
    addressof(c_char.in_dll(lib, "test_evaluate_dm_missing_warning"))
).decode()


class ExtensionCallbackMocks:
    """
    this class collects the mock callback functions the Python interface tests of the
    extensions share
    """

    def mock_hess_x(self, x, hess_x):
        """
        this function is a mock function for the Hessian linear transformation
        function
        """
        hess_x[:] = 4 * x

    def mock_update_orbs(self, kappa, grad, h_diag):
        """
        this function is a mock function for the orbital update function
        """
        func = np.sum(kappa)
        grad[:] = 2 * kappa
        h_diag[:] = 3 * kappa

        return func, self.mock_hess_x

    def mock_obj_func(self, kappa):
        """
        this function is a mock function for the objective function
        """
        return np.sum(kappa)


# multiples of the density matrix the mock density matrix evaluating function returns
# for the Fock matrix and the response
evaluate_dm_factors = list((c_real * 2).in_dll(lib, "test_evaluate_dm_factors"))


class OrbitalCallbackMocks:
    """
    this class collects the mock exact-evaluation callback functions the Python
    interface tests of the orbital bases share
    """

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


def ao_overlap_pattern():
    """
    this function returns an AO overlap matrix whose values encode the pair of its AO
    indices regardless of their order, as ao_overlap_pattern of the Fortran test
    reference does, so that it is symmetric while its values reveal a matrix of the
    wrong shape
    """
    ao_index = np.arange(1, n_ao + 1)
    return (
        10 * np.maximum.outer(ao_index, ao_index) + np.minimum.outer(ao_index, ao_index)
    ).astype(np.float64)


def release_mock_extension(extension):
    """
    this function frees an extension a mock factory returned through the mock release
    of the library, so that the owners the Python factories set up in their tests,
    whose finalizers keep the release they were set up with, never call the release of
    the library, also when they are finalized after the test
    """
    lib.mock_release_extension.restype = None
    lib.mock_release_extension.argtypes = [c_void_p]
    lib.mock_release_extension(extension)


def extension_reference(lookup, name, field_type):
    """
    this function returns the reference value of a field of the settings of an
    extension by its name from the lookup function of the Fortran test reference module
    of the extension
    """
    value = c_real()
    keyword = dict(SolverSettings.c_struct._fields_)["subsystem_solver"]()
    lookup(name.encode(), byref(value), keyword)
    if issubclass(field_type, Array):
        return keyword.value.decode()
    if field_type == c_bool:
        return bool(value.value)
    if field_type == c_int:
        return int(value.value)
    return value.value


def assign_reference_settings(settings, lookup):
    """
    this function assigns the reference values, which the lookup function of the
    Fortran test reference module of the extension returns by name, to every field of
    the settings of an extension but the logging function, the host context and the
    initialization flag
    """
    for field_name, field_type in settings.c_struct._fields_:
        if field_type == c_void_p or field_name == "initialized":
            continue
        setattr(
            settings, field_name, extension_reference(lookup, field_name, field_type)
        )


def setup_wired_solver_settings():
    """
    this function returns solver settings for the Python test of a factory that wires
    them, with a projection and logging function and a convergence check of their
    stability check, whose C callback functions the factory has to be handed, together
    with these Python functions
    """
    solver_settings = SolverSettings()
    solver_settings.project = lambda vector: None
    solver_settings.logger = lambda message: None
    solver_settings.stability_settings.conv_check = lambda: False
    python_functions = (
        solver_settings.project,
        solver_settings.logger,
        solver_settings.stability_settings.conv_check,
    )

    return solver_settings, python_functions


def check_wired_solver_settings(
    solver_settings, python_functions, test_name, returned, extension_settings
):
    """
    this function checks that solver settings a mock factory wired still hold their
    Python functions, that the preconditioners the mock factory wired into both levels
    are called with the contexts of their levels, that a Python function assigned
    after the wiring is called with the context the mock factory wrote into the
    settings, and that the functions the factory returned keep the solver settings and
    the settings of the extension alive, whose callback functions its object calls
    """
    test_passed = True

    if (
        solver_settings.project,
        solver_settings.logger,
        solver_settings.stability_settings.conv_check,
    ) != python_functions:
        print(f" {test_name} failed: Python functions of the settings not kept.")
        test_passed = False

    residual = np.full(n_param, 1.0, dtype=np.float64)
    precond_residual = np.empty(n_param, dtype=np.float64)
    for level, settings in (
        ("solver", solver_settings),
        ("stability check", solver_settings.stability_settings),
    ):
        precond_residual[:] = 0.0
        try:
            settings.precond(residual, 2.0, precond_residual)
        except RuntimeError:
            print(
                f" {test_name} failed: Wired preconditioner of the {level} "
                "raises error."
            )
            test_passed = False
            continue
        if not np.allclose(precond_residual, 2.0):
            print(
                f" {test_name} failed: Wired preconditioner of the {level} "
                "returns wrong result."
            )
            test_passed = False

    # assign a Hessian linear transformation for the stability check after the wiring
    # and call it through the C settings with the context the mock factory wrote and
    # the number of parameters and exception collector the solver would hand it
    x = np.ones(n_param, dtype=np.float64)
    hess_x = np.zeros(n_param, dtype=np.float64)

    def stability_hess_x(x, hess_x):
        hess_x[:] = 4 * x

    solver_settings.stability_hess_x = stability_hess_x
    solver_settings.prepare_call(n_param, {})
    error = hess_x_interface_type(solver_settings.settings_c.stability_hess_x)(
        x.ctypes.data_as(POINTER(c_real)),
        hess_x.ctypes.data_as(POINTER(c_real)),
        solver_settings.settings_c.context,
    )
    if error != 0 or not np.allclose(hess_x, 4.0):
        print(
            f" {test_name} failed: Python function assigned after the wiring "
            "not called."
        )
        test_passed = False
    solver_settings.stability_hess_x = None

    # the returned functions keep the solver settings alive and the extension's
    # settings directly or through the owner of the extension, which the solver
    # settings keep alive
    owner = getattr(solver_settings, "_extension_owner", None)
    for function in returned:
        kept_settings = function.saved_objects.get("settings")
        if kept_settings is None and owner is not None:
            kept_settings = owner.saved_objects.get("settings")
        if (
            function.saved_objects.get("solver_settings") is not solver_settings
            or kept_settings is not extension_settings
        ):
            print(
                f" {test_name} failed: Returned function does not keep the "
                "settings alive."
            )
            test_passed = False

    return test_passed


def settings_hold_reference(test_name, settings, lookup):
    """
    this function checks that the logging function of the settings of an extension is
    unset and that every other field but the host context holds the reference value the
    lookup function of the Fortran test reference module of the extension returns by
    name, and that assigning a logging function writes the C callback function the
    settings built for it into the C settings, which calls it whatever context it
    receives, while reading returns the logging function
    """
    test_passed = True
    for field_name, field_type in settings.c_struct._fields_:
        if field_name == "context":
            continue
        if field_type == c_void_p:
            if getattr(settings, field_name) is not None:
                print(f" {test_name} failed: Function {field_name} not unset.")
                test_passed = False
            continue
        value = getattr(settings, field_name)
        ref_value = extension_reference(lookup, field_name, field_type)
        if field_type == c_real:
            match = np.isclose(value, ref_value)
        else:
            match = value == ref_value
        if not match:
            print(f" {test_name} failed: Field {field_name} not initialized correctly.")
            test_passed = False

    # assign a logging function and call the C callback function the settings built
    messages = []
    settings.logger = messages.append
    logger_interface_type(settings.settings_c.logger)(b"test", 12345)
    if settings.logger != messages.append or messages != ["test"]:
        print(f" {test_name} failed: Logging function not set correctly.")
        test_passed = False

    return test_passed


def check_settings_layout(test_name, settings, export, lookup, size):
    """
    this function checks that the ctypes structure of the settings of an extension has
    the size of the bind(C) settings type and that every field is read back under its
    own name from C settings whose fields Fortran sets one by one to the reference
    values, compared against the value the lookup function returns for the field's
    name, so that a stray field is caught as well as a misplaced one
    """
    if sizeof(settings.c_struct) != size:
        print(
            f" {test_name} failed: Size of {settings.__name__} structure differs from "
            "the library's."
        )
        return False
    settings_c = settings.c_struct()
    export(byref(settings_c))
    test_passed = True
    for name, field_type in settings.c_struct._fields_:
        value = getattr(settings_c, name)
        if field_type == c_void_p:
            value = value or 0
        elif issubclass(field_type, Array):
            value = value.decode()
        if value != extension_reference(lookup, name, field_type):
            print(f" {test_name} failed: Field {name} misplaced.")
            test_passed = False
    return test_passed


def check_extension_owner(call_factory, test_name):
    """
    this function checks that a mock factory, which call_factory calls with solver
    settings and a solver state, received the record of the state and handed the
    extension it returned to an owner, which the solver settings keep alive and which
    keeps the state and exactly the callback interfaces whose C functions the factory
    received alive
    """
    test_passed = True
    received_solver_state = c_void_p.in_dll(lib, "received_solver_state")
    returned_extension = c_void_p.in_dll(lib, "returned_extension")
    received_callbacks = (c_void_p * 4).in_dll(lib, "received_callbacks")

    # call the factory
    solver_settings = setup_wired_solver_settings()[0]
    solver_state = SolverState()
    call_factory(solver_settings, solver_state)

    # the factory received the record of the state, and the solver settings keep the
    # owner of the extension it returned, which keeps the state alive
    if received_solver_state.value != addressof(solver_state.state_c):
        print(f" {test_name} failed: Record of the solver state not passed.")
        test_passed = False
    owner = getattr(solver_settings, "_extension_owner", None)
    if (
        owner is None
        or owner.extension != returned_extension.value
        or owner.saved_objects.get("solver_state") is not solver_state
    ):
        print(f" {test_name} failed: Returned extension not owned.")
        test_passed = False
    elif {
        cast(interface, c_void_p).value
        for key, interface in owner.saved_objects.items()
        if key.endswith("_interface") and interface is not None
    } != {callback for callback in received_callbacks if callback}:
        print(
            f" {test_name} failed: Callback interfaces whose C functions the "
            "factory received not kept alive."
        )
        test_passed = False

    return test_passed


@functools.cache
def toy_model():
    """
    this function reads the toy model of common_test_reference.f90, which the
    integration tests of the orbital bases and density QN solve: a mean-field energy in
    a non-orthogonal AO basis with the core Hamiltonian, the AO overlap matrix, the
    two-electron kernel g, whose Coulomb and exchange matrices are J[D]_pq = sum_rs
    g_pqrs D_rs and K[D]_pq = sum_rs g_prqs D_rs, a non-linear term c sum_p n_p^3 of
    the total diagonal density n_p, starting orbitals orthonormal in the AO metric and
    the occupations of both particle channels
    """
    h, s, mo_coeff = (np.empty((n_ao, n_ao), order="F") for _ in range(3))
    g = np.empty(4 * (n_ao,), order="F")
    nonlinear_factor = c_real()
    n_occ = (c_int * 2)()
    lib.get_toy_model(
        *(array.ctypes.data_as(POINTER(c_real)) for array in (h, s, g)),
        byref(nonlinear_factor),
        mo_coeff.ctypes.data_as(POINTER(c_real)),
        n_occ,
    )
    return SimpleNamespace(
        h=h, s=s, g=g, c=nonlinear_factor.value, mo_coeff=mo_coeff, n_occ=list(n_occ)
    )


class ToyRun:
    """
    this class prepares a solve of the toy model through the Python interface in the
    given orbital basis and shell, with solver settings requesting a stability check
    and holding a projection of the host, which counts its calls, the solver state and
    the orbitals the orbital basis moves, MO coefficients, the closed-shell ones
    row-major, so that the library rotates a copy of them, and the open-shell ones
    column-major, so that it rotates them in place, or density matrices, together with
    the density matrix evaluating functions of the orbital basis and density QN, which
    count their calls and record whether the solver state was current in them
    """

    def __init__(self, mo_basis, n_particle):
        model = toy_model()
        self.mo_basis = mo_basis
        self.n_particle = n_particle
        self.n_occ = model.n_occ[:n_particle]
        self.settings = SolverSettings()
        self.settings.stability = True
        self.settings.project = self.project
        self.state = SolverState()
        self.calls = {
            "evaluate_dm": 0,
            "get_response": 0,
            "density_qn_evaluate_dm": 0,
            "project": 0,
        }
        self.n_updates_seen = {"evaluate_dm": 0, "density_qn_evaluate_dm": 0}
        self.state_not_current = False

        # MO coefficients of fewer MOs than AOs for every channel, so that the numbers
        # of AOs and MOs cannot be confused, or the density matrices of their occupied
        # orbitals for every channel
        if mo_basis:
            self.n_mo = n_ao - 1
            mo_coeff = model.mo_coeff[:, : self.n_mo]
            self.orbitals = (
                np.ascontiguousarray(mo_coeff)
                if n_particle == 1
                else np.ascontiguousarray(np.array(2 * [mo_coeff.T])).swapaxes(-1, -2)
            )
            self.n_param = sum(occ * (self.n_mo - occ) for occ in self.n_occ)
        else:
            dm = np.array(
                [
                    model.mo_coeff[:, :occ] @ model.mo_coeff[:, :occ].T
                    for occ in self.n_occ
                ]
            )
            self.orbitals = np.ascontiguousarray(dm[0] if n_particle == 1 else dm)
            self.n_param = n_particle * n_ao * (n_ao - 1) // 2
        self.start_orbitals = self.orbitals.copy()

    def called(self, name):
        """
        this function counts a call of a density matrix evaluating or response function
        and records whether the number of orbital updates and objective function
        evaluations in the solver state ever decreased between the calls of a density
        matrix evaluating function, which it does once the state is not kept current
        """
        self.calls[name] += 1
        if name in self.n_updates_seen:
            n_updates = self.state.n_update_orbs + self.state.n_obj_func
            if n_updates < self.n_updates_seen[name]:
                self.state_not_current = True
            self.n_updates_seen[name] = max(self.n_updates_seen[name], n_updates)

    @staticmethod
    def outputs(dm, fock=None, v_coulomb=None, v_exchange=None, v_nonlinear=None):
        """
        this function returns the energy and writes, where wanted, the Fock matrix and
        the Coulomb, exact-exchange and non-linear potentials of the closed-shell
        density matrix of one spin channel or the open-shell density matrices of both
        channels
        """
        model = toy_model()
        closed_shell = dm.ndim == 2
        dms = dm.reshape(-1, n_ao, n_ao)
        dm_total = 2.0 * dm if closed_shell else dm.sum(axis=0)
        n_diag = np.diag(dm_total)
        j = np.einsum("pqrs,irs->ipq", model.g, dms)
        k = np.einsum("prqs,irs->ipq", model.g, dms)
        j_total = np.einsum("pqrs,rs->pq", model.g, dm_total)
        nonlinear = np.diag(3.0 * model.c * n_diag**2)
        for output, value in (
            (fock, model.h + j_total - k + nonlinear),
            (v_coulomb, j_total if closed_shell else j),
            (v_exchange, k),
            (v_nonlinear, nonlinear),
        ):
            if output is not None:
                output[...] = np.broadcast_to(value, dms.shape).reshape(output.shape)

        # closed shell: E = 2 tr(h D) + tr(D J[2 D]) - tr(D K[D]) + c sum_p n_p^3, open
        # shell: E = sum_s tr(h D_s) + tr(D_s J[D_a + D_b]) / 2 - tr(D_s K[D_s]) / 2 +
        # c sum_p n_p^3
        weight = 2.0 if closed_shell else 1.0
        return weight * np.sum(
            dms * (model.h + 0.5 * j_total - 0.5 * k)
        ) + model.c * np.sum(n_diag**3)

    def evaluate_dm(self, dm, fock, get_response):
        """
        this function is the density matrix evaluating function of the orbital basis,
        whose response function belongs to the density matrices it was requested for,
        which the library may overwrite after the call
        """
        self.called("evaluate_dm")
        energy = self.outputs(dm, fock=fock)
        if not get_response:
            return energy, None
        model = toy_model()
        dm_total = 2.0 * dm if dm.ndim == 2 else dm.sum(axis=0)
        n_diag = np.diag(dm_total)

        def response(ddm, response):
            self.called("get_response")
            ddm_total = 2.0 * ddm if ddm.ndim == 2 else ddm.sum(axis=0)
            response[...] = (
                np.einsum("pqrs,rs->pq", model.g, ddm_total)
                - np.einsum("prqs,...rs->...pq", model.g, ddm)
                + np.diag(6.0 * model.c * n_diag * np.diag(ddm_total))
            )

        return energy, response

    def project(self, vector):
        """
        this function is the projection of the host, the identity, which counts its
        calls, so that the orbital bases are seen to apply it, the OAO object after its
        own
        """
        self.called("project")

    def density_qn_evaluate_dm(self, dm, fock, v_coulomb, v_exchange, v_nonlinear):
        """
        this function is the density matrix evaluating function of density QN
        """
        self.called("density_qn_evaluate_dm")
        return self.outputs(dm, fock, v_coulomb, v_exchange, v_nonlinear)

    def solve(self, update_orbs, obj_func, test_name, case_name, density_qn=False):
        """
        this function solves the prepared run with the functions a factory returned and
        checks that the solve converged, that the density matrix evaluating functions
        were called with the current solver state, that the exact Hessian of the
        orbital basis and, with density QN, its density matrix evaluating function were
        used, that the projection of the host was applied and that the orbitals of the
        caller moved and stayed orthonormal
        """
        try:
            solver(update_orbs, obj_func, self.n_param, self.settings, self.state)
        except RuntimeError as e:
            print(f" {test_name} failed: Solver raised for the {case_name} case: {e}")
            return False
        test_passed = True
        if not (
            self.state.grad_rms < self.settings.conv_tol
            or self.state.max_precision_reached
        ):
            print(
                f" {test_name} failed: Solver did not converge for the {case_name} case."
            )
            test_passed = False

        # the exact Hessian of the orbital basis is used for the solve or, with density
        # QN, for the stability check, so the response function is called, and the
        # projection of the host is applied, also behind the OAO projection
        names = ["evaluate_dm", "get_response", "project"] + (
            ["density_qn_evaluate_dm"] if density_qn else []
        )
        for name in names:
            if self.calls[name] == 0:
                print(
                    f" {test_name} failed: {name} not called for the {case_name} case."
                )
                test_passed = False
        n_updates = self.state.n_update_orbs + self.state.n_obj_func
        if self.state_not_current or not all(
            0 < self.n_updates_seen[name] <= n_updates
            for name in names
            if name in self.n_updates_seen
        ):
            print(
                f" {test_name} failed: Solver state not current in the density matrix "
                f"evaluating functions for the {case_name} case."
            )
            test_passed = False

        # the orbitals of the caller moved and stayed orthonormal: C^T S C = I in the MO
        # basis, D S D = D in the OAO basis
        if np.max(np.abs(self.orbitals - self.start_orbitals)) < 1e-3:
            print(
                f" {test_name} failed: Orbitals of the caller not moved for the "
                f"{case_name} case."
            )
            test_passed = False
        s = toy_model().s
        for x in self.orbitals.reshape((self.n_particle,) + self.orbitals.shape[-2:]):
            target = np.eye(x.shape[1]) if self.mo_basis else x
            if np.max(np.abs(x.T @ s @ x - target)) > 1e-8:
                print(
                    f" {test_name} failed: Orbitals of the caller not orthonormal for "
                    f"the {case_name} case."
                )
                test_passed = False
        return test_passed


def check_released(releases, test_name, case_name):
    """
    this function checks that the extension objects whose owners have the given
    finalizers were released once nothing keeps them alive any longer
    """
    gc.collect()
    if any(release.alive for release in releases):
        print(f" {test_name} failed: Extension not released for the {case_name} case.")
        return False
    return True


@add_tests
class CommonTests(unittest.TestCase):
    """
    this class contains unit tests for the routines shared between extensions
    """

    tests = fortran_tests["common_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for shared extension routines...")
        return super().setUpClass()


@add_tests
class CommonCInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the C interface shared between extensions
    """

    tests = fortran_tests["common_c_interface_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for shared extension C interface...")
        return super().setUpClass()


class CommonPyInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the Python interface shared between extensions
    """

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for shared extension Python interface...")
        return super().setUpClass()

    def test_check_wrapped_functions(self):
        """
        this function tests the function which checks that a factory wraps the orbital
        updating function and the objective function the factory that wired the solver
        settings last returned, or the host's own for settings no factory wired
        """
        # initialize test flag
        test_passed = True

        # host's functions and functions returned by factories for the levels at two
        # addresses
        def host_update_orbs(kappa, grad, h_diag):
            return 0.0, None

        def host_obj_func(kappa):
            return 0.0

        def marked_wrapper(function):
            @functools.wraps(function)
            def wrapper(*args):
                return function(*args)

            return wrapper

        level, other_level = 0x1000, 0x2000
        returned_update_orbs = {
            address: UpdateOrbsPyInterface(
                update_orbs_funptr=None, saved_objects={}, context=address
            )
            for address in (level, other_level)
        }
        returned_obj_func = {
            address: ObjFuncPyInterface(
                obj_func_funptr=None, saved_objects={}, context=address
            )
            for address in (level, other_level)
        }

        # cases of the context of the settings, None for settings no factory wired,
        # the functions passed and whether they have to be accepted
        cases = (
            (
                "host's functions for settings no factory wired",
                None,
                host_update_orbs,
                host_obj_func,
                True,
            ),
            (
                "returned functions for settings no factory wired",
                None,
                returned_update_orbs[level],
                returned_obj_func[level],
                False,
            ),
            (
                "a returned orbital updating function for settings no factory wired",
                None,
                returned_update_orbs[level],
                host_obj_func,
                False,
            ),
            (
                "a returned objective function for settings no factory wired",
                None,
                host_update_orbs,
                returned_obj_func[level],
                False,
            ),
            (
                "functions returned for the wired settings",
                level,
                returned_update_orbs[level],
                returned_obj_func[level],
                True,
            ),
            (
                "host's functions for wired settings",
                level,
                host_update_orbs,
                host_obj_func,
                False,
            ),
            (
                "the host's orbital updating function for wired settings",
                level,
                host_update_orbs,
                returned_obj_func[level],
                False,
            ),
            (
                "the host's objective function for wired settings",
                level,
                returned_update_orbs[level],
                host_obj_func,
                False,
            ),
            (
                "functions returned for settings wired before",
                level,
                returned_update_orbs[other_level],
                returned_obj_func[other_level],
                False,
            ),
            (
                "wrappers made with functools.wraps of the functions returned for the "
                "wired settings",
                level,
                marked_wrapper(returned_update_orbs[level]),
                marked_wrapper(returned_obj_func[level]),
                True,
            ),
            (
                "wrappers made with functools.wraps of functions returned for settings "
                "wired before",
                level,
                marked_wrapper(returned_update_orbs[other_level]),
                marked_wrapper(returned_obj_func[other_level]),
                False,
            ),
            (
                "unmarked wrappers of the functions returned for the wired settings",
                level,
                lambda *args: returned_update_orbs[level](*args),
                lambda *args: returned_obj_func[level](*args),
                False,
            ),
        )
        for case, context, update_orbs, obj_func, accepted in cases:
            solver_settings = SolverSettings()
            solver_settings.settings_c.context = context
            try:
                check_wrapped_functions(update_orbs, obj_func, solver_settings)
                raised = False
            except ValueError:
                raised = True
            if raised == accepted:
                print(
                    f" test_check_wrapped_functions failed: {case} "
                    f"{'refused' if accepted else 'accepted'}."
                )
                test_passed = False

        self.assertTrue(test_passed, "test_check_wrapped_functions failed")
        print(" test_check_wrapped_functions PASSED")

    def test_evaluate_dm_py_interface(self):
        """
        this function tests the Python interface to the density matrix evaluating
        function
        """
        # initialize test flag
        test_passed = True

        # check if a density matrix evaluating function which returns no response
        # function although one is requested produces an error
        def evaluate_dm_without_response(dm_ao, fock, get_response):
            return np.sum(dm_ao), None

        exception = {}
        evaluate_dm = EvaluateDMInterface(
            evaluate_dm_without_response, n_ao, 1, True, exception
        )
        dm_ao = np.full(2 * (n_ao,), 1.0, dtype=np.float64)
        energy = (c_real * 1)()
        get_response_funptr = (c_void_p * 1)()
        error = evaluate_dm(
            dm_ao.ctypes.data_as(POINTER(c_real)),
            energy,
            None,
            get_response_funptr,
            None,
        )
        if error == 0 or not isinstance(exception.get("exc"), RuntimeError):
            print(
                " test_evaluate_dm_py_interface failed: Missing response function "
                "requested from density matrix evaluating function does not produce "
                "an error."
            )
            test_passed = False

        self.assertTrue(test_passed, "test_evaluate_dm_py_interface failed")
        print(" test_evaluate_dm_py_interface PASSED")

    def test_extension_functions(self):
        """
        this function tests the orbital updating function and the objective function a
        factory returns around the C functions the C factory returned, which have to
        call them with the context of the solver settings, keep these alive, share the
        collector of the factory call, copy the orbitals back after every orbital update
        and hand over and receive arrays of any layout and kind
        """
        test_passed = True
        solver_settings = SolverSettings()
        exception = {}
        copies = []
        update_orbs, obj_func = extension_functions(
            update_orbs_interface_type(lib.mock_update_orbs),
            obj_func_interface_type(lib.mock_obj_func),
            solver_settings,
            exception,
            lambda: copies.append(True),
        )

        for function in (update_orbs, obj_func):
            if (
                function.saved_objects["solver_settings"] is not solver_settings
                or function._otr_exception is not exception
                or function.context != solver_settings.settings_c.context
            ):
                print(
                    " test_extension_functions failed: Returned function does not "
                    "keep the solver settings alive, share the collector or call with "
                    "their context."
                )
                test_passed = False

        # the mock orbital update returns the sum, twice and three times the parameters
        # and a Hessian linear transformation applying four times the identity, the
        # mock objective function the sum, for a non-contiguous view of the parameters,
        # a float32 gradient and a non-contiguous Hessian diagonal
        kappa = np.random.random((n_param, 2))[:, 1]
        grad = np.zeros(n_param, dtype=np.float32)
        h_diag = np.zeros((n_param, 2))[:, 0]
        func, hess_x = update_orbs(kappa, grad, h_diag)
        if not (
            np.isclose(func, kappa.sum())
            and np.allclose(grad, 2.0 * kappa, rtol=1e-6)
            and np.allclose(h_diag, 3.0 * kappa)
        ):
            print(
                " test_extension_functions failed: Arrays of the orbital update not "
                "handed over or received."
            )
            test_passed = False
        if copies != [True]:
            print(
                " test_extension_functions failed: Orbitals not copied back after the "
                "orbital update."
            )
            test_passed = False

        # the returned Hessian linear transformation keeps the orbital updating function
        # alive, which keeps what the object of the factory calls alive
        if hess_x._owner is not update_orbs:
            print(
                " test_extension_functions failed: Returned Hessian linear "
                "transformation does not keep the orbital updating function alive."
            )
            test_passed = False
        hess_x_out = np.zeros(n_param, dtype=np.float32)
        hess_x(kappa, hess_x_out)
        if not np.allclose(hess_x_out, 4.0 * kappa, rtol=1e-6):
            print(
                " test_extension_functions failed: Arrays of the Hessian linear "
                "transformation not handed over or received."
            )
            test_passed = False
        if not np.isclose(obj_func(kappa), kappa.sum()):
            print(
                " test_extension_functions failed: Parameters of the objective "
                "function not handed over."
            )
            test_passed = False

        self.assertTrue(test_passed, "test_extension_functions failed")
        print(" test_extension_functions PASSED")

    @patch(
        "pyopentrustregion.python_interface.lib.release_extension",
        lib.mock_release_extension,
    )
    def test_extension_owner(self):
        """
        this function tests the owner of the extension object a factory returned for
        the solver settings it wired
        """
        # initialize test flag
        test_passed = True
        released_extension = c_void_p.in_dll(lib, "released_extension")
        n_released_extensions = c_int.in_dll(lib, "n_released_extensions")
        lib.mock_extension.restype = c_void_p
        lib.mock_extension.argtypes = [c_void_p]

        # let two extensions be owned one after the other for the same solver settings
        solver_settings = SolverSettings()
        kept = object()
        extensions = (lib.mock_extension(None), lib.mock_extension(None))
        first = ExtensionOwner(extensions[0], solver_settings, {"kept": kept})
        second = ExtensionOwner(extensions[1], solver_settings, {})

        # the settings keep the latest owner alive, which keeps the earlier one alive,
        # which keeps what it was given alive
        if (
            getattr(solver_settings, "_extension_owner", None) is not second
            or second.saved_objects.get("wired_before") is not first
            or first.saved_objects.get("wired_before") is not None
            or first.saved_objects.get("kept") is not kept
            or first.extension != extensions[0]
        ):
            print(" test_extension_owner failed: Owners not kept alive.")
            test_passed = False

        # dropping the earlier owner releases nothing, since the latest keeps it alive,
        # after releasing the owners earlier tests left in reference cycles
        gc.collect()
        n_released_extensions.value = 0
        del first
        gc.collect()
        if n_released_extensions.value != 0:
            print(" test_extension_owner failed: Owner wired over released too early.")
            test_passed = False

        # dropping the settings and the latest owner releases both, the latest first
        del second, solver_settings
        gc.collect()
        if (
            n_released_extensions.value != 2
            or released_extension.value != extensions[0]
        ):
            print(
                " test_extension_owner failed: Extensions not released in reverse order "
                "of wiring."
            )
            test_passed = False

        self.assertTrue(test_passed, "test_extension_owner failed")
        print(" test_extension_owner PASSED")
