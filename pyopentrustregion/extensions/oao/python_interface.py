# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import annotations

import operator
import numpy as np
from ctypes import POINTER, byref, c_bool, c_void_p, Structure
from typing import TYPE_CHECKING
from pyopentrustregion.python_interface import (
    lib,
    c_int,
    c_real,
    obj_func_interface_type,
    update_orbs_interface_type,
    logger_interface_type,
    LoggerInterface,
    LoggerPyInterface,
    raise_on_failure,
    Settings,
    SolverSettings,
    SolverSettingsC,
    SolverState,
    SolverStateC,
    auto_bind_fields,
    bind_callback_fields,
)
from pyopentrustregion.extensions.common.python_interface import (
    evaluate_dm_interface_type,
    EvaluateDMInterface,
    OrbitalOwner,
    check_callback_arguments,
    extension_functions,
)

if TYPE_CHECKING:
    from typing import Tuple, Callable, Optional
    from pyopentrustregion.extensions.common.python_interface import EvaluateDMType


# define classes corresponding to C structs for settings
class OAOSettingsC(Structure):
    _fields_ = [
        ("logger", c_void_p),
        ("initialized", c_bool),
        ("verbose", c_int),
        ("context", c_void_p),
    ]


class OAOSettings(Settings):

    c_struct = OAOSettingsC
    try:
        init_c_struct = lib.init_oao_settings
    except AttributeError:
        raise AttributeError(
            "Please reinstall the package with: "
            "pip install -C cmake.define.ENABLE_OAO=ON ."
        )

    callback_fields = {
        "logger": (LoggerInterface, logger_interface_type, LoggerPyInterface)
    }

    logger: Optional[Callable[[str], None]]


# ensure that appropriate fields are automatically set in settings_c object
auto_bind_fields(OAOSettings)
bind_callback_fields(OAOSettings)


def check_dm_ao(
    dm_ao: np.ndarray, ao_overlap: np.ndarray, n_particle: int, n_ao: int
) -> np.ndarray:
    """
    this function checks that an AO density matrix has the given dimensions, as an
    (n_ao, n_ao) array for the closed-shell case and as a (2, n_ao, n_ao) array for the
    open-shell case, and that it can be updated in place, since the factories keep
    pointing at it, and returns the AO overlap matrix of these dimensions as a
    contiguous array
    """
    n_particle, n_ao = operator.index(n_particle), operator.index(n_ao)
    shape = (n_ao, n_ao) if n_particle == 1 else (n_particle, n_ao, n_ao)
    if dm_ao.shape != shape:
        raise ValueError(
            f"The AO density matrix has to be of shape {shape} for {n_particle} "
            f"particle(s) and {n_ao} AOs, got shape {dm_ao.shape}."
        )
    if (
        dm_ao.dtype != np.float64
        or not dm_ao.flags.c_contiguous
        or not dm_ao.flags.writeable
    ):
        raise ValueError(
            "The AO density matrix has to be a writeable and C-contiguous float64 "
            "array, since it is updated in place."
        )

    # the AO overlap matrix is symmetric, so it needs no transposition
    ao_overlap = np.ascontiguousarray(ao_overlap, dtype=np.float64)
    if ao_overlap.shape != (n_ao, n_ao):
        raise ValueError(
            f"The AO overlap matrix has to be of shape ({n_ao}, {n_ao}), got shape "
            f"{ao_overlap.shape}."
        )

    return ao_overlap


def oao_factory(
    dm_ao: np.ndarray,
    ao_overlap: np.ndarray,
    n_particle: int,
    n_ao: int,
    evaluate_dm: Optional[EvaluateDMType],
    solver_settings: SolverSettings,
    solver_state: SolverState,
    settings: Optional[OAOSettings] = None,
) -> Tuple[
    Callable[
        [np.ndarray, np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ],
    Callable[[np.ndarray], float],
]:
    # check the density and overlap matrices against the dimensions and determine if
    # closed-shell or open-shell formalism is used
    ao_overlap = check_dm_ao(dm_ao, ao_overlap, n_particle, n_ao)
    closed_shell = n_particle == 1

    # get pointers to arrays
    dm_ao_ptr = dm_ao.ctypes.data_as(POINTER(c_real))
    ao_overlap_ptr = ao_overlap.ctypes.data_as(POINTER(c_real))

    # collector for exceptions raised inside the callback functions the object calls,
    # which the solver adopts from the returned functions and the density QN factory
    # shares
    exception = {}

    # define interfaces for callback functions, whose arguments the C interface cannot
    # check, without a density matrix evaluating function for an object only the
    # density QN factory builds on
    evaluate_dm_interface = evaluate_dm_interface_type()
    if evaluate_dm is not None:
        check_callback_arguments(
            evaluate_dm, "evaluate_dm", ("dm", "fock", "get_response")
        )
        evaluate_dm_interface = evaluate_dm_interface_type(
            EvaluateDMInterface(evaluate_dm, n_ao, n_particle, closed_shell, exception)
        )

    # hand the logging function of the settings the exception collector of the call
    if settings is not None:
        settings.prepare_call(None, exception)

    # define result and argument types
    lib.oao_factory.restype = c_int
    lib.oao_factory.argtypes = [
        POINTER(c_real),
        POINTER(c_real),
        c_int,
        c_int,
        evaluate_dm_interface_type,
        POINTER(update_orbs_interface_type),
        POINTER(obj_func_interface_type),
        POINTER(c_void_p),
        POINTER(SolverSettingsC),
        POINTER(SolverStateC),
        POINTER(OAOSettingsC),
    ]

    # call Fortran function
    update_orbs_oao_funptr = update_orbs_interface_type()
    obj_func_oao_funptr = obj_func_interface_type()
    extension = c_void_p()
    error = lib.oao_factory(
        dm_ao_ptr,
        ao_overlap_ptr,
        n_particle,
        n_ao,
        evaluate_dm_interface,
        byref(update_orbs_oao_funptr),
        byref(obj_func_oao_funptr),
        byref(extension),
        byref(solver_settings.settings_c),
        byref(solver_state.state_c),
        None if settings is None else byref(settings.settings_c),
    )

    raise_on_failure(error, exception, "OAO factory")

    # the owner of the OAO object, which the solver settings keep alive, keeps alive
    # what the object calls and the density matrix it updates in place and holds the
    # dimensions of the orbitals for the density QN factory, and the returned functions
    # are called with the context the factory has written into the solver settings,
    # through which they reach the object, and keep the solver settings alive
    OrbitalOwner(
        extension.value,
        solver_settings,
        {
            "settings": settings,
            "solver_state": solver_state,
            "evaluate_dm_interface": evaluate_dm_interface,
            "dm_ao": dm_ao,
        },
        n_ao,
        n_particle,
        exception,
    )
    return extension_functions(
        update_orbs_oao_funptr, obj_func_oao_funptr, solver_settings, exception
    )
