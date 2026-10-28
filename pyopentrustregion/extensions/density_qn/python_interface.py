# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import annotations

import numpy as np
from ctypes import POINTER, c_bool, c_void_p, c_char, Structure, byref
from dataclasses import dataclass
from typing import TYPE_CHECKING
from pyopentrustregion.python_interface import (
    lib,
    c_int,
    c_real,
    kw_len,
    callback_type,
    CallbackInterface,
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
    ExtensionOwner,
    OrbitalOwner,
    check_callback_arguments,
    extension_functions,
)

if TYPE_CHECKING:
    from typing import Tuple, Callable, Optional, Dict

    DensityQNEvaluateDMType = Callable[
        [
            np.ndarray,
            Optional[np.ndarray],
            Optional[np.ndarray],
            Optional[np.ndarray],
            Optional[np.ndarray],
        ],
        float,
    ]


# callback function ctypes specifications, ctypes can only deal with simple return
# types so we interface to Fortran subroutines by creating pointers to the relevant
# data
density_qn_evaluate_dm_interface_type = callback_type(
    c_int,
    POINTER(c_real),
    POINTER(c_real),
    POINTER(c_real),
    POINTER(c_real),
    POINTER(c_real),
    POINTER(c_real),
)


# define classes corresponding to C structs for settings
class DensityQNSettingsC(Structure):
    _fields_ = [
        ("logger", c_void_p),
        ("initialized", c_bool),
        ("verbose", c_int),
        ("qn_method", c_char * (kw_len + 1)),
        ("context", c_void_p),
    ]


class DensityQNSettings(Settings):

    c_struct = DensityQNSettingsC
    try:
        init_c_struct = lib.init_density_qn_settings
    except AttributeError:
        raise AttributeError(
            "Please reinstall the package with: "
            "pip install -C cmake.define.ENABLE_DENSITY_QN=ON ."
        )

    callback_fields = {
        "logger": (LoggerInterface, logger_interface_type, LoggerPyInterface)
    }

    logger: Optional[Callable[[str], None]]


# ensure that appropriate fields are automatically set in settings_c object
auto_bind_fields(DensityQNSettings)
bind_callback_fields(DensityQNSettings)


# define interface factories
@dataclass
class DensityQNEvaluateDMInterface(CallbackInterface):
    """
    this class provides the interface to the density matrix evaluating function with
    Coulomb, exact-exchange and non-linear potential contributions
    """

    evaluate_dm: DensityQNEvaluateDMType
    n_ao: int
    n_particle: int
    closed_shell: bool
    exception: Dict[str, BaseException]

    def call(
        self,
        dm_ao_ptr,
        energy_ptr,
        fock_ptr,
        v_coulomb_ptr,
        v_exchange_ptr,
        v_nonlinear_ptr,
    ) -> int:
        # convert matrix pointers to numpy arrays and make the density matrix read-only
        shape = (
            2 * (self.n_ao,)
            if self.closed_shell
            else (self.n_particle, self.n_ao, self.n_ao)
        )
        dm_ao = np.ctypeslib.as_array(dm_ao_ptr, shape=shape)
        dm_ao.flags.writeable = False
        fock = np.ctypeslib.as_array(fock_ptr, shape=shape) if fock_ptr else None
        v_coulomb = (
            np.ctypeslib.as_array(v_coulomb_ptr, shape=shape) if v_coulomb_ptr else None
        )
        v_exchange = (
            np.ctypeslib.as_array(v_exchange_ptr, shape=shape)
            if v_exchange_ptr
            else None
        )
        v_nonlinear = (
            np.ctypeslib.as_array(v_nonlinear_ptr, shape=shape)
            if v_nonlinear_ptr
            else None
        )

        # get energy, and the Fock matrix, Coulomb, exact-exchange and non-linear
        # potentials where wanted
        energy_ptr[0] = self.evaluate_dm(
            dm_ao, fock, v_coulomb, v_exchange, v_nonlinear
        )

        return 0


def density_qn_factory(
    evaluate_dm: DensityQNEvaluateDMType,
    solver_settings: SolverSettings,
    solver_state: SolverState,
    settings: Optional[DensityQNSettings] = None,
) -> Tuple[
    Callable[
        [np.ndarray, np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ],
    Callable[[np.ndarray], float],
]:
    # the density QN object builds on the orbital basis object of the MO or OAO factory
    # which wired the solver settings last, whose owner holds the dimensions of the
    # orbitals the callback interfaces need, which the C interface does not receive
    orbitals = getattr(solver_settings, "_extension_owner", None)
    if not isinstance(orbitals, OrbitalOwner):
        raise ValueError(
            "The solver settings have to be those an MO or OAO factory wired last."
        )
    n_ao, n_particle = orbitals.n_ao, orbitals.n_particle

    # collector for exceptions raised inside the callback functions the object calls,
    # shared with the orbital basis object, whose functions the density QN object
    # calls, which the solver adopts from the returned functions
    exception = orbitals.exception

    # define interfaces for callback functions, whose arguments the C interface cannot
    # check
    check_callback_arguments(
        evaluate_dm,
        "evaluate_dm",
        ("dm", "fock", "v_coulomb", "v_exchange", "v_nonlinear"),
    )
    evaluate_dm_interface = density_qn_evaluate_dm_interface_type(
        DensityQNEvaluateDMInterface(
            evaluate_dm, n_ao, n_particle, n_particle == 1, exception
        )
    )

    # hand the logging function of the settings the exception collector of the call
    if settings is not None:
        settings.prepare_call(None, exception)

    # define result and argument types
    lib.density_qn_factory.restype = c_int
    lib.density_qn_factory.argtypes = [
        density_qn_evaluate_dm_interface_type,
        POINTER(update_orbs_interface_type),
        POINTER(obj_func_interface_type),
        POINTER(c_void_p),
        POINTER(SolverSettingsC),
        POINTER(SolverStateC),
        POINTER(DensityQNSettingsC),
    ]

    # call Fortran function
    update_orbs_density_qn_funptr = update_orbs_interface_type()
    obj_func_density_qn_funptr = obj_func_interface_type()
    extension = c_void_p()
    error = lib.density_qn_factory(
        evaluate_dm_interface,
        byref(update_orbs_density_qn_funptr),
        byref(obj_func_density_qn_funptr),
        byref(extension),
        byref(solver_settings.settings_c),
        byref(solver_state.state_c),
        None if settings is None else byref(settings.settings_c),
    )

    raise_on_failure(error, exception, "density QN factory")

    # the owner of the density QN object, which the solver settings keep alive, keeps
    # alive what the object calls and the owner of the orbital basis object, which the
    # settings were wired with before, and the returned functions are called with the
    # context the factory has written into the solver settings, through which they
    # reach the object, keep the solver settings alive and copy the orbitals back like
    # those of the orbital basis object
    ExtensionOwner(
        extension.value,
        solver_settings,
        {
            "settings": settings,
            "solver_state": solver_state,
            "evaluate_dm_interface": evaluate_dm_interface,
        },
    )
    return extension_functions(
        update_orbs_density_qn_funptr,
        obj_func_density_qn_funptr,
        solver_settings,
        exception,
        orbitals.copy_back,
    )
