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
    adopt_collector,
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
    ObjFuncPyInterface,
    OrbitalOwner,
    UpdateOrbsPyInterface,
    check_callback_arguments,
)

if TYPE_CHECKING:
    from typing import Tuple, Callable, Optional, Any, Union, Sequence
    from pyopentrustregion.extensions.common.python_interface import EvaluateDMType


# define classes corresponding to C structs for settings
class MOSettingsC(Structure):
    _fields_ = [
        ("logger", c_void_p),
        ("initialized", c_bool),
        ("verbose", c_int),
        ("context", c_void_p),
    ]


class MOSettings(Settings):

    c_struct = MOSettingsC
    try:
        init_c_struct = lib.init_mo_settings
    except AttributeError:
        raise AttributeError(
            "Please reinstall the package with: "
            "CMAKE_FLAGS='-DENABLE_MO=ON' pip install ."
        )

    callback_fields = {
        "logger": (LoggerInterface, logger_interface_type, LoggerPyInterface)
    }

    logger: Optional[Callable[[str], None]]


# ensure that appropriate fields are automatically set in settings_c object
auto_bind_fields(MOSettings)
bind_callback_fields(MOSettings)


def check_mo_coeff(
    mo_coeff: np.ndarray,
    ao_overlap: np.ndarray,
    n_occ: Union[int, Sequence[int]],
    n_particle: int,
    n_ao: int,
    n_mo: int,
) -> Tuple[np.ndarray, bool, Any, np.ndarray]:
    """
    this function checks that MO coefficients have the given dimensions, as an
    (n_ao, n_mo) array for the closed-shell case and as a (2, n_ao, n_mo) array for the
    open-shell case, with one number of occupied orbitals per particle channel, and
    that they can be rotated in place, and returns the column-major buffer the library
    rotates, which is the array of the caller itself if it is stored like that and a
    copy of it otherwise, whether it is the array of the caller, the occupations as a C
    array and the AO overlap matrix as a contiguous array
    """
    n_particle, n_ao, n_mo = (operator.index(n) for n in (n_particle, n_ao, n_mo))
    shape = (n_ao, n_mo) if n_particle == 1 else (n_particle, n_ao, n_mo)
    if mo_coeff.shape != shape:
        raise ValueError(
            f"The MO coefficients have to be of shape {shape} for {n_particle} "
            f"particle(s), {n_ao} AOs and {n_mo} MOs, got shape {mo_coeff.shape}."
        )
    if not np.issubdtype(mo_coeff.dtype, np.floating) or not mo_coeff.flags.writeable:
        raise ValueError(
            "The MO coefficients have to be a real floating-point and writeable array, "
            "since they are rotated in place."
        )
    n_occ_list = [operator.index(occ) for occ in np.atleast_1d(n_occ)]
    if len(n_occ_list) != n_particle:
        raise ValueError(
            "The number of occupied orbitals has to be given for every particle "
            f"channel ({n_particle}), got {len(n_occ_list)}."
        )
    n_occ_c = (c_int * n_particle)(*n_occ_list)

    # the AO overlap matrix is symmetric, so it needs no transposition
    ao_overlap = np.ascontiguousarray(ao_overlap, dtype=np.float64)
    if ao_overlap.shape != (n_ao, n_ao):
        raise ValueError(
            f"The AO overlap matrix has to be of shape ({n_ao}, {n_ao}), got shape "
            f"{ao_overlap.shape}."
        )

    # pass the MO coefficients column-major, rotating the caller's array in place if it
    # is stored like that and a copy otherwise
    mo_coeff_t = mo_coeff.swapaxes(-1, -2)
    in_place = mo_coeff.dtype == np.float64 and mo_coeff_t.flags.c_contiguous
    mo_coeff_buffer = (
        mo_coeff_t if in_place else np.ascontiguousarray(mo_coeff_t, dtype=np.float64)
    )

    return mo_coeff_buffer, in_place, n_occ_c, ao_overlap


def check_orbsym(
    orbsym: Optional[Union[Sequence[int], np.ndarray]], n_particle: int, n_mo: int
) -> Any:
    """
    this function checks that the irreps of the MOs, if given, have the given
    dimensions, as an (n_mo,) array for the closed-shell case and as an
    (n_particle, n_mo) array for the open-shell case, and returns them as a C array
    holding the irreps of every particle channel one after another, or None if they
    are not given
    """
    if orbsym is None:
        return None
    n_particle, n_mo = operator.index(n_particle), operator.index(n_mo)
    orbsym = np.asarray(orbsym)
    shape = (n_mo,) if n_particle == 1 else (n_particle, n_mo)
    if orbsym.shape != shape:
        raise ValueError(
            f"The MO irreps have to be of shape {shape} for {n_particle} particle(s) "
            f"and {n_mo} MOs, got shape {orbsym.shape}."
        )
    if not np.issubdtype(orbsym.dtype, np.integer):
        raise ValueError("The MO irreps have to be integers.")

    return (c_int * orbsym.size)(*(int(irrep) for irrep in orbsym.ravel()))


def mo_factory(
    mo_coeff: np.ndarray,
    ao_overlap: np.ndarray,
    n_occ: Union[int, Sequence[int]],
    n_particle: int,
    n_ao: int,
    n_mo: int,
    evaluate_dm: Optional[EvaluateDMType],
    solver_settings: SolverSettings,
    solver_state: SolverState,
    settings: Optional[MOSettings] = None,
    orbsym: Optional[Union[Sequence[int], np.ndarray]] = None,
) -> Tuple[
    Callable[
        [np.ndarray, np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ],
    Callable[[np.ndarray], float],
]:
    # check the MO irreps, MO coefficients and overlap matrix against the dimensions
    # and determine if closed-shell or open-shell formalism is used
    orbsym_c = check_orbsym(orbsym, n_particle, n_mo)
    mo_coeff_buffer, in_place, n_occ_c, ao_overlap = check_mo_coeff(
        mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo
    )
    closed_shell = n_particle == 1

    # get pointers to arrays
    mo_coeff_ptr = mo_coeff_buffer.ctypes.data_as(POINTER(c_real))
    ao_overlap_ptr = ao_overlap.ctypes.data_as(POINTER(c_real))

    # collector for exceptions raised inside the wrapped user callbacks; adopted from
    # evaluate_dm when it is itself factory-produced so a whole chain of factories
    # shares one, and handed to solver through the returned object
    exception = adopt_collector(evaluate_dm)

    # define interfaces for callback functions, whose arguments the C interface cannot
    # check, without a density matrix evaluating function for an object only the ARH
    # factory builds on
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

    if not hasattr(lib, "mo_factory"):
        raise RuntimeError(
            "Please reinstall the package with: "
            "CMAKE_FLAGS='-DENABLE_MO=ON' pip install ."
        )

    # define result and argument types
    lib.mo_factory.restype = c_int
    lib.mo_factory.argtypes = [
        POINTER(c_real),
        POINTER(c_real),
        POINTER(c_int),
        c_int,
        c_int,
        c_int,
        evaluate_dm_interface_type,
        POINTER(update_orbs_interface_type),
        POINTER(obj_func_interface_type),
        POINTER(c_void_p),
        POINTER(SolverSettingsC),
        POINTER(SolverStateC),
        POINTER(MOSettingsC),
        POINTER(c_int),
    ]

    # call Fortran function
    update_orbs_mo_funptr = update_orbs_interface_type()
    obj_func_mo_funptr = obj_func_interface_type()
    extension = c_void_p()
    error = lib.mo_factory(
        mo_coeff_ptr,
        ao_overlap_ptr,
        n_occ_c,
        n_particle,
        n_ao,
        n_mo,
        evaluate_dm_interface,
        byref(update_orbs_mo_funptr),
        byref(obj_func_mo_funptr),
        byref(extension),
        byref(solver_settings.settings_c),
        byref(solver_state.state_c),
        None if settings is None else byref(settings.settings_c),
        orbsym_c,
    )

    if error:
        if "exc" in exception:
            raise RuntimeError(
                f"OpenTrustRegion MO factory produced error (code {error})."
            ) from exception["exc"]
        else:
            raise RuntimeError(
                f"OpenTrustRegion MO factory produced error (code {error})."
            )

    # the owner of the MO object, which the solver settings keep alive, keeps alive
    # what the object calls and the MO coefficients it rotates and holds the dimensions
    # of the orbitals and, if the library rotates a copy of the MO coefficients, the
    # function copying them back for the ARH factory, and the returned functions are
    # called with the context the factory has written into the solver settings,
    # through which they reach the object, and keep the solver settings alive
    owner = OrbitalOwner(
        extension.value,
        solver_settings,
        {
            "settings": settings,
            "solver_state": solver_state,
            "evaluate_dm_interface": evaluate_dm_interface,
            "mo_coeff_buffer": mo_coeff_buffer,
        },
        n_ao,
        n_particle,
        (
            None
            if in_place
            else lambda: np.copyto(mo_coeff, mo_coeff_buffer.swapaxes(-1, -2))
        ),
    )
    saved_objects = {"solver_settings": solver_settings}
    return (
        UpdateOrbsPyInterface(
            update_orbs_funptr=update_orbs_mo_funptr,
            saved_objects=saved_objects,
            _otr_exception=exception,
            context=solver_settings.settings_c.context,
            copy_back=owner.copy_back,
        ),
        ObjFuncPyInterface(
            obj_func_funptr=obj_func_mo_funptr,
            saved_objects=saved_objects,
            _otr_exception=exception,
            context=solver_settings.settings_c.context,
        ),
    )
