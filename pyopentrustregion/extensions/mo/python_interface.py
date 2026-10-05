# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import annotations

import operator
import weakref
import numpy as np
from ctypes import POINTER, byref, c_bool, c_void_p, Structure
from dataclasses import dataclass
from typing import TYPE_CHECKING
from pyopentrustregion.python_interface import (
    lib,
    c_int,
    c_real,
    obj_func_interface_type,
    update_orbs_interface_type,
    logger_interface_type,
    LoggerInterface,
    adopt_collector,
    Settings,
    SolverSettings,
    SolverSettingsC,
    auto_bind_fields,
)
from pyopentrustregion.extensions.common.python_interface import (
    evaluate_dm_interface_type,
    EvaluateDMInterface,
    ObjFuncPyInterface,
    UpdateOrbsPyInterface,
    attach_wired_callbacks,
)

if TYPE_CHECKING:
    from typing import Tuple, Callable, Optional, Any, Union, Sequence, Dict
    from pyopentrustregion.extensions.common.python_interface import EvaluateDMType


# define classes corresponding to C structs for settings
class MOSettingsC(Structure):
    _fields_ = [
        ("logger", c_void_p),
        ("initialized", c_bool),
        ("verbose", c_int),
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

    logger: Optional[Callable[[str], None]]
    logger_interface: Any


# ensure that appropriate fields are automatically set in settings_c object
auto_bind_fields(MOSettings)


# column-major copies of the MO coefficient arrays of the callers which cannot be
# rotated in place, which the factories of every MO-basis extension share per array,
# keyed by the identity of the array of the caller, which is only referenced weakly
_mo_coeff_buffers: Dict[int, Tuple[weakref.ref, np.ndarray]] = {}


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
    rotates, whether it is the array of the caller itself, the occupations as a C array
    and the AO overlap matrix as a contiguous array; the buffer is shared by every
    factory called with the same array object, since all of them point the same MO
    object at it, and takes over the current values of the array of the caller
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
    # is stored like that, which every factory then does through its own view, and a
    # buffer copied back after every update otherwise, taking over the buffer of a
    # previous factory called with the same array
    mo_coeff_t = mo_coeff.swapaxes(-1, -2)
    in_place = mo_coeff.dtype == np.float64 and mo_coeff_t.flags.c_contiguous
    if in_place:
        mo_coeff_buffer = mo_coeff_t
    else:
        entry = _mo_coeff_buffers.get(id(mo_coeff))
        if (
            entry is not None
            and entry[0]() is mo_coeff
            and entry[1].shape == mo_coeff_t.shape
        ):
            mo_coeff_buffer = entry[1]
            np.copyto(mo_coeff_buffer, mo_coeff_t)
        else:
            mo_coeff_buffer = np.ascontiguousarray(mo_coeff_t, dtype=np.float64)
            key = id(mo_coeff)
            _mo_coeff_buffers[key] = (
                weakref.ref(mo_coeff, lambda _: _mo_coeff_buffers.pop(key, None)),
                mo_coeff_buffer,
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


# define interface factories
@dataclass
class ObjFuncMOPyInterface(ObjFuncPyInterface):
    """
    this class provides the Python interface to the objective function for orbitals
    parameterized in the MO basis, mo_coeff_buffer is stored to ensure that the MO
    coefficients the library works on are not garbage collected
    """

    mo_coeff_buffer: Optional[np.ndarray] = None


@dataclass
class UpdateOrbsMOPyInterface(UpdateOrbsPyInterface):
    """
    this class provides the Python interface to the orbital updating function for
    orbitals parameterized in the MO basis, which additionally copies the rotated MO
    coefficients back into the array of the caller whenever the library could not
    rotate that array in place
    """

    mo_coeff: Optional[np.ndarray] = None
    mo_coeff_buffer: Optional[np.ndarray] = None

    def __call__(
        self, kappa: np.ndarray, grad: np.ndarray, h_diag: np.ndarray
    ) -> Tuple[float, Callable[[np.ndarray, np.ndarray], None]]:
        try:
            return super().__call__(kappa, grad, h_diag)
        finally:
            # copy rotated MO coefficients back, also on failure since the buffer may
            # have been rotated by then
            if self.mo_coeff is not None and self.mo_coeff_buffer is not None:
                np.copyto(self.mo_coeff, self.mo_coeff_buffer.swapaxes(-1, -2))


def mo_factory(
    mo_coeff: np.ndarray,
    ao_overlap: np.ndarray,
    n_occ: Union[int, Sequence[int]],
    n_particle: int,
    n_ao: int,
    n_mo: int,
    evaluate_dm: EvaluateDMType,
    solver_settings: SolverSettings,
    settings: MOSettings,
    orbsym: Optional[Union[Sequence[int], np.ndarray]] = None,
) -> Tuple[
    Callable[[np.ndarray], float],
    Callable[
        [np.ndarray, np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ],
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

    # define interfaces for callback functions
    evaluate_dm_interface = evaluate_dm_interface_type(
        EvaluateDMInterface(evaluate_dm, n_ao, n_particle, closed_shell, exception)
    )

    # set interfaces for optional callback functions, these need to be set here since
    # the interface might need parameters that are not known when the attribute to
    # settings is set (e.g. n_param)
    settings.set_optional_callback(
        "logger", settings.logger, LoggerInterface, logger_interface_type
    )

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
        POINTER(obj_func_interface_type),
        POINTER(update_orbs_interface_type),
        POINTER(SolverSettingsC),
        POINTER(MOSettingsC),
        POINTER(c_int),
    ]

    # call Fortran function
    obj_func_mo_funptr = obj_func_interface_type()
    update_orbs_mo_funptr = update_orbs_interface_type()
    error = lib.mo_factory(
        mo_coeff_ptr,
        ao_overlap_ptr,
        n_occ_c,
        n_particle,
        n_ao,
        n_mo,
        evaluate_dm_interface,
        byref(obj_func_mo_funptr),
        byref(update_orbs_mo_funptr),
        byref(solver_settings.settings_c),
        byref(settings.settings_c),
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

    # attach the routines the factory has wired into the solver settings
    attach_wired_callbacks(solver_settings, exception)

    return (
        ObjFuncMOPyInterface(
            obj_func_funptr=obj_func_mo_funptr,
            evaluate_dm_interface=evaluate_dm_interface,
            _otr_exception=exception,
            mo_coeff_buffer=mo_coeff_buffer,
        ),
        UpdateOrbsMOPyInterface(
            update_orbs_funptr=update_orbs_mo_funptr,
            _otr_exception=exception,
            saved_objects={
                "evaluate_dm_interface": evaluate_dm_interface,
                "mo_coeff_buffer": mo_coeff_buffer,
            },
            mo_coeff=None if in_place else mo_coeff,
            mo_coeff_buffer=None if in_place else mo_coeff_buffer,
        ),
    )


def mo_deconstructor():
    if not hasattr(lib, "mo_deconstructor"):
        raise RuntimeError(
            "Please reinstall the package with: "
            "CMAKE_FLAGS='-DENABLE_MO=ON' pip install ."
        )

    # define result and argument types
    lib.mo_deconstructor.restype = None
    lib.mo_deconstructor.argtypes = []

    # call Fortran function
    lib.mo_deconstructor()

    return
