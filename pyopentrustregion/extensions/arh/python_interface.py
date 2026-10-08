# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import annotations

import numpy as np
from ctypes import CFUNCTYPE, POINTER, c_bool, c_void_p, c_char, Structure, byref
from dataclasses import dataclass
from inspect import signature
from typing import TYPE_CHECKING
from pyopentrustregion.python_interface import (
    lib,
    c_int,
    c_real,
    kw_len,
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
    ObjFuncPyInterface,
    UpdateOrbsPyInterface,
    attach_wired_callbacks,
    check_callback_arguments,
)
from pyopentrustregion.extensions.mo.python_interface import (
    check_mo_coeff,
    check_orbsym,
    ObjFuncMOPyInterface,
    UpdateOrbsMOPyInterface,
)
from pyopentrustregion.extensions.oao.python_interface import check_dm_ao

if TYPE_CHECKING:
    from typing import Tuple, Callable, Optional, Any, Union, Dict, Sequence

    ARHEvaluateDMType = Callable[
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
arh_evaluate_dm_interface_type = CFUNCTYPE(
    c_int,
    POINTER(c_real),
    POINTER(c_real),
    POINTER(c_real),
    POINTER(c_real),
    POINTER(c_real),
    POINTER(c_real),
)


# define classes corresponding to C structs for settings
class ARHSettingsC(Structure):
    _fields_ = [
        ("logger", c_void_p),
        ("initialized", c_bool),
        ("verbose", c_int),
        ("arh_type", c_char * (kw_len + 1)),
    ]


class ARHSettings(Settings):

    c_struct = ARHSettingsC
    try:
        init_c_struct = lib.init_arh_settings
    except AttributeError:
        raise AttributeError(
            "Please reinstall the package with: "
            "CMAKE_FLAGS='-DENABLE_ARH=ON' pip install ."
        )

    logger: Optional[Callable[[str], None]]
    logger_interface: Any


# ensure that appropriate fields are automatically set in settings_c object
auto_bind_fields(ARHSettings)


# define interface factories
@dataclass
class ARHEvaluateDMInterface:
    """
    this class provides the interface to the density matrix evaluating function with
    Coulomb, exact-exchange and non-linear potential contributions
    """

    evaluate_dm: ARHEvaluateDMType
    n_ao: int
    n_particle: int
    closed_shell: bool
    exception: Dict[str, Exception]

    def __call__(
        self,
        dm_ao_ptr,
        energy_ptr,
        fock_ptr,
        v_coulomb_ptr,
        v_exchange_ptr,
        v_nonlinear_ptr,
    ) -> int:
        # convert matrix pointers to numpy arrays
        shape = (
            2 * (self.n_ao,)
            if self.closed_shell
            else (self.n_particle, self.n_ao, self.n_ao)
        )
        dm_ao = np.ctypeslib.as_array(dm_ao_ptr, shape=shape)
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
        try:
            energy_ptr[0] = self.evaluate_dm(
                dm_ao, fock, v_coulomb, v_exchange, v_nonlinear
            )
        except Exception as e:
            self.exception["exc"] = e
            return 1

        return 0


def arh_factory_mo(
    mo_coeff: np.ndarray,
    ao_overlap: np.ndarray,
    n_occ: Union[int, Sequence[int]],
    n_particle: int,
    n_ao: int,
    n_mo: int,
    evaluate_dm: ARHEvaluateDMType,
    solver_settings: SolverSettings,
    settings: ARHSettings,
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

    # define interfaces for callback functions, whose arguments the C interface cannot
    # check
    check_callback_arguments(
        evaluate_dm,
        "evaluate_dm",
        ("dm", "fock", "v_coulomb", "v_exchange", "v_nonlinear"),
    )
    evaluate_dm_interface = arh_evaluate_dm_interface_type(
        ARHEvaluateDMInterface(evaluate_dm, n_ao, n_particle, closed_shell, exception)
    )

    # set interfaces for optional callback functions, these need to be set here since
    # the interface might need parameters that are not known when the attribute to
    # settings is set (e.g. n_param)
    settings.set_optional_callback(
        "logger", settings.logger, LoggerInterface, logger_interface_type
    )

    if not hasattr(lib, "arh_factory_mo"):
        raise RuntimeError(
            "Please reinstall the package with: "
            "CMAKE_FLAGS='-DENABLE_ARH=ON' pip install ."
        )

    # define result and argument types
    lib.arh_factory_mo.restype = c_int
    lib.arh_factory_mo.argtypes = [
        POINTER(c_real),
        POINTER(c_real),
        POINTER(c_int),
        c_int,
        c_int,
        c_int,
        arh_evaluate_dm_interface_type,
        POINTER(obj_func_interface_type),
        POINTER(update_orbs_interface_type),
        POINTER(SolverSettingsC),
        POINTER(ARHSettingsC),
        POINTER(c_int),
    ]

    # call Fortran function
    obj_func_arh_funptr = obj_func_interface_type()
    update_orbs_arh_funptr = update_orbs_interface_type()
    error = lib.arh_factory_mo(
        mo_coeff_ptr,
        ao_overlap_ptr,
        n_occ_c,
        n_particle,
        n_ao,
        n_mo,
        evaluate_dm_interface,
        byref(obj_func_arh_funptr),
        byref(update_orbs_arh_funptr),
        byref(solver_settings.settings_c),
        byref(settings.settings_c),
        orbsym_c,
    )

    if error:
        if "exc" in exception:
            raise RuntimeError(
                f"OpenTrustRegion ARH MO factory produced error (code {error})."
            ) from exception["exc"]
        else:
            raise RuntimeError(
                f"OpenTrustRegion ARH MO factory produced error (code {error})."
            )

    # attach the routines the factory has wired into the solver settings
    attach_wired_callbacks(solver_settings, exception)

    return (
        ObjFuncMOPyInterface(
            obj_func_funptr=obj_func_arh_funptr,
            evaluate_dm_interface=evaluate_dm_interface,
            _otr_exception=exception,
            mo_coeff_buffer=mo_coeff_buffer,
        ),
        UpdateOrbsMOPyInterface(
            update_orbs_funptr=update_orbs_arh_funptr,
            _otr_exception=exception,
            saved_objects={
                "evaluate_dm_interface": evaluate_dm_interface,
                "mo_coeff_buffer": mo_coeff_buffer,
            },
            mo_coeff=None if in_place else mo_coeff,
            mo_coeff_buffer=None if in_place else mo_coeff_buffer,
        ),
    )


def arh_factory_oao(
    dm_ao: np.ndarray,
    ao_overlap: np.ndarray,
    n_particle: int,
    n_ao: int,
    evaluate_dm: ARHEvaluateDMType,
    solver_settings: SolverSettings,
    settings: ARHSettings,
) -> Tuple[
    Callable[[np.ndarray], float],
    Callable[
        [np.ndarray, np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ],
]:
    # check the density and overlap matrices against the dimensions and determine if
    # closed-shell or open-shell formalism is used
    ao_overlap = check_dm_ao(dm_ao, ao_overlap, n_particle, n_ao)
    closed_shell = n_particle == 1

    # get pointers to arrays
    dm_ao_ptr = dm_ao.ctypes.data_as(POINTER(c_real))
    ao_overlap_ptr = ao_overlap.ctypes.data_as(POINTER(c_real))

    # collector for exceptions raised inside the wrapped user callbacks; adopted from
    # evaluate_dm when it is itself factory-produced so a whole chain of factories
    # shares one, and handed to solver through the returned object
    exception = adopt_collector(evaluate_dm)

    # define interfaces for callback functions, whose arguments the C interface cannot
    # check
    check_callback_arguments(
        evaluate_dm,
        "evaluate_dm",
        ("dm", "fock", "v_coulomb", "v_exchange", "v_nonlinear"),
    )
    evaluate_dm_interface = arh_evaluate_dm_interface_type(
        ARHEvaluateDMInterface(evaluate_dm, n_ao, n_particle, closed_shell, exception)
    )

    # set interfaces for optional callback functions, these need to be set here since
    # the interface might need parameters that are not known when the attribute to
    # settings is set (e.g. n_param)
    settings.set_optional_callback(
        "logger", settings.logger, LoggerInterface, logger_interface_type
    )

    if not hasattr(lib, "arh_factory_oao"):
        raise RuntimeError(
            "Please reinstall the package with: "
            "CMAKE_FLAGS='-DENABLE_ARH=ON' pip install ."
        )

    # define result and argument types
    lib.arh_factory_oao.restype = c_int
    lib.arh_factory_oao.argtypes = [
        POINTER(c_real),
        POINTER(c_real),
        c_int,
        c_int,
        arh_evaluate_dm_interface_type,
        POINTER(obj_func_interface_type),
        POINTER(update_orbs_interface_type),
        POINTER(SolverSettingsC),
        POINTER(ARHSettingsC),
    ]

    # call Fortran function
    obj_func_arh_funptr = obj_func_interface_type()
    update_orbs_arh_funptr = update_orbs_interface_type()
    error = lib.arh_factory_oao(
        dm_ao_ptr,
        ao_overlap_ptr,
        n_particle,
        n_ao,
        evaluate_dm_interface,
        byref(obj_func_arh_funptr),
        byref(update_orbs_arh_funptr),
        byref(solver_settings.settings_c),
        byref(settings.settings_c),
    )

    if error:
        if "exc" in exception:
            raise RuntimeError(
                f"OpenTrustRegion ARH OAO factory produced error (code {error})."
            ) from exception["exc"]
        else:
            raise RuntimeError(
                f"OpenTrustRegion ARH OAO factory produced error (code {error})."
            )

    # attach the routines the factory has wired into the solver settings
    attach_wired_callbacks(solver_settings, exception)

    return (
        ObjFuncPyInterface(
            obj_func_funptr=obj_func_arh_funptr,
            evaluate_dm_interface=evaluate_dm_interface,
            _otr_exception=exception,
        ),
        UpdateOrbsPyInterface(
            update_orbs_funptr=update_orbs_arh_funptr,
            _otr_exception=exception,
            saved_objects={
                "evaluate_dm_interface": evaluate_dm_interface,
            },
        ),
    )


# signatures of the factories the ARH factory selects from
arh_factory_signatures = (
    ("arh_factory_mo", signature(arh_factory_mo)),
    ("arh_factory_oao", signature(arh_factory_oao)),
)


def arh_factory(*args: Any, **kwargs: Any) -> Tuple[
    Callable[[np.ndarray], float],
    Callable[
        [np.ndarray, np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ],
]:
    # select the orbital basis from the arguments, like the generic arh_factory of the
    # Fortran interface: the arguments of arh_factory_oao (an AO density matrix) select
    # the OAO basis and those of arh_factory_mo (MO coefficients with their
    # occupations and number of MOs) the MO basis
    matches = []
    for name, factory_signature in arh_factory_signatures:
        try:
            factory_signature.bind(*args, **kwargs)
        except TypeError:
            continue
        matches.append(name)
    if len(matches) != 1:
        raise TypeError(
            "The arguments have to match those of exactly one of arh_factory_mo and "
            "arh_factory_oao."
        )
    factory = arh_factory_mo if matches[0] == "arh_factory_mo" else arh_factory_oao
    return factory(*args, **kwargs)


def arh_deconstructor():
    if not hasattr(lib, "arh_deconstructor"):
        raise RuntimeError(
            "Please reinstall the package with: "
            "CMAKE_FLAGS='-DENABLE_ARH=ON' pip install ."
        )

    # define result and argument types
    lib.arh_deconstructor.restype = None
    lib.arh_deconstructor.argtypes = []

    # call Fortran function
    lib.arh_deconstructor()

    return
