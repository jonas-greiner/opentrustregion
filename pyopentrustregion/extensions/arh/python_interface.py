# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import annotations

import sys
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
    project_interface_type,
    precond_interface_type,
    precond_pd_interface_type,
    logger_interface_type,
    LoggerInterface,
    adopt_collector,
    Settings,
    auto_bind_fields,
)
from pyopentrustregion.extensions.common.python_interface import UpdateOrbsPyInterface
from pyopentrustregion.extensions.oao.python_interface import (
    ObjFuncPyInterface,
    PrecondPyInterface,
    PrecondPDPyInterface,
    ProjectPyInterface,
)

if TYPE_CHECKING:
    from typing import Tuple, Callable, Optional, Any, Union, TypeGuard, Dict

    EvaluateDMCSType = Callable[
        [np.ndarray, Optional[np.ndarray], Optional[np.ndarray]], float
    ]
    EvaluateDMOSType = Callable[
        [
            np.ndarray,
            Optional[np.ndarray],
            Optional[np.ndarray],
            Optional[np.ndarray],
            Optional[np.ndarray],
        ],
        float,
    ]


# type guards
def is_evaluate_dm_cs(func: Any) -> TypeGuard[EvaluateDMCSType]:
    try:
        sig = signature(func)
        return len(sig.parameters) == 3
    except (ValueError, TypeError):
        return False


def is_evaluate_dm_os(func: Any) -> TypeGuard[EvaluateDMOSType]:
    try:
        sig = signature(func)
        return len(sig.parameters) == 5
    except (ValueError, TypeError):
        return False


# callback function ctypes specifications, ctypes can only deal with simple return
# types so we interface to Fortran subroutines by creating pointers to the relevant
# data
evaluate_dm_os_interface_type = CFUNCTYPE(
    c_int,
    POINTER(c_real),
    POINTER(c_real),
    POINTER(c_real),
    POINTER(c_real),
    POINTER(c_real),
    POINTER(c_real),
)
evaluate_dm_cs_interface_type = CFUNCTYPE(
    c_int,
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
class EvaluateDMCSInterface:
    """
    this class provides the interface to the density matrix evaluating function with a
    separate non-linear potential contribution for the closed-shell case
    """

    evaluate_dm_cs: EvaluateDMCSType
    n_ao: int
    n_particle: int
    closed_shell: bool
    exception: Dict[str, Exception]

    def __call__(self, dm_ao_ptr, energy_ptr, fock_ptr, v_nonlinear_ptr) -> int:
        # convert matrix pointers to numpy arrays
        shape = 2 * (self.n_ao,)
        dm_ao = np.ctypeslib.as_array(dm_ao_ptr, shape=shape)
        fock = np.ctypeslib.as_array(fock_ptr, shape=shape) if fock_ptr else None
        v_nonlinear = (
            np.ctypeslib.as_array(v_nonlinear_ptr, shape=shape)
            if v_nonlinear_ptr
            else None
        )

        # get energy, and the Fock matrix, and non-linear potential where wanted
        try:
            energy_ptr[0] = self.evaluate_dm_cs(dm_ao, fock, v_nonlinear)
        except Exception as e:
            self.exception["exc"] = e
            return 1

        return 0


@dataclass
class EvaluateDMOSInterface:
    """
    this class provides the interface density matrix evaluating function with same- and
    opposite-spin and non-linear potential contributions for the open-shell case
    """

    evaluate_dm_os: EvaluateDMOSType
    n_ao: int
    n_particle: int
    closed_shell: bool
    exception: Dict[str, Exception]

    def __call__(
        self,
        dm_ao_ptr,
        energy_ptr,
        fock_ptr,
        v_same_spin_ptr,
        v_opposite_spin_ptr,
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
        v_same_spin = (
            np.ctypeslib.as_array(v_same_spin_ptr, shape=shape)
            if v_same_spin_ptr
            else None
        )
        v_opposite_spin = (
            np.ctypeslib.as_array(v_opposite_spin_ptr, shape=shape)
            if v_opposite_spin_ptr
            else None
        )
        v_nonlinear = (
            np.ctypeslib.as_array(v_nonlinear_ptr, shape=shape)
            if v_nonlinear_ptr
            else None
        )

        # get energy, and the Fock matrix, and same-, opposite-spin, and non-linear
        # potentials where wanted
        try:
            energy_ptr[0] = self.evaluate_dm_os(
                dm_ao, fock, v_same_spin, v_opposite_spin, v_nonlinear
            )
        except Exception as e:
            self.exception["exc"] = e
            return 1

        return 0


def arh_factory(
    dm_ao: np.ndarray,
    ao_overlap: np.ndarray,
    n_particle: int,
    n_ao: int,
    evaluate_dm: Union[EvaluateDMCSType, EvaluateDMOSType],
    settings: ARHSettings,
) -> Tuple[
    Callable[[np.ndarray], float],
    Callable[
        [np.ndarray, np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ],
    Callable[[np.ndarray, float, np.ndarray], None],
    Callable[[np.ndarray, np.ndarray], None],
    Callable[[np.ndarray], None],
]:
    # get pointers to arrays
    dm_ao_ptr = dm_ao.ctypes.data_as(POINTER(c_real))
    ao_overlap_ptr = ao_overlap.ctypes.data_as(POINTER(c_real))

    # determine if closed-shell or open-shell formalism is used
    closed_shell = dm_ao.ndim == 2

    # collector for exceptions raised inside the wrapped user callbacks; adopted from
    # evaluate_dm when it is itself factory-produced so a whole chain of factories
    # shares one, and handed to solver through the returned object
    exception = adopt_collector(evaluate_dm)

    # define interfaces for callback functions
    if is_evaluate_dm_cs(evaluate_dm):
        evaluate_dm_cs_interface = evaluate_dm_cs_interface_type(
            EvaluateDMCSInterface(
                evaluate_dm, n_ao, n_particle, closed_shell, exception
            )
        )
    elif is_evaluate_dm_os(evaluate_dm):
        evaluate_dm_os_interface = evaluate_dm_os_interface_type(
            EvaluateDMOSInterface(
                evaluate_dm, n_ao, n_particle, closed_shell, exception
            )
        )

    # set interfaces for optional callback functions, these need to be set here since
    # the interface might need parameters that are not known when the attribute to
    # settings is set (e.g. n_param)
    settings.set_optional_callback(
        "logger", settings.logger, LoggerInterface, logger_interface_type
    )

    if not hasattr(lib, "arh_factory"):
        raise RuntimeError(
            "Please reinstall the package with: "
            "CMAKE_FLAGS='-DENABLE_ARH=ON' pip install ."
        )

    # define result and argument types
    lib.arh_factory.restype = c_int
    lib.arh_factory.argtypes = [
        POINTER(c_real),
        POINTER(c_real),
        c_int,
        c_int,
        (
            evaluate_dm_cs_interface_type
            if closed_shell
            else evaluate_dm_os_interface_type
        ),
        POINTER(obj_func_interface_type),
        POINTER(update_orbs_interface_type),
        POINTER(precond_interface_type),
        POINTER(precond_pd_interface_type),
        POINTER(project_interface_type),
        POINTER(ARHSettingsC),
    ]

    # call Fortran function
    obj_func_arh_funptr = obj_func_interface_type()
    update_orbs_arh_funptr = update_orbs_interface_type()
    precond_arh_funptr = precond_interface_type()
    precond_pd_arh_funptr = precond_pd_interface_type()
    project_arh_funptr = project_interface_type()
    error = lib.arh_factory(
        dm_ao_ptr,
        ao_overlap_ptr,
        n_particle,
        n_ao,
        evaluate_dm_cs_interface if closed_shell else evaluate_dm_os_interface,
        byref(obj_func_arh_funptr),
        byref(update_orbs_arh_funptr),
        byref(precond_arh_funptr),
        byref(precond_pd_arh_funptr),
        byref(project_arh_funptr),
        byref(settings.settings_c),
    )

    if error:
        if "exc" in exception:
            raise RuntimeError(
                f"OpenTrustRegion ARH factory produced error (code {error})."
            ) from exception["exc"]
        else:
            raise RuntimeError(
                f"OpenTrustRegion ARH factory produced error (code {error})."
            )

    return (
        ObjFuncPyInterface(
            obj_func_funptr=obj_func_arh_funptr,
            evaluate_dm_interface=(
                evaluate_dm_cs_interface if closed_shell else evaluate_dm_os_interface
            ),
            _otr_exception=exception,
        ),
        UpdateOrbsPyInterface(
            update_orbs_funptr=update_orbs_arh_funptr,
            _otr_exception=exception,
            saved_objects={
                "evaluate_dm_interface": (
                    evaluate_dm_cs_interface
                    if closed_shell
                    else evaluate_dm_os_interface
                ),
            },
        ),
        PrecondPyInterface(precond_funptr=precond_arh_funptr, _otr_exception=exception),
        PrecondPDPyInterface(
            precond_pd_funptr=precond_pd_arh_funptr, _otr_exception=exception
        ),
        ProjectPyInterface(project_funptr=project_arh_funptr, _otr_exception=exception),
    )


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
