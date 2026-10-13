# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import annotations

import numpy as np
from ctypes import POINTER, byref, c_bool, c_void_p, c_char, Structure
from dataclasses import dataclass
from typing import TYPE_CHECKING
from pyopentrustregion.python_interface import (
    lib,
    c_int,
    c_real,
    kw_len,
    callback_type,
    CallbackInterface,
    update_orbs_interface_type,
    obj_func_interface_type,
    logger_interface_type,
    UpdateOrbsInterface,
    ObjFuncInterface,
    LoggerInterface,
    LoggerPyInterface,
    adopt_collector,
    Settings,
    SolverSettings,
    SolverSettingsC,
    auto_bind_fields,
    bind_callback_fields,
)
from pyopentrustregion.extensions.common.python_interface import (
    ObjFuncPyInterface,
    UpdateOrbsPyInterface,
    check_callback_arguments,
)

if TYPE_CHECKING:
    from typing import Tuple, Callable, Optional, Dict

    TransportType = Callable[[np.ndarray, np.ndarray], None]
    InitHessType = Callable[[np.ndarray], None]


# callback function ctypes specifications, ctypes can only deal with simple return
# types so we interface to Fortran subroutines by creating pointers to the relevant
# data
transport_interface_type = callback_type(c_int, POINTER(c_real), POINTER(c_real))
init_hess_interface_type = callback_type(c_int, POINTER(c_real))


# define classes corresponding to C structs for settings
class QNSettingsC(Structure):
    _fields_ = [
        ("logger", c_void_p),
        ("initialized", c_bool),
        ("verbose", c_int),
        ("max_points", c_int),
        ("hess_update_scheme", c_char * (kw_len + 1)),
        ("context", c_void_p),
    ]


class QNSettings(Settings):

    c_struct = QNSettingsC
    try:
        init_c_struct = lib.init_qn_settings
    except AttributeError:
        raise AttributeError(
            "Please reinstall the package with: "
            "CMAKE_FLAGS='-DENABLE_QUASI_NEWTON=ON' pip install ."
        )

    callback_fields = {
        "logger": (LoggerInterface, logger_interface_type, LoggerPyInterface)
    }

    logger: Optional[Callable[[str], None]]


# ensure that appropriate fields are automatically set in settings_c object
auto_bind_fields(QNSettings)
bind_callback_fields(QNSettings)


# define interface factories
@dataclass
class TransportInterface(CallbackInterface):
    """
    this class provides the interface to transport a tangent vector along a
    geodesic and to write it to the memory provided through pointers
    """

    transport: TransportType
    n_param: int
    exception: Dict[str, Exception]

    def call(self, geodesic_ptr, tangent_vector_ptr) -> int:
        # convert matrix pointers to numpy arrays
        geodesic = np.ctypeslib.as_array(geodesic_ptr, shape=(self.n_param,))
        tangent_vector = np.ctypeslib.as_array(
            tangent_vector_ptr, shape=(self.n_param,)
        )

        # transport and retrieve tangent vector
        try:
            self.transport(geodesic, tangent_vector)
        except Exception as e:
            self.exception["exc"] = e
            return 1

        return 0


@dataclass
class InitHessInterface(CallbackInterface):
    """
    this class provides the interface to initialize the Hessian and to write it
    to the memory provided through pointers
    """

    init_hess: InitHessType
    n_param: int
    exception: Dict[str, Exception]

    def call(self, vector_ptr) -> int:
        # convert matrix pointers to numpy arrays
        vector = np.ctypeslib.as_array(vector_ptr, shape=(self.n_param,))

        # initialize Hessian and retrieve vector
        try:
            self.init_hess(vector)
        except Exception as e:
            self.exception["exc"] = e
            return 1

        return 0


def qn_factory(
    n_param: int,
    update_orbs: Callable[
        [np.ndarray, np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ],
    obj_func: Callable[[np.ndarray], float],
    transport: TransportType,
    init_hess: InitHessType,
    solver_settings: SolverSettings,
    settings: QNSettings,
) -> Tuple[
    Callable[
        [np.ndarray, np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ],
    Callable[[np.ndarray], float],
]:
    # collector for exceptions raised inside the wrapped user callbacks; adopted from
    # update_orbs or obj_func when these are themselves factory-produced so a whole
    # chain of factories shares one, and handed to solver through the returned objects
    exception = adopt_collector(update_orbs, obj_func)

    # define interfaces for callback functions, whose arguments the C interface cannot
    # check
    check_callback_arguments(transport, "transport", ("geodesic", "tangent_vector"))
    check_callback_arguments(init_hess, "init_hess", ("vector",))
    update_orbs_interface = update_orbs_interface_type(
        UpdateOrbsInterface(update_orbs, n_param, exception)
    )
    obj_func_interface = obj_func_interface_type(
        ObjFuncInterface(obj_func, n_param, exception)
    )
    transport_interface = transport_interface_type(
        TransportInterface(transport, n_param, exception)
    )
    init_hess_interface = init_hess_interface_type(
        InitHessInterface(init_hess, n_param, exception)
    )

    # hand the logging function of the settings the exception collector of the call
    settings.prepare_call(None, exception)

    if not hasattr(lib, "qn_factory"):
        raise RuntimeError(
            "Please reinstall the package with: "
            "CMAKE_FLAGS='-DENABLE_QUASI_NEWTON=ON' pip install ."
        )

    # define result and argument types
    lib.qn_factory.restype = c_int
    lib.qn_factory.argtypes = [
        c_int,
        update_orbs_interface_type,
        obj_func_interface_type,
        transport_interface_type,
        init_hess_interface_type,
        POINTER(update_orbs_interface_type),
        POINTER(obj_func_interface_type),
        POINTER(SolverSettingsC),
        POINTER(QNSettingsC),
    ]

    # call Fortran function
    update_orbs_qn_funptr = update_orbs_interface_type()
    obj_func_qn_funptr = obj_func_interface_type()
    error = lib.qn_factory(
        n_param,
        update_orbs_interface,
        obj_func_interface,
        transport_interface,
        init_hess_interface,
        byref(update_orbs_qn_funptr),
        byref(obj_func_qn_funptr),
        byref(solver_settings.settings_c),
        byref(settings.settings_c),
    )

    if error:
        if "exc" in exception:
            raise RuntimeError(
                f"OpenTrustRegion quasi-Newton factory produced error (code {error})."
            ) from exception["exc"]
        else:
            raise RuntimeError(
                f"OpenTrustRegion quasi-Newton factory produced error (code {error})."
            )

    # the returned functions are called with the context the factory has written into
    # the solver settings, through which they reach the quasi-Newton object, and keep
    # the settings alive, whose callback functions the object calls
    saved_objects = {
        "solver_settings": solver_settings,
        "settings": settings,
        "update_orbs_interface": update_orbs_interface,
        "obj_func_interface": obj_func_interface,
        "transport_interface": transport_interface,
        "init_hess_interface": init_hess_interface,
    }
    return (
        UpdateOrbsPyInterface(
            update_orbs_funptr=update_orbs_qn_funptr,
            saved_objects=saved_objects,
            _otr_exception=exception,
            context=solver_settings.settings_c.context,
        ),
        ObjFuncPyInterface(
            obj_func_funptr=obj_func_qn_funptr,
            saved_objects=saved_objects,
            _otr_exception=exception,
            context=solver_settings.settings_c.context,
        ),
    )


def update_orbs_qn_deconstructor():
    if not hasattr(lib, "update_orbs_qn_deconstructor"):
        raise RuntimeError(
            "Please reinstall the package with: "
            "CMAKE_FLAGS='-DENABLE_QUASI_NEWTON=ON' pip install ."
        )

    # define result and argument types
    lib.update_orbs_qn_deconstructor.restype = None
    lib.update_orbs_qn_deconstructor.argtypes = []

    # call Fortran function
    lib.update_orbs_qn_deconstructor()

    return
