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
    logger_interface_type,
    UpdateOrbsInterface,
    LoggerInterface,
    adopt_collector,
    Settings,
    auto_bind_fields,
)
from pyopentrustregion.extensions.common.python_interface import (
    UpdateOrbsPyInterface,
    check_callback_arguments,
)

if TYPE_CHECKING:
    from typing import Tuple, Callable, Optional, Any, Dict

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

    logger: Optional[Callable[[str], None]]
    logger_interface: Any


# ensure that appropriate fields are automatically set in settings_c object
auto_bind_fields(QNSettings)


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
    update_orbs: Callable[
        [np.ndarray, np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ],
    transport: TransportType,
    init_hess: InitHessType,
    n_param: int,
    settings: QNSettings,
) -> Callable[
    [np.ndarray, np.ndarray, np.ndarray],
    Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
]:
    # collector for exceptions raised inside the wrapped user callback; adopted from
    # update_orbs when it is itself factory-produced so a whole chain of factories
    # shares one, and handed to solver through the returned object
    exception = adopt_collector(update_orbs)

    # define interfaces for callback functions, whose arguments the C interface cannot
    # check
    check_callback_arguments(transport, "transport", ("geodesic", "tangent_vector"))
    check_callback_arguments(init_hess, "init_hess", ("vector",))
    update_orbs_interface = update_orbs_interface_type(
        UpdateOrbsInterface(update_orbs, n_param, exception)
    )
    transport_interface = transport_interface_type(
        TransportInterface(transport, n_param, exception)
    )
    init_hess_interface = init_hess_interface_type(
        InitHessInterface(init_hess, n_param, exception)
    )

    # set interfaces for optional callback functions, these need to be set here since
    # the interface might need parameters that are not known when the attribute to
    # settings is set (e.g. n_param)
    settings.set_optional_callback(
        "logger",
        settings.logger,
        LoggerInterface,
        logger_interface_type,
        exception,
    )

    if not hasattr(lib, "qn_factory"):
        raise RuntimeError(
            "Please reinstall the package with: "
            "CMAKE_FLAGS='-DENABLE_QUASI_NEWTON=ON' pip install ."
        )

    # define result and argument types
    lib.qn_factory.restype = c_int
    lib.qn_factory.argtypes = [
        update_orbs_interface_type,
        transport_interface_type,
        init_hess_interface_type,
        c_int,
        POINTER(QNSettingsC),
        POINTER(update_orbs_interface_type),
    ]

    # call Fortran function
    update_orbs_qn_funptr = update_orbs_interface_type()
    error = lib.qn_factory(
        update_orbs_interface,
        transport_interface,
        init_hess_interface,
        n_param,
        byref(settings.settings_c),
        update_orbs_qn_funptr,
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

    return UpdateOrbsPyInterface(
        update_orbs_funptr=update_orbs_qn_funptr,
        _otr_exception=exception,
        saved_objects={
            "update_orbs_interface": update_orbs_interface,
            "transport_interface": transport_interface,
            "init_hess_interface": init_hess_interface,
        },
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
