# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import annotations

import numpy as np
from ctypes import CFUNCTYPE, POINTER, c_bool, c_void_p, c_char, Structure
from typing import TYPE_CHECKING
from pyopentrustregion.python_interface import (
    lib,
    c_int,
    c_real,
    kw_len,
    update_orbs_interface_type,
    hess_x_interface_type,
    update_orbs_interface_factory,
    logger_interface_factory,
    Settings,
    auto_bind_fields,
)

if TYPE_CHECKING:
    from typing import Tuple, Callable, Optional, Any


# callback function ctypes specifications, ctypes can only deal with simple return
# types so we interface to Fortran subroutines by creating pointers to the relevant
# data
transport_interface_type = CFUNCTYPE(c_int, POINTER(c_real), POINTER(c_real))


# define classes corresponding to C structs for settings
class QNSettingsC(Structure):
    _fields_ = [
        ("logger", c_void_p),
        ("initialized", c_bool),
        ("verbose", c_int),
        ("max_points", c_int),
        ("hess_update_scheme", c_char * (kw_len + 1)),
    ]


class QNSettings(Settings):

    c_struct = QNSettingsC
    try:
        init_c_struct = lib.init_qn_settings
    except AttributeError:
        raise AttributeError(
            "Please reinstall the package with: "
            "CMAKE_FLAGS='-DENABLE_HESSIAN_UPDATES=ON' pip install ."
        )

    logger: Optional[Callable[[str], None]]
    logger_interface: Any


# ensure that appropriate fields are automatically set in settings_c object
auto_bind_fields(QNSettings)


def update_orbs_qn_factory(
    update_orbs: Callable[
        [np.ndarray, np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ],
    transport: Callable[[np.ndarray, np.ndarray], None],
    n_param: int,
    settings: QNSettings,
) -> Callable[
    [np.ndarray, np.ndarray, np.ndarray],
    Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
]:
    # define interfaces for callback functions
    update_orbs_interface = update_orbs_interface_factory(update_orbs, n_param)

    @transport_interface_type
    def transport_interface(geodesic_ptr, tangent_vector_ptr):
        """
        this function provides the interface to transport a tangent vector along a
        geodesic and to write it to the memory provided through pointers
        """
        # convert pointers to numpy arrays and deal with row- vs column-major order
        geodesic = np.ctypeslib.as_array(geodesic_ptr, shape=(n_param,))
        tangent_vector = np.ctypeslib.as_array(tangent_vector_ptr, shape=(n_param,))

        # transport and retrieve tangent vector
        try:
            transport(geodesic, tangent_vector)
        except RuntimeError:
            return 1

        return 0

    # set interfaces for optional callback functions, these need to be set here since
    # the interface might need parameters that are not known when the attribute to
    # settings is set (e.g. n_param)
    settings.set_optional_callback("logger", logger_interface_factory(settings.logger))

    if not hasattr(lib, "update_orbs_qn_factory"):
        raise RuntimeError(
            "Please reinstall the package with: "
            "CMAKE_FLAGS='-DENABLE_HESSIAN_UPDATES=ON' pip install ."
        )

    # define result and argument types
    lib.update_orbs_qn_factory.restype = c_int
    lib.update_orbs_qn_factory.argtypes = [
        update_orbs_interface_type,
        transport_interface_type,
        c_int,
        QNSettingsC,
        POINTER(update_orbs_interface_type),
    ]

    # call Fortran function
    update_orbs_qn_funptr = update_orbs_interface_type()
    error = lib.update_orbs_qn_factory(
        update_orbs_interface,
        transport_interface,
        n_param,
        settings.settings_c,
        update_orbs_qn_funptr,
    )

    if error:
        raise RuntimeError(
            "OpenTrustRegion quasi-Newton update factory produced error."
        )

    # python wrapper which converts numpy arrays to ctypes and calls the C function
    def update_orbs_qn_interface(
        kappa: np.ndarray, grad: np.ndarray, h_diag: np.ndarray
    ) -> Tuple[float, Callable[[np.ndarray, np.ndarray], None]]:
        # initialize real
        func = c_real()

        # get pointers to arrays
        kappa_ptr = kappa.ctypes.data_as(POINTER(c_real))
        grad_ptr = grad.ctypes.data_as(POINTER(c_real))
        h_diag_ptr = h_diag.ctypes.data_as(POINTER(c_real))

        # initialize Hessian linear transformation function pointer
        hess_x_qn_funptr = hess_x_interface_type()

        # update orbital function
        error = update_orbs_qn_funptr(
            kappa_ptr, func, grad_ptr, h_diag_ptr, hess_x_qn_funptr
        )
        if error != 0:
            raise RuntimeError("Quasi-Newton orbital updating function raised error.")

        def hess_x_qn_interface(x: np.ndarray, hess_x: np.ndarray):
            # get pointers to arrays
            x_ptr = x.ctypes.data_as(POINTER(c_real))
            hess_x_ptr = hess_x.ctypes.data_as(POINTER(c_real))

            # update orbital function
            error = hess_x_qn_funptr(x_ptr, hess_x_ptr)
            if error != 0:
                raise RuntimeError(
                    "Quasi-Newton Hessian linear transformation function raised error."
                )

        return func.value, hess_x_qn_interface

    # keep all functions alive that are involved in the quasi-Newton orbital update
    update_orbs_qn_interface.update_orbs_interface = update_orbs_interface
    update_orbs_qn_interface.transport_interface = transport_interface
    update_orbs_qn_interface.update_orbs_qn_funptr = update_orbs_qn_funptr

    return update_orbs_qn_interface


def update_orbs_qn_deconstructor():
    if not hasattr(lib, "update_orbs_qn_deconstructor"):
        raise RuntimeError(
            "Please reinstall the package with: "
            "CMAKE_FLAGS='-DENABLE_HESSIAN_UPDATES=ON' pip install ."
        )

    # define result and argument types
    lib.update_orbs_qn_deconstructor.restype = None
    lib.update_orbs_qn_deconstructor.argtypes = []

    # call Fortran function
    lib.update_orbs_qn_deconstructor()

    return
