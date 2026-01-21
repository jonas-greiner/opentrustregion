# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import annotations

import numpy as np
from ctypes import POINTER, c_bool, c_void_p, Structure
from typing import TYPE_CHECKING
from pyopentrustregion.python_interface import (
    lib,
    c_int,
    c_real,
    update_orbs_interface_type,
    hess_x_interface_type,
    update_orbs_interface_factory,
    logger_interface_factory,
    Settings,
    auto_bind_fields,
)
from pyopentrustregion.extensions.common.python_interface import (
    change_reference_interface_type,
    change_reference_interface_factory,
)

if TYPE_CHECKING:
    from typing import Tuple, Callable, Optional, Any


# define classes corresponding to C structs for settings
class SGEKSettingsC(Structure):
    _fields_ = [
        ("logger", c_void_p),
        ("initialized", c_bool),
        ("use_subspace", c_bool),
        ("verbose", c_int),
        ("max_points", c_int),
    ]


class SGEKSettings(Settings):

    c_struct = SGEKSettingsC
    try:
        init_c_struct = lib.init_s_gek_settings
    except AttributeError:
        raise AttributeError(
            "Please reinstall the package with: "
            "CMAKE_FLAGS='-DENABLE_HESSIAN_UPDATES=ON' pip install ."
        )

    logger: Optional[Callable[[str], None]]
    logger_interface: Any


# ensure that appropriate fields are automatically set in settings_c object
auto_bind_fields(SGEKSettings)


def update_orbs_s_gek_factory(
    update_orbs: Callable[
        [np.ndarray, np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ],
    change_reference: Callable[[np.ndarray, np.ndarray, np.ndarray, np.ndarray], None],
    n_param: int,
    settings: SGEKSettings,
) -> Callable[
    [np.ndarray, np.ndarray, np.ndarray],
    Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
]:
    # define interfaces for callback functions
    update_orbs_interface = update_orbs_interface_factory(update_orbs, n_param)
    change_reference_interface = change_reference_interface_factory(
        change_reference, n_param
    )

    # set interfaces for optional callback functions, these need to be set here since
    # the interface might need parameters that are not known when the attribute to
    # settings is set (e.g. n_param)
    settings.set_optional_callback("logger", logger_interface_factory(settings.logger))

    if not hasattr(lib, "update_orbs_s_gek_factory"):
        raise RuntimeError(
            "Please reinstall the package with: "
            "CMAKE_FLAGS='-DENABLE_HESSIAN_UPDATES=ON' pip install ."
        )

    # define result and argument types
    lib.update_orbs_s_gek_factory.restype = c_int
    lib.update_orbs_s_gek_factory.argtypes = [
        update_orbs_interface_type,
        change_reference_interface_type,
        c_int,
        SGEKSettingsC,
        POINTER(update_orbs_interface_type),
    ]

    # call Fortran function
    update_orbs_s_gek_funptr = update_orbs_interface_type()
    error = lib.update_orbs_s_gek_factory(
        update_orbs_interface,
        change_reference_interface,
        n_param,
        settings.settings_c,
        update_orbs_s_gek_funptr,
    )

    if error:
        raise RuntimeError("OpenTrustRegion S-GEK update factory produced error.")

    # python wrapper which converts numpy arrays to ctypes and calls the C function
    def update_orbs_s_gek_interface(
        kappa: np.ndarray, grad: np.ndarray, h_diag: np.ndarray
    ) -> Tuple[float, Callable[[np.ndarray, np.ndarray], None]]:
        # initialize real
        func = c_real()

        # get pointers to arrays
        kappa_ptr = kappa.ctypes.data_as(POINTER(c_real))
        grad_ptr = grad.ctypes.data_as(POINTER(c_real))
        h_diag_ptr = h_diag.ctypes.data_as(POINTER(c_real))

        # initialize Hessian linear transformation function pointer
        hess_x_s_gek_funptr = hess_x_interface_type()

        # update orbital function
        error = update_orbs_s_gek_funptr(
            kappa_ptr, func, grad_ptr, h_diag_ptr, hess_x_s_gek_funptr
        )
        if error != 0:
            raise RuntimeError("S-GEK orbital updating function raised error.")

        def hess_x_s_gek_interface(x: np.ndarray, hess_x: np.ndarray):
            # get pointers to arrays
            x_ptr = x.ctypes.data_as(POINTER(c_real))
            hess_x_ptr = hess_x.ctypes.data_as(POINTER(c_real))

            # update orbital function
            error = hess_x_s_gek_funptr(x_ptr, hess_x_ptr)
            if error != 0:
                raise RuntimeError(
                    "S-GEK Hessian linear transformation function raised error."
                )

        return func.value, hess_x_s_gek_interface

    # keep all functions alive that are involved in the S-GEK orbital update
    update_orbs_s_gek_interface.update_orbs_interface = update_orbs_interface
    update_orbs_s_gek_interface.change_reference_interface = change_reference_interface
    update_orbs_s_gek_interface.update_orbs_s_gek_funptr = update_orbs_s_gek_funptr

    return update_orbs_s_gek_interface


def update_orbs_s_gek_deconstructor():
    if not hasattr(lib, "update_orbs_s_gek_deconstructor"):
        raise RuntimeError(
            "Please reinstall the package with: "
            "CMAKE_FLAGS='-DENABLE_HESSIAN_UPDATES=ON' pip install ."
        )

    # define result and argument types
    lib.update_orbs_s_gek_deconstructor.restype = None
    lib.update_orbs_s_gek_deconstructor.argtypes = []

    # call Fortran function
    lib.update_orbs_s_gek_deconstructor()

    return
