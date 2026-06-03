# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import annotations

import numpy as np
from ctypes import CFUNCTYPE, POINTER, byref, c_bool, c_void_p, Structure
from dataclasses import dataclass
from typing import TYPE_CHECKING
from pyopentrustregion.python_interface import (
    lib,
    c_int,
    c_real,
    update_orbs_interface_type,
    logger_interface_type,
    UpdateOrbsInterface,
    LoggerInterface,
    Settings,
    auto_bind_fields,
)
from pyopentrustregion.extensions.common.python_interface import UpdateOrbsPyInterface

if TYPE_CHECKING:
    from typing import Tuple, Callable, Optional, Any


# callback function ctypes specifications, ctypes can only deal with simple return
# types so we interface to Fortran subroutines by creating pointers to the relevant
# data
change_reference_interface_type = CFUNCTYPE(
    c_int, POINTER(c_real), c_int, POINTER(c_real), POINTER(c_real), POINTER(c_real)
)


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


# define interface factories
@dataclass
class ChangeReferenceInterface:
    """
    this class provides the interface to change the reference and to write the
    parameter, local gradient, and gradient lists to the memory provided through
    pointers
    """

    change_reference: Callable[[np.ndarray, np.ndarray, np.ndarray, np.ndarray], None]
    n_param: int

    def __call__(
        self, new_ref_ptr, n_points, kappa_list_ptr, local_grad_list_ptr, grad_list_ptr
    ) -> int:
        # convert matrix pointers to numpy arrays
        new_ref = np.ctypeslib.as_array(new_ref_ptr, shape=(self.n_param,))
        kappa_list = np.ctypeslib.as_array(
            kappa_list_ptr, shape=(n_points, self.n_param)
        )
        local_grad_list = np.ctypeslib.as_array(
            local_grad_list_ptr, shape=(n_points, self.n_param)
        )
        grad_list = np.ctypeslib.as_array(grad_list_ptr, shape=(n_points, self.n_param))

        # change reference and retrieve parameter, local gradient, and gradient lists
        try:
            self.change_reference(new_ref, kappa_list, local_grad_list, grad_list)
        except RuntimeError:
            return 1

        return 0


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
    update_orbs_interface = update_orbs_interface_type(
        UpdateOrbsInterface(update_orbs, n_param)
    )
    change_reference_interface = change_reference_interface_type(
        ChangeReferenceInterface(change_reference, n_param)
    )

    # set interfaces for optional callback functions, these need to be set here since
    # the interface might need parameters that are not known when the attribute to
    # settings is set (e.g. n_param)
    settings.set_optional_callback(
        "logger", settings.logger, LoggerInterface, logger_interface_type
    )

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

    return UpdateOrbsPyInterface(
        update_orbs_funptr=update_orbs_s_gek_funptr,
        saved_objects={
            "update_orbs_interface": update_orbs_interface,
            "change_reference_interface": change_reference_interface,
        },
    )


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
