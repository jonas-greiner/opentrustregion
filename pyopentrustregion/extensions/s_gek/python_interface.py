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
    check_struct_sizes,
    lib,
    c_int,
    c_real,
    kw_len,
    callback_type,
    CallbackInterface,
    update_orbs_interface_type,
    obj_func_interface_type,
    logger_interface_type,
    hess_x_interface_type,
    UpdateOrbsInterface,
    ObjFuncInterface,
    LoggerInterface,
    LoggerPyInterface,
    SettingsHessXInterface,
    HessXPyInterface,
    LibraryCall,
    Settings,
    SolverSettings,
    SolverSettingsC,
    SolverState,
    SolverStateC,
    auto_bind_fields,
    bind_callback_fields,
)
from pyopentrustregion.extensions.common.python_interface import (
    obj_func_grad_interface_type,
    ObjFuncGradInterface,
    ObjFuncGradPyInterface,
)

if TYPE_CHECKING:
    from typing import Tuple, Callable, Optional
    from pyopentrustregion.extensions.common.python_interface import ObjFuncGradType

    ChangeReferenceType = Callable[
        [np.ndarray, np.ndarray, np.ndarray, np.ndarray], None
    ]


# callback function ctypes specifications, ctypes can only deal with simple return
# types so we interface to Fortran subroutines by creating pointers to the relevant
# data, while the common interface of the extensions holds that of the objective
# function that also returns the gradient at the trial point
change_reference_interface_type = callback_type(
    c_int, POINTER(c_real), c_int, POINTER(c_real), POINTER(c_real), POINTER(c_real)
)


@dataclass
class ChangeReferenceInterface(CallbackInterface):
    """
    this class provides the interface for the chart transformation, which re-expresses
    the positions, local gradients and gradients of the points, one point per row, in
    the chart of a new reference in place
    """

    change_reference: ChangeReferenceType

    def call(
        self, new_ref_ptr, n_points, kappa_list_ptr, local_grad_list_ptr, grad_list_ptr
    ) -> int:
        # convert pointers to numpy arrays and make the new reference read-only
        new_ref = np.ctypeslib.as_array(new_ref_ptr, shape=(self.n_param,))
        new_ref.flags.writeable = False
        kappa_list = np.ctypeslib.as_array(
            kappa_list_ptr, shape=(n_points, self.n_param)
        )
        local_grad_list = np.ctypeslib.as_array(
            local_grad_list_ptr, shape=(n_points, self.n_param)
        )
        grad_list = np.ctypeslib.as_array(grad_list_ptr, shape=(n_points, self.n_param))

        self.change_reference(new_ref, kappa_list, local_grad_list, grad_list)

        return 0


# define classes corresponding to C structs for settings
class SGEKSettingsC(Structure):
    _fields_ = [
        ("obj_func_grad", c_void_p),
        ("init_hess", c_void_p),
        ("logger", c_void_p),
        ("initialized", c_bool),
        ("excess_dispersion", c_bool),
        ("verbose", c_int),
        ("inner_solver", c_char * (kw_len + 1)),
        ("prior_update", c_char * (kw_len + 1)),
        ("acceptance", c_char * (kw_len + 1)),
        ("context", c_void_p),
    ]


class SGEKSettings(Settings):

    c_struct = SGEKSettingsC
    try:
        init_c_struct = lib.init_s_gek_settings
    except AttributeError:
        raise AttributeError(
            "Please reinstall the package with: "
            "pip install -C cmake.define.ENABLE_S_GEK=ON ."
        )

    callback_fields = {
        "obj_func_grad": (
            ObjFuncGradInterface,
            obj_func_grad_interface_type,
            ObjFuncGradPyInterface,
        ),
        "init_hess": (SettingsHessXInterface, hess_x_interface_type, HessXPyInterface),
        "logger": (LoggerInterface, logger_interface_type, LoggerPyInterface),
    }

    obj_func_grad: Optional[ObjFuncGradType]
    init_hess: Optional[Callable[[np.ndarray, np.ndarray], None]]
    logger: Optional[Callable[[str], None]]


# ensure that appropriate fields are automatically set in settings_c object
auto_bind_fields(SGEKSettings)
bind_callback_fields(SGEKSettings)

# check that the settings structure has the size the library exports for it
check_struct_sizes(lib, [("s_gek_settings", SGEKSettingsC)])


def s_gek_solver(
    update_orbs: Callable[
        [np.ndarray, np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ],
    obj_func: Callable[[np.ndarray], float],
    change_reference: ChangeReferenceType,
    n_param: int,
    settings: Optional[SolverSettings] = None,
    state: Optional[SolverState] = None,
    s_gek_settings: Optional[SGEKSettings] = None,
):
    # call into the library, which keeps what the callback functions raise
    library_call = LibraryCall("S-GEK solver", n_param)

    # define interfaces for callback functions
    obj_func_interface = obj_func_interface_type(ObjFuncInterface(obj_func))
    update_orbs_interface = update_orbs_interface_type(UpdateOrbsInterface(update_orbs))
    change_reference_interface = change_reference_interface_type(
        ChangeReferenceInterface(change_reference)
    )

    # define result and argument types
    lib.s_gek_solver.restype = c_int
    lib.s_gek_solver.argtypes = [
        update_orbs_interface_type,
        obj_func_interface_type,
        change_reference_interface_type,
        c_int,
        POINTER(SolverSettingsC),
        POINTER(SolverStateC),
        POINTER(SGEKSettingsC),
    ]

    # call Fortran function
    library_call.run(
        lib.s_gek_solver,
        update_orbs_interface,
        obj_func_interface,
        change_reference_interface,
        n_param,
        byref(settings.settings_c) if settings is not None else None,
        byref(state.state_c) if state is not None else None,
        byref(s_gek_settings.settings_c) if s_gek_settings is not None else None,
    )
