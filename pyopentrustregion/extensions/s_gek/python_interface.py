# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import annotations

import numpy as np
from ctypes import POINTER, byref, c_bool, c_void_p, Structure
from dataclasses import dataclass
from typing import TYPE_CHECKING
from pyopentrustregion.python_interface import (
    lib,
    c_int,
    c_real,
    callback_type,
    CallbackInterface,
    update_orbs_interface_type,
    obj_func_interface_type,
    logger_interface_type,
    UpdateOrbsInterface,
    ObjFuncInterface,
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
    detach_python_callbacks,
)

if TYPE_CHECKING:
    from typing import Tuple, Callable, Optional, Any, Dict

    ChangeReferenceType = Callable[
        [np.ndarray, np.ndarray, np.ndarray, np.ndarray], None
    ]


# callback function ctypes specifications, ctypes can only deal with simple return
# types so we interface to Fortran subroutines by creating pointers to the relevant
# data
change_reference_interface_type = callback_type(
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
        ("context", c_void_p),
    ]


class SGEKSettings(Settings):

    c_struct = SGEKSettingsC
    try:
        init_c_struct = lib.init_s_gek_settings
    except AttributeError:
        raise AttributeError(
            "Please reinstall the package with: "
            "CMAKE_FLAGS='-DENABLE_S_GEK=ON' pip install ."
        )

    logger: Optional[Callable[[str], None]]
    logger_interface: Any


# ensure that appropriate fields are automatically set in settings_c object
auto_bind_fields(SGEKSettings)


# define interface factories
@dataclass
class ChangeReferenceInterface(CallbackInterface):
    """
    this class provides the interface to change the reference and to write the
    parameter, local gradient, and gradient lists to the memory provided through
    pointers
    """

    change_reference: ChangeReferenceType
    n_param: int
    exception: Dict[str, Exception]

    def call(
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
        except Exception as e:
            self.exception["exc"] = e
            return 1

        return 0


def s_gek_factory(
    n_param: int,
    update_orbs: Callable[
        [np.ndarray, np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ],
    obj_func: Callable[[np.ndarray], float],
    change_reference: ChangeReferenceType,
    solver_settings: SolverSettings,
    settings: SGEKSettings,
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
    check_callback_arguments(
        change_reference,
        "change_reference",
        ("new_ref", "kappa_list", "local_grad_list", "grad_list"),
    )
    update_orbs_interface = update_orbs_interface_type(
        UpdateOrbsInterface(update_orbs, n_param, exception)
    )
    obj_func_interface = obj_func_interface_type(
        ObjFuncInterface(obj_func, n_param, exception)
    )
    change_reference_interface = change_reference_interface_type(
        ChangeReferenceInterface(change_reference, n_param, exception)
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

    # the factory keeps calling the callback functions of the solver settings it wires,
    # so those the solver built from the Python functions of the settings are removed
    # first, the solver builds them anew
    detach_python_callbacks(solver_settings)

    if not hasattr(lib, "s_gek_factory"):
        raise RuntimeError(
            "Please reinstall the package with: "
            "CMAKE_FLAGS='-DENABLE_S_GEK=ON' pip install ."
        )

    # define result and argument types
    lib.s_gek_factory.restype = c_int
    lib.s_gek_factory.argtypes = [
        c_int,
        update_orbs_interface_type,
        obj_func_interface_type,
        change_reference_interface_type,
        POINTER(update_orbs_interface_type),
        POINTER(obj_func_interface_type),
        POINTER(SolverSettingsC),
        POINTER(SGEKSettingsC),
    ]

    # call Fortran function
    update_orbs_s_gek_funptr = update_orbs_interface_type()
    obj_func_s_gek_funptr = obj_func_interface_type()
    error = lib.s_gek_factory(
        n_param,
        update_orbs_interface,
        obj_func_interface,
        change_reference_interface,
        byref(update_orbs_s_gek_funptr),
        byref(obj_func_s_gek_funptr),
        byref(solver_settings.settings_c),
        byref(settings.settings_c),
    )

    if error:
        if "exc" in exception:
            raise RuntimeError(
                f"OpenTrustRegion S-GEK factory produced error (code {error})."
            ) from exception["exc"]
        else:
            raise RuntimeError(
                f"OpenTrustRegion S-GEK factory produced error (code {error})."
            )

    # attach the routines the factory has wired into the solver settings
    attach_wired_callbacks(solver_settings, exception)

    # the returned functions are called with the context the factory has written into
    # the solver settings, through which they reach the S-GEK object
    saved_objects = {
        "update_orbs_interface": update_orbs_interface,
        "obj_func_interface": obj_func_interface,
        "change_reference_interface": change_reference_interface,
    }
    return (
        UpdateOrbsPyInterface(
            update_orbs_funptr=update_orbs_s_gek_funptr,
            saved_objects=saved_objects,
            _otr_exception=exception,
            context=solver_settings.settings_c.context,
        ),
        ObjFuncPyInterface(
            obj_func_funptr=obj_func_s_gek_funptr,
            saved_objects=saved_objects,
            _otr_exception=exception,
            context=solver_settings.settings_c.context,
        ),
    )


def update_orbs_s_gek_deconstructor():
    if not hasattr(lib, "update_orbs_s_gek_deconstructor"):
        raise RuntimeError(
            "Please reinstall the package with: "
            "CMAKE_FLAGS='-DENABLE_S_GEK=ON' pip install ."
        )

    # define result and argument types
    lib.update_orbs_s_gek_deconstructor.restype = None
    lib.update_orbs_s_gek_deconstructor.argtypes = []

    # call Fortran function
    lib.update_orbs_s_gek_deconstructor()

    return
