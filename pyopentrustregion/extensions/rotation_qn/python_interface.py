# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import annotations

import numpy as np
from ctypes import POINTER, byref, c_bool, c_void_p, c_char, Structure
from typing import TYPE_CHECKING
from pyopentrustregion.python_interface import (
    check_struct_sizes,
    lib,
    c_int,
    c_real,
    kw_len,
    callback_type,
    SettingsCallbackInterface,
    CFunctionInterface,
    c_buffer,
    write_back,
    update_orbs_interface_type,
    obj_func_interface_type,
    logger_interface_type,
    hess_x_interface_type,
    precond_interface_type,
    UpdateOrbsInterface,
    ObjFuncInterface,
    LoggerInterface,
    LoggerPyInterface,
    SettingsHessXInterface,
    HessXPyInterface,
    PrecondInterface,
    PrecondPyInterface,
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
    ExtensionOwner,
    extension_functions,
    check_wrapped_functions,
)

if TYPE_CHECKING:
    from typing import Tuple, Callable, Optional

    TransportType = Callable[[np.ndarray, np.ndarray], None]
    ObjFuncGradType = Callable[[np.ndarray, np.ndarray], float]


# callback function ctypes specifications, ctypes can only deal with simple return
# types so we interface to Fortran subroutines by creating pointers to the relevant
# data
obj_func_grad_interface_type = callback_type(
    c_int, POINTER(c_real), POINTER(c_real), POINTER(c_real)
)
transport_interface_type = callback_type(c_int, POINTER(c_real), POINTER(c_real))


class ObjFuncGradInterface(SettingsCallbackInterface):
    """
    this class provides the interface to the objective function of the rotation QN
    settings, which also writes the gradient at the trial point to the memory provided
    through a pointer
    """

    def call(self, kappa_ptr, grad_ptr, func_ptr) -> int:
        # convert pointers to numpy arrays and make the step read-only
        kappa = np.ctypeslib.as_array(kappa_ptr, shape=(self.n_param,))
        kappa.flags.writeable = False
        grad = np.ctypeslib.as_array(grad_ptr, shape=(self.n_param,))

        func_ptr[0] = self.func(kappa, grad)

        return 0


class TransportInterface(SettingsCallbackInterface):
    """
    this class provides the interface to the transport function of the rotation QN
    settings, which transports a tangent vector along a geodesic in place
    """

    def call(self, geodesic_ptr, tangent_vector_ptr) -> int:
        # convert pointers to numpy arrays and make the geodesic read-only
        geodesic = np.ctypeslib.as_array(geodesic_ptr, shape=(self.n_param,))
        geodesic.flags.writeable = False
        tangent_vector = np.ctypeslib.as_array(
            tangent_vector_ptr, shape=(self.n_param,)
        )

        self.func(geodesic, tangent_vector)

        return 0


class ObjFuncGradPyInterface(CFunctionInterface):
    """
    this class provides the Python interface to an objective function C function that
    also returns the gradient at the trial point
    """

    def __call__(self, kappa: np.ndarray, grad: np.ndarray) -> float:
        func = c_real()
        kappa_buffer, grad_buffer = c_buffer(kappa), c_buffer(grad)
        LibraryCall("objective function", kappa_buffer.size).run(
            self.funptr,
            kappa_buffer.ctypes.data_as(POINTER(c_real)),
            grad_buffer.ctypes.data_as(POINTER(c_real)),
            byref(func),
            self.context,
        )
        write_back(grad, grad_buffer)
        return func.value


class TransportPyInterface(CFunctionInterface):
    """
    this class provides the Python interface to a transport C function
    """

    def __call__(self, geodesic: np.ndarray, tangent_vector: np.ndarray):
        geodesic_buffer = c_buffer(geodesic)
        tangent_vector_buffer = c_buffer(tangent_vector)
        LibraryCall("transport function", geodesic_buffer.size).run(
            self.funptr,
            geodesic_buffer.ctypes.data_as(POINTER(c_real)),
            tangent_vector_buffer.ctypes.data_as(POINTER(c_real)),
            self.context,
        )
        write_back(tangent_vector, tangent_vector_buffer)


# define classes corresponding to C structs for settings
class RotationQNSettingsC(Structure):
    _fields_ = [
        ("init_hess", c_void_p),
        ("init_precond", c_void_p),
        ("obj_func_grad", c_void_p),
        ("transport", c_void_p),
        ("logger", c_void_p),
        ("initialized", c_bool),
        ("verbose", c_int),
        ("max_points", c_int),
        ("qn_method", c_char * (kw_len + 1)),
        ("context", c_void_p),
    ]


class RotationQNSettings(Settings):

    c_struct = RotationQNSettingsC
    try:
        init_c_struct = lib.init_rotation_qn_settings
    except AttributeError:
        raise AttributeError(
            "Please reinstall the package with: "
            "pip install -C cmake.define.ENABLE_ROTATION_QN=ON ."
        )

    callback_fields = {
        "init_hess": (SettingsHessXInterface, hess_x_interface_type, HessXPyInterface),
        "init_precond": (PrecondInterface, precond_interface_type, PrecondPyInterface),
        "obj_func_grad": (
            ObjFuncGradInterface,
            obj_func_grad_interface_type,
            ObjFuncGradPyInterface,
        ),
        "transport": (
            TransportInterface,
            transport_interface_type,
            TransportPyInterface,
        ),
        "logger": (LoggerInterface, logger_interface_type, LoggerPyInterface),
    }

    init_hess: Optional[Callable[[np.ndarray, np.ndarray], None]]
    init_precond: Optional[Callable[[np.ndarray, float, np.ndarray], None]]
    obj_func_grad: Optional[ObjFuncGradType]
    transport: Optional[TransportType]
    logger: Optional[Callable[[str], None]]


# ensure that appropriate fields are automatically set in settings_c object
auto_bind_fields(RotationQNSettings)
bind_callback_fields(RotationQNSettings)

# check that the settings structure has the size the library exports for it
check_struct_sizes(lib, [("rotation_qn_settings", RotationQNSettingsC)])


def rotation_qn_factory(
    n_param: int,
    update_orbs: Callable[
        [np.ndarray, np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ],
    obj_func: Callable[[np.ndarray], float],
    solver_settings: SolverSettings,
    solver_state: SolverState,
    settings: Optional[RotationQNSettings] = None,
) -> Tuple[
    Callable[
        [np.ndarray, np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ],
    Callable[[np.ndarray], float],
]:
    # the C interface cannot check whether the orbital updating function and the
    # objective function belong to the solver settings
    check_wrapped_functions(update_orbs, obj_func, solver_settings)

    # define interfaces for callback functions, whose arguments the C interface cannot
    # check
    update_orbs_interface = update_orbs_interface_type(UpdateOrbsInterface(update_orbs))
    obj_func_interface = obj_func_interface_type(ObjFuncInterface(obj_func))

    # define result and argument types
    lib.rotation_qn_factory.restype = c_int
    lib.rotation_qn_factory.argtypes = [
        c_int,
        update_orbs_interface_type,
        obj_func_interface_type,
        POINTER(update_orbs_interface_type),
        POINTER(obj_func_interface_type),
        POINTER(c_void_p),
        POINTER(SolverSettingsC),
        POINTER(SolverStateC),
        POINTER(RotationQNSettingsC),
    ]

    # call Fortran function
    update_orbs_rotation_qn_funptr = update_orbs_interface_type()
    obj_func_rotation_qn_funptr = obj_func_interface_type()
    extension = c_void_p()
    LibraryCall("rotation QN factory", n_param).run(
        lib.rotation_qn_factory,
        n_param,
        update_orbs_interface,
        obj_func_interface,
        byref(update_orbs_rotation_qn_funptr),
        byref(obj_func_rotation_qn_funptr),
        byref(extension),
        byref(solver_settings.settings_c),
        byref(solver_state.state_c),
        None if settings is None else byref(settings.settings_c),
    )

    # the owner of the rotation QN object, which the solver settings keep alive, keeps
    # alive what the object calls, and the returned functions are called with the
    # context the factory has written into the solver settings, through which they
    # reach the object, and keep the solver settings alive
    ExtensionOwner(
        extension.value,
        solver_settings,
        {
            "settings": settings,
            "solver_state": solver_state,
            "update_orbs_interface": update_orbs_interface,
            "obj_func_interface": obj_func_interface,
        },
    )
    return extension_functions(
        update_orbs_rotation_qn_funptr,
        obj_func_rotation_qn_funptr,
        solver_settings,
    )
