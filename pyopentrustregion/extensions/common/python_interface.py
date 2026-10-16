# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import annotations

import numpy as np
import weakref
from ctypes import POINTER, byref, c_void_p
from dataclasses import dataclass
from inspect import signature, unwrap
from typing import TYPE_CHECKING
from pyopentrustregion.python_interface import (
    lib,
    c_int,
    c_real,
    callback_type,
    CallbackInterface,
    hess_x_interface_type,
    HessXPyInterface,
    c_buffer,
    write_back,
    raise_callback_error,
)

if TYPE_CHECKING:
    from typing import Tuple, Callable, Optional, Any, Dict, Sequence
    from pyopentrustregion.python_interface import SolverSettings

    GetResponseType = Callable[[np.ndarray, np.ndarray], None]
    EvaluateDMType = Callable[
        [np.ndarray, Optional[np.ndarray], bool],
        Tuple[float, Optional[GetResponseType]],
    ]


# callback function ctypes specifications, ctypes can only deal with simple return
# types so we interface to Fortran subroutines by creating pointers to the relevant data
get_response_interface_type = callback_type(c_int, POINTER(c_real), POINTER(c_real))
evaluate_dm_interface_type = callback_type(
    c_int,
    POINTER(c_real),
    POINTER(c_real),
    POINTER(c_real),
    POINTER(get_response_interface_type),
)


# define interface factories
@dataclass
class GetResponseInterface(CallbackInterface):
    """
    this class provides the interface to write the response to the memory provided for
    a given density matrix
    """

    get_response: GetResponseType
    n_ao: int
    n_particle: int
    closed_shell: bool
    exception: Dict[str, Exception]

    def call(self, dm_ao_ptr, response_ptr) -> int:
        # convert matrix pointers to numpy arrays
        if self.closed_shell:
            dm_ao = np.ctypeslib.as_array(dm_ao_ptr, shape=2 * (self.n_ao,))
            response = np.ctypeslib.as_array(response_ptr, shape=2 * (self.n_ao,))
        else:
            dm_ao = np.ctypeslib.as_array(
                dm_ao_ptr, shape=(self.n_particle, self.n_ao, self.n_ao)
            )
            response = np.ctypeslib.as_array(
                response_ptr, shape=(self.n_particle, self.n_ao, self.n_ao)
            )

        # get response
        try:
            self.get_response(dm_ao, response)
        except Exception as e:
            self.exception["exc"] = e
            return 1

        return 0


@dataclass
class EvaluateDMInterface(CallbackInterface):
    """
    this class provides the interface to the density matrix evaluating function
    """

    evaluate_dm: EvaluateDMType
    n_ao: int
    n_particle: int
    closed_shell: bool
    exception: Dict[str, Exception]
    get_response_funptr: Optional[Any] = None

    def call(self, dm_ao_ptr, energy_ptr, fock_ptr, get_response_funptr) -> int:
        # convert matrix pointers to numpy arrays
        shape = (
            2 * (self.n_ao,)
            if self.closed_shell
            else (self.n_particle, self.n_ao, self.n_ao)
        )
        dm_ao = np.ctypeslib.as_array(dm_ao_ptr, shape=shape)
        fock = np.ctypeslib.as_array(fock_ptr, shape=shape) if fock_ptr else None

        # get energy, and the Fock matrix and response function where wanted
        try:
            energy_ptr[0], get_response = self.evaluate_dm(
                dm_ao, fock, bool(get_response_funptr)
            )
        except Exception as e:
            self.exception["exc"] = e
            return 1

        if get_response_funptr:
            if get_response is None:
                self.exception["exc"] = RuntimeError(
                    "evaluate_dm returned no response function although one was "
                    "requested."
                )
                return 1

            # attach the response interface to the object so that it persists in Python
            # to ensure that it is not garbage collected when the factory completes
            self.get_response_funptr = get_response_interface_type(
                GetResponseInterface(
                    get_response,
                    self.n_ao,
                    self.n_particle,
                    self.closed_shell,
                    self.exception,
                )
            )
            get_response_funptr[0] = self.get_response_funptr

        return 0


@dataclass
class ObjFuncPyInterface:
    """
    this class provides the Python interface to the objective function a factory
    returns, which is called with the context the factory has written into the solver
    settings and keeps these alive, which keep the owner of the factory's object alive
    """

    obj_func_funptr: Any
    saved_objects: Dict[str, Any]
    _otr_exception: Optional[Dict[str, Exception]] = None
    context: Optional[int] = None

    def __call__(self, kappa: np.ndarray) -> float:
        # initialize real
        func = c_real()

        # get pointers to arrays
        kappa_ptr = c_buffer(kappa).ctypes.data_as(POINTER(c_real))

        # objective function
        error = self.obj_func_funptr(kappa_ptr, byref(func), self.context)
        raise_callback_error(
            error, self._otr_exception, "Objective function raised error."
        )

        return func.value


def extension_functions(
    update_orbs_funptr: Any,
    obj_func_funptr: Any,
    solver_settings: SolverSettings,
    exception: Dict[str, Exception],
    copy_back: Optional[Callable[[], None]] = None,
) -> Tuple[UpdateOrbsPyInterface, ObjFuncPyInterface]:
    """
    this function returns the orbital updating function and the objective function a
    factory returns, the Python interfaces to the C functions the C factory returned,
    which are called with the context the factory has written into the solver settings
    and keep these alive, and the orbital update of which copies the orbitals back if
    the library moves a copy of them
    """
    saved_objects = {"solver_settings": solver_settings}
    context = solver_settings.settings_c.context
    return (
        UpdateOrbsPyInterface(
            update_orbs_funptr=update_orbs_funptr,
            saved_objects=saved_objects,
            _otr_exception=exception,
            context=context,
            copy_back=copy_back,
        ),
        ObjFuncPyInterface(
            obj_func_funptr=obj_func_funptr,
            saved_objects=saved_objects,
            _otr_exception=exception,
            context=context,
        ),
    )


def release_extension(extension: Optional[int]) -> None:
    """
    this function frees the extension object a C factory returned, looking up the
    release of the library when it runs, so that tests can replace it
    """
    lib.release_extension.restype = None
    lib.release_extension.argtypes = [c_void_p]
    lib.release_extension(extension)


class ExtensionOwner:
    """
    this class owns the extension object a factory created for a Python host and the
    solver settings it wired, which keep their latest owner alive in place of the owner
    of the object they were wired with before, which the new owner keeps alive, so that
    the functions the factory returned only have to keep the solver settings alive; it
    releases the object once nothing keeps it alive any longer and keeps alive what the
    object calls: the callback interfaces of the factory call, the settings of the
    extension and the record of the solver state
    """

    def __init__(
        self,
        extension: Optional[int],
        solver_settings: SolverSettings,
        saved_objects: Dict[str, Any],
    ):
        self.extension = extension
        self.saved_objects = {
            **saved_objects,
            "wired_before": getattr(solver_settings, "_extension_owner", None),
        }
        self.release = weakref.finalize(self, release_extension, extension)
        solver_settings._extension_owner = self


class OrbitalOwner(ExtensionOwner):
    """
    this class owns the orbital basis object an MO or OAO factory created for a Python
    host like ExtensionOwner and also holds the dimensions of its orbitals, the
    collector of the exceptions of the factory call and, if the library moves a copy of
    the orbitals of the caller, the function that copies them back, which the ARH
    factory, which builds on the object, needs for its own callback interfaces and
    orbital updating function and shares with them
    """

    def __init__(
        self,
        extension: Optional[int],
        solver_settings: SolverSettings,
        saved_objects: Dict[str, Any],
        n_ao: int,
        n_particle: int,
        exception: Dict[str, Exception],
        copy_back: Optional[Callable[[], None]] = None,
    ):
        super().__init__(extension, solver_settings, saved_objects)
        self.n_ao = n_ao
        self.n_particle = n_particle
        self.exception = exception
        self.copy_back = copy_back


def check_callback_arguments(
    callback: Any, name: str, arguments: Sequence[str]
) -> None:
    """
    this function raises if a callback cannot be called with the given positional
    arguments, which the C interface cannot check since it only receives a function
    pointer; callbacks without an inspectable signature are not checked
    """
    try:
        sig = signature(callback)
    except (ValueError, TypeError):
        return
    try:
        sig.bind(*arguments)
    except TypeError:
        raise TypeError(f"{name} has to take ({', '.join(arguments)}).") from None


def check_wrapped_functions(
    update_orbs: Any, obj_func: Any, solver_settings: SolverSettings
) -> None:
    """
    this function raises if the orbital updating function or the objective function a
    factory wraps is not the one the factory that wired the solver settings last
    returned, which reaches its object through the context of the settings, or a
    wrapper of it marked with functools.wraps, or, for settings no factory wired, the
    host's own, which the C interface cannot check, since it cannot tell the level of an
    extension object from a host's context
    """
    context = solver_settings.settings_c.context
    for function, returned_type in (
        (unwrap(update_orbs), UpdateOrbsPyInterface),
        (unwrap(obj_func), ObjFuncPyInterface),
    ):
        returned = isinstance(function, returned_type)
        if (returned and function.context != context) or (
            context is not None and not returned
        ):
            raise ValueError(
                "update_orbs and obj_func have to be the functions the factory that "
                "wired the solver settings last returned, or wrappers of these made "
                "with functools.wraps, or the host's own for settings no factory "
                "wired."
            )


@dataclass
class UpdateOrbsPyInterface:
    """
    this class provides the Python interface to the orbital updating function a factory
    returns, which is called with the context the factory has written into the solver
    settings and keeps these alive, which keep the owner of the factory's object alive
    """

    update_orbs_funptr: Any
    saved_objects: Dict[str, Any]
    hess_x: Optional[Callable[[np.ndarray, np.ndarray], None]] = None
    # collector of the exceptions of the factory call, which the solver adopts, so that
    # it chains the exceptions raised inside the callback functions the object calls
    _otr_exception: Optional[Dict[str, Exception]] = None
    context: Optional[int] = None
    # function copying the orbitals back to the array of the caller after every update
    # if the library moves a copy of them
    copy_back: Optional[Callable[[], None]] = None

    def __call__(
        self, kappa: np.ndarray, grad: np.ndarray, h_diag: np.ndarray
    ) -> Tuple[float, Callable[[np.ndarray, np.ndarray], None]]:
        # initialize real
        func = c_real()

        # get pointers to arrays
        grad_buffer, h_diag_buffer = c_buffer(grad), c_buffer(h_diag)
        kappa_ptr = c_buffer(kappa).ctypes.data_as(POINTER(c_real))
        grad_ptr = grad_buffer.ctypes.data_as(POINTER(c_real))
        h_diag_ptr = h_diag_buffer.ctypes.data_as(POINTER(c_real))

        # initialize Hessian linear transformation function pointer
        hess_x_funptr = hess_x_interface_type()

        # update orbital function and copy the orbitals back, also on failure, since
        # they may have been moved by then
        try:
            error = self.update_orbs_funptr(
                kappa_ptr,
                byref(func),
                grad_ptr,
                h_diag_ptr,
                byref(hess_x_funptr),
                self.context,
            )
        finally:
            if self.copy_back is not None:
                self.copy_back()
        write_back(grad, grad_buffer)
        write_back(h_diag, h_diag_buffer)
        raise_callback_error(
            error, self._otr_exception, "Orbital updating function raised error."
        )

        # keep the interface of the Hessian linear transformation available as an
        # attribute and let it keep this object alive
        self.hess_x = HessXPyInterface(
            hess_x_funptr, self._otr_exception, self.context, _owner=self
        )

        return func.value, self.hess_x
