# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import annotations

import numpy as np
from ctypes import POINTER, byref
from dataclasses import dataclass
from inspect import signature
from typing import TYPE_CHECKING
from pyopentrustregion.python_interface import (
    c_int,
    c_real,
    callback_type,
    CallbackInterface,
    hess_x_interface_type,
    HessXPyInterface,
)

if TYPE_CHECKING:
    from typing import Tuple, Callable, Optional, Any, Dict, Sequence

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
    this class provides the Python interface to the objective function, which is
    called with the context the factory has written into the solver settings, the
    callback interfaces and the settings whose callback functions the factory's object
    calls are stored to ensure that they are not garbage collected when the factory
    completes
    """

    obj_func_funptr: Any
    saved_objects: Dict[str, Any]
    _otr_exception: Optional[Dict[str, Exception]] = None
    context: Optional[int] = None

    def __call__(self, kappa: np.ndarray) -> float:
        # initialize real
        func = c_real()

        # get pointers to arrays
        kappa_ptr = kappa.ctypes.data_as(POINTER(c_real))

        # objective function
        error = self.obj_func_funptr(kappa_ptr, byref(func), self.context)
        if error != 0:
            if self._otr_exception is not None and "exc" in self._otr_exception:
                raise RuntimeError(
                    "Objective function raised error."
                ) from self._otr_exception["exc"]
            raise RuntimeError("Objective function raised error.")

        return func.value


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


@dataclass
class UpdateOrbsPyInterface:
    """
    this class provides the Python interface to the orbital updating function, which is
    called with the context the factory has written into the solver settings, the
    callback interfaces and the settings whose callback functions the factory's object
    calls are stored to ensure that they are not garbage collected when the factory
    completes
    """

    update_orbs_funptr: Any
    saved_objects: Dict[str, Any]
    hess_x: Optional[Callable[[np.ndarray, np.ndarray], None]] = None
    # collector shared with the wrapped user callback, so that solver can adopt it and
    # chain exceptions raised inside that callback; None when this object does not wrap
    # a user callback of its own
    _otr_exception: Optional[Dict[str, Exception]] = None
    context: Optional[int] = None

    def __call__(
        self, kappa: np.ndarray, grad: np.ndarray, h_diag: np.ndarray
    ) -> Tuple[float, Callable[[np.ndarray, np.ndarray], None]]:
        # initialize real
        func = c_real()

        # get pointers to arrays
        kappa_ptr = kappa.ctypes.data_as(POINTER(c_real))
        grad_ptr = grad.ctypes.data_as(POINTER(c_real))
        h_diag_ptr = h_diag.ctypes.data_as(POINTER(c_real))

        # initialize Hessian linear transformation function pointer
        hess_x_funptr = hess_x_interface_type()

        # update orbital function
        error = self.update_orbs_funptr(
            kappa_ptr,
            byref(func),
            grad_ptr,
            h_diag_ptr,
            byref(hess_x_funptr),
            self.context,
        )
        if error != 0:
            msg = "Orbital updating function raised error."
            if self._otr_exception is not None and "exc" in self._otr_exception:
                raise RuntimeError(msg) from self._otr_exception["exc"]
            raise RuntimeError(msg)

        # attach the Hessian linear transformation interface to the object so that it
        # persists in Python to ensure that it is not garbage collected when the
        # factory completes, and let it keep this object alive
        self.hess_x = HessXPyInterface(
            hess_x_funptr, self._otr_exception, self.context, _owner=self
        )

        return func.value, self.hess_x
