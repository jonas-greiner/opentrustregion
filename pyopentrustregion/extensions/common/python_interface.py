# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import annotations

import numpy as np
from ctypes import CFUNCTYPE, POINTER, byref, c_void_p, cast
from dataclasses import dataclass
from inspect import signature
from typing import TYPE_CHECKING
from pyopentrustregion.python_interface import (
    c_int,
    c_real,
    hess_x_interface_type,
    precond_interface_type,
    precond_pd_interface_type,
    project_interface_type,
    get_extra_trial_vectors_interface_type,
    SolverSettings,
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
get_response_interface_type = CFUNCTYPE(c_int, POINTER(c_real), POINTER(c_real))
evaluate_dm_interface_type = CFUNCTYPE(
    c_int,
    POINTER(c_real),
    POINTER(c_real),
    POINTER(c_real),
    POINTER(get_response_interface_type),
)


# define interface factories
@dataclass
class GetResponseInterface:
    """
    this class provides the interface to write the response to the memory provided for
    a given density matrix
    """

    get_response: GetResponseType
    n_ao: int
    n_particle: int
    closed_shell: bool
    exception: Dict[str, Exception]

    def __call__(self, dm_ao_ptr, response_ptr) -> int:
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
class EvaluateDMInterface:
    """
    this class provides the interface to the density matrix evaluating function
    """

    evaluate_dm: EvaluateDMType
    n_ao: int
    n_particle: int
    closed_shell: bool
    exception: Dict[str, Exception]
    get_response_funptr: Optional[Any] = None

    def __call__(self, dm_ao_ptr, energy_ptr, fock_ptr, get_response_funptr) -> int:
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
    this class provides the Python interface to the objective function,
    evaluate_dm_interface is stored to ensure that it is not garbage collected when the
    factory completes
    """

    obj_func_funptr: Any
    evaluate_dm_interface: Any
    _otr_exception: Optional[Dict[str, Exception]] = None

    def __call__(self, kappa: np.ndarray) -> float:
        # initialize real
        func = c_real()

        # get pointers to arrays
        kappa_ptr = kappa.ctypes.data_as(POINTER(c_real))

        # objective function
        error = self.obj_func_funptr(kappa_ptr, byref(func))
        if error != 0:
            if self._otr_exception is not None and "exc" in self._otr_exception:
                raise RuntimeError(
                    "Objective function raised error."
                ) from self._otr_exception["exc"]
            raise RuntimeError("Objective function raised error.")

        return func.value


@dataclass
class ProjectPyInterface:
    """
    this class provides the Python interface to the projection function
    """

    project_funptr: Any
    _otr_exception: Optional[Dict[str, Exception]] = None

    def __call__(self, vector: np.ndarray):
        # get pointers to arrays
        vector_ptr = vector.ctypes.data_as(POINTER(c_real))

        # projection function
        error = self.project_funptr(vector_ptr)
        if error != 0:
            if self._otr_exception is not None and "exc" in self._otr_exception:
                raise RuntimeError(
                    "Projection function raised error."
                ) from self._otr_exception["exc"]
            raise RuntimeError("Projection function raised error.")

        return


@dataclass
class GetExtraTrialVectorsPyInterface:
    """
    this class provides the Python interface to the extra trial vector function
    """

    get_extra_trial_vectors_funptr: Any
    _otr_exception: Optional[Dict[str, Exception]] = None

    def __call__(self, trial_vectors: np.ndarray):
        # get pointers to arrays
        trial_vectors_ptr = trial_vectors.ctypes.data_as(POINTER(c_real))

        # extra trial vector function
        error = self.get_extra_trial_vectors_funptr(
            trial_vectors_ptr, trial_vectors.shape[0]
        )
        if error != 0:
            if self._otr_exception is not None and "exc" in self._otr_exception:
                raise RuntimeError(
                    "Extra trial vector function raised error."
                ) from self._otr_exception["exc"]
            raise RuntimeError("Extra trial vector function raised error.")

        return


@dataclass
class PrecondPyInterface:
    """
    this class provides the Python interface to the level-shifted preconditioner
    function
    """

    precond_funptr: Any
    _otr_exception: Optional[Dict[str, Exception]] = None

    def __call__(self, residual: np.ndarray, mu: float, precond_residual: np.ndarray):
        # get pointers to arrays
        residual_ptr = residual.ctypes.data_as(POINTER(c_real))
        precond_residual_ptr = precond_residual.ctypes.data_as(POINTER(c_real))

        # call preconditioner function
        error = self.precond_funptr(residual_ptr, c_real(mu), precond_residual_ptr)
        if error != 0:
            if self._otr_exception is not None and "exc" in self._otr_exception:
                raise RuntimeError(
                    "Preconditioner function raised error."
                ) from self._otr_exception["exc"]
            raise RuntimeError("Preconditioner function raised error.")

        return


@dataclass
class PrecondPDPyInterface:
    """
    this class provides the Python interface to the positive-definite preconditioner
    function
    """

    precond_pd_funptr: Any
    _otr_exception: Optional[Dict[str, Exception]] = None

    def __call__(self, residual: np.ndarray, precond_residual: np.ndarray):
        # get pointers to arrays
        residual_ptr = residual.ctypes.data_as(POINTER(c_real))
        precond_residual_ptr = precond_residual.ctypes.data_as(POINTER(c_real))

        # call preconditioner function
        error = self.precond_pd_funptr(residual_ptr, precond_residual_ptr)
        if error != 0:
            msg = "Positive-definite preconditioner function raised error."
            if self._otr_exception is not None and "exc" in self._otr_exception:
                raise RuntimeError(msg) from self._otr_exception["exc"]
            raise RuntimeError(msg)

        return


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


def attach_wired_callbacks(
    solver_settings: SolverSettings, exception: Dict[str, Exception]
) -> None:
    """
    this function attaches the routines a factory has wired into the C solver
    settings, skipping the C callbacks the solver built from the Python functions of
    the settings, which the settings already hold and which a wrapper around them would
    outlive once the solver replaces them
    """
    wrappers = {
        "precond": (precond_interface_type, PrecondPyInterface),
        "precond_pd": (precond_pd_interface_type, PrecondPDPyInterface),
        "project": (project_interface_type, ProjectPyInterface),
        "get_extra_trial_vectors": (
            get_extra_trial_vectors_interface_type,
            GetExtraTrialVectorsPyInterface,
        ),
    }
    for settings in (solver_settings, solver_settings.stability_settings):
        for name, (interface_type, py_interface) in wrappers.items():
            if not hasattr(settings.settings_c, name):
                continue
            settings_c_funptr = getattr(settings.settings_c, name)
            own_interface = getattr(settings.settings_c, name + "_interface", None)
            if not settings_c_funptr or (
                own_interface is not None
                and settings_c_funptr == cast(own_interface, c_void_p).value
            ):
                continue
            setattr(
                settings,
                name,
                py_interface(
                    interface_type(settings_c_funptr), _otr_exception=exception
                ),
            )


@dataclass
class HessXPyInterface:
    """
    this class provides the Python interface to the Hessian linear transformation
    function
    """

    hess_x_funptr: Any
    _otr_exception: Optional[Dict[str, Exception]] = None

    def __call__(self, x: np.ndarray, hess_x: np.ndarray):
        # get pointers to arrays
        x_ptr = x.ctypes.data_as(POINTER(c_real))
        hess_x_ptr = hess_x.ctypes.data_as(POINTER(c_real))

        # update orbital function
        error = self.hess_x_funptr(x_ptr, hess_x_ptr)
        if error != 0:
            msg = "Hessian linear transformation function raised error."
            if self._otr_exception is not None and "exc" in self._otr_exception:
                raise RuntimeError(msg) from self._otr_exception["exc"]
            raise RuntimeError(msg)


@dataclass
class UpdateOrbsPyInterface:
    """
    this class provides the Python interface to the orbital updating function, the
    callback interfaces are stored to ensure that they are not garbage collected when
    the factory completes
    """

    update_orbs_funptr: Any
    saved_objects: Dict[str, Any]
    hess_x: Optional[Callable[[np.ndarray, np.ndarray], None]] = None
    # collector shared with the wrapped user callback, so that solver can adopt it and
    # chain exceptions raised inside that callback; None when this object does not wrap
    # a user callback of its own
    _otr_exception: Optional[Dict[str, Exception]] = None

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
        )
        if error != 0:
            msg = "Orbital updating function raised error."
            if self._otr_exception is not None and "exc" in self._otr_exception:
                raise RuntimeError(msg) from self._otr_exception["exc"]
            raise RuntimeError(msg)

        # attach the Hessian linear transformation interface to the object so that it
        # persists in Python to ensure that it is not garbage collected when the
        # factory completes
        self.hess_x = HessXPyInterface(hess_x_funptr, self._otr_exception)

        return func.value, self.hess_x
