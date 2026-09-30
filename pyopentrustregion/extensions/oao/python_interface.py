# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import annotations

import operator
import numpy as np
from ctypes import CFUNCTYPE, POINTER, byref, c_bool, c_void_p, cast, Structure
from dataclasses import dataclass
from typing import TYPE_CHECKING
from pyopentrustregion.python_interface import (
    lib,
    c_int,
    c_real,
    obj_func_interface_type,
    update_orbs_interface_type,
    precond_interface_type,
    precond_pd_interface_type,
    project_interface_type,
    get_extra_trial_vectors_interface_type,
    logger_interface_type,
    LoggerInterface,
    adopt_collector,
    Settings,
    SolverSettings,
    SolverSettingsC,
    auto_bind_fields,
)
from pyopentrustregion.extensions.common.python_interface import UpdateOrbsPyInterface

if TYPE_CHECKING:
    from typing import Tuple, Callable, Optional, Any, Dict

    GetResponseType = Callable[[np.ndarray, np.ndarray], None]
    EvaluateDMType = Callable[
        [np.ndarray, Optional[np.ndarray], bool],
        Tuple[float, Optional[GetResponseType]],
    ]


# callback function ctypes specifications, ctypes can only deal with simple return
# types so we interface to Fortran subroutines by creating pointers to the relevant
# data
get_response_interface_type = CFUNCTYPE(c_int, POINTER(c_real), POINTER(c_real))
evaluate_dm_interface_type = CFUNCTYPE(
    c_int,
    POINTER(c_real),
    POINTER(c_real),
    POINTER(c_real),
    POINTER(get_response_interface_type),
)


# define classes corresponding to C structs for settings
class OAOSettingsC(Structure):
    _fields_ = [
        ("logger", c_void_p),
        ("initialized", c_bool),
        ("verbose", c_int),
    ]


class OAOSettings(Settings):

    c_struct = OAOSettingsC
    try:
        init_c_struct = lib.init_oao_settings
    except AttributeError:
        raise AttributeError(
            "Please reinstall the package with: "
            "CMAKE_FLAGS='-DENABLE_OAO=ON' pip install ."
        )

    logger: Optional[Callable[[str], None]]
    logger_interface: Any


# ensure that appropriate fields are automatically set in settings_c object
auto_bind_fields(OAOSettings)


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


def check_dm_ao(
    dm_ao: np.ndarray, ao_overlap: np.ndarray, n_particle: int, n_ao: int
) -> np.ndarray:
    """
    this function checks that an AO density matrix has the given dimensions, as an
    (n_ao, n_ao) array for the closed-shell case and as a (2, n_ao, n_ao) array for the
    open-shell case, and that it can be updated in place, since the factories keep
    pointing at it, and returns the AO overlap matrix of these dimensions as a
    contiguous array
    """
    n_particle, n_ao = operator.index(n_particle), operator.index(n_ao)
    shape = (n_ao, n_ao) if n_particle == 1 else (n_particle, n_ao, n_ao)
    if dm_ao.shape != shape:
        raise ValueError(
            f"The AO density matrix has to be of shape {shape} for {n_particle} "
            f"particle(s) and {n_ao} AOs, got shape {dm_ao.shape}."
        )
    if (
        dm_ao.dtype != np.float64
        or not dm_ao.flags.c_contiguous
        or not dm_ao.flags.writeable
    ):
        raise ValueError(
            "The AO density matrix has to be a writeable and C-contiguous float64 "
            "array, since it is updated in place."
        )

    # the AO overlap matrix is symmetric, so it needs no transposition
    ao_overlap = np.ascontiguousarray(ao_overlap, dtype=np.float64)
    if ao_overlap.shape != (n_ao, n_ao):
        raise ValueError(
            f"The AO overlap matrix has to be of shape ({n_ao}, {n_ao}), got shape "
            f"{ao_overlap.shape}."
        )

    return ao_overlap


def oao_factory(
    dm_ao: np.ndarray,
    ao_overlap: np.ndarray,
    n_particle: int,
    n_ao: int,
    evaluate_dm: EvaluateDMType,
    solver_settings: SolverSettings,
    settings: OAOSettings,
) -> Tuple[
    Callable[[np.ndarray], float],
    Callable[
        [np.ndarray, np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ],
]:
    # check the density and overlap matrices against the dimensions and determine if
    # closed-shell or open-shell formalism is used
    ao_overlap = check_dm_ao(dm_ao, ao_overlap, n_particle, n_ao)
    closed_shell = n_particle == 1

    # get pointers to arrays
    dm_ao_ptr = dm_ao.ctypes.data_as(POINTER(c_real))
    ao_overlap_ptr = ao_overlap.ctypes.data_as(POINTER(c_real))

    # collector for exceptions raised inside the wrapped user callbacks; adopted from
    # evaluate_dm when it is itself factory-produced so a whole chain of factories
    # shares one, and handed to solver through the returned object
    exception = adopt_collector(evaluate_dm)

    # define interfaces for callback functions
    evaluate_dm_interface = evaluate_dm_interface_type(
        EvaluateDMInterface(evaluate_dm, n_ao, n_particle, closed_shell, exception)
    )

    # set interfaces for optional callback functions, these need to be set here since
    # the interface might need parameters that are not known when the attribute to
    # settings is set (e.g. n_param)
    settings.set_optional_callback(
        "logger", settings.logger, LoggerInterface, logger_interface_type
    )

    if not hasattr(lib, "oao_factory"):
        raise RuntimeError(
            "Please reinstall the package with: "
            "CMAKE_FLAGS='-DENABLE_OAO=ON' pip install ."
        )

    # define result and argument types
    lib.oao_factory.restype = c_int
    lib.oao_factory.argtypes = [
        POINTER(c_real),
        POINTER(c_real),
        c_int,
        c_int,
        evaluate_dm_interface_type,
        POINTER(obj_func_interface_type),
        POINTER(update_orbs_interface_type),
        POINTER(SolverSettingsC),
        POINTER(OAOSettingsC),
    ]

    # call Fortran function
    obj_func_oao_funptr = obj_func_interface_type()
    update_orbs_oao_funptr = update_orbs_interface_type()
    error = lib.oao_factory(
        dm_ao_ptr,
        ao_overlap_ptr,
        n_particle,
        n_ao,
        evaluate_dm_interface,
        byref(obj_func_oao_funptr),
        byref(update_orbs_oao_funptr),
        byref(solver_settings.settings_c),
        byref(settings.settings_c),
    )

    if error:
        if "exc" in exception:
            raise RuntimeError(
                f"OpenTrustRegion OAO factory produced error (code {error})."
            ) from exception["exc"]
        else:
            raise RuntimeError(
                f"OpenTrustRegion OAO factory produced error (code {error})."
            )

    # attach the routines the factory has wired into the solver settings
    attach_wired_callbacks(solver_settings, exception)

    return (
        ObjFuncPyInterface(
            obj_func_funptr=obj_func_oao_funptr,
            evaluate_dm_interface=evaluate_dm_interface,
            _otr_exception=exception,
        ),
        UpdateOrbsPyInterface(
            update_orbs_funptr=update_orbs_oao_funptr,
            _otr_exception=exception,
            saved_objects={
                "evaluate_dm_interface": evaluate_dm_interface,
            },
        ),
    )


def oao_deconstructor():
    if not hasattr(lib, "oao_deconstructor"):
        raise RuntimeError(
            "Please reinstall the package with: "
            "CMAKE_FLAGS='-DENABLE_OAO=ON' pip install ."
        )

    # define result and argument types
    lib.oao_deconstructor.restype = None
    lib.oao_deconstructor.argtypes = []

    # call Fortran function
    lib.oao_deconstructor()

    return
