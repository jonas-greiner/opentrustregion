# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import annotations

import os
import sys
import warnings
from pathlib import Path
import numpy as np
from importlib import resources
from ctypes import (
    CDLL,
    CFUNCTYPE,
    POINTER,
    byref,
    string_at,
    cast,
    addressof,
    memmove,
    c_double,
    c_int64,
    c_int32,
    c_bool,
    c_void_p,
    c_char_p,
    c_char,
    Structure,
)
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Tuple, Callable, Optional, Any, Dict, Type, Union


# load the opentrustregion library, fallback to testsuite in case opentrustregion was
# statically compiled
if sys.platform == "darwin":
    ext = "dylib"
elif sys.platform == "win32":
    ext = "dll"
else:
    ext = "so"
lib_candidates = [
    f"libopentrustregion.{ext}",
    f"libopentrustregion_32.{ext}",
    f"libopentrustregion_64.{ext}",
]
if sys.platform != "win32":
    # libopentrustregion symbols not exported through libotrtestsuite on windows
    lib_candidates.append(f"libotrtestsuite.{ext}")
conda_prefix = Path("/opt/anaconda1anaconda2anaconda3")


def _load_lib() -> CDLL:
    # try to load from installed package (conda)
    for lib_name in lib_candidates:
        conda_path_unix = conda_prefix / "lib" / lib_name
        conda_path_wind = conda_prefix / "Library" / "bin" / lib_name
        if conda_path_unix.exists():
            return CDLL(str(conda_path_unix))
        elif conda_path_wind.exists():
            return CDLL(str(conda_path_wind))

    # fallback: try to load from installed package (site-packages)
    for lib_name in lib_candidates:
        lib_path = resources.files("pyopentrustregion") / lib_name
        if lib_path.is_file():
            return CDLL(str(lib_path))

    # fallback: try to load from same directory (editable install)
    for lib_name in lib_candidates:
        local_path = os.path.join(os.path.dirname(__file__), lib_name)
        if os.path.exists(local_path):
            return CDLL(local_path)

    # fallback: try to load from ../build (development build)
    for lib_name in lib_candidates:
        build_path = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "../build", lib_name)
        )
        if os.path.exists(build_path):
            return CDLL(build_path)

    # if all failed
    raise FileNotFoundError(
        f"Cannot find either opentrustregion or testsuite library ({lib_candidates})"
    )


lib = _load_lib()

# determine integer size used in library
ilp64 = c_bool.in_dll(lib, "ilp64")
c_int: Type[Union[c_int32, c_int64]]
if ilp64.value:
    c_int = c_int64
else:
    c_int = c_int32

# define real type
c_real = c_double

# fixed size strings for keywords
kw_len = c_int.in_dll(lib, "kw_len_c").value


def callback_type(restype, *argtypes):
    """
    this function builds the ctypes signature of a callback function, every callback
    function of the C interface receives the opaque host context as its last argument
    so it is appended here instead of at every declaration
    """
    return CFUNCTYPE(restype, *argtypes, c_void_p)


# callback function ctypes specifications, ctypes can only deal with simple return
# types so we interface to Fortran subroutines by creating pointers to the relevant data
hess_x_interface_type = callback_type(c_int, POINTER(c_real), POINTER(c_real))
update_orbs_interface_type = callback_type(
    c_int,
    POINTER(c_real),
    POINTER(c_real),
    POINTER(c_real),
    POINTER(c_real),
    POINTER(hess_x_interface_type),
)
obj_func_interface_type = callback_type(c_int, POINTER(c_real), POINTER(c_real))
precond_interface_type = callback_type(
    c_int, POINTER(c_real), POINTER(c_real), POINTER(c_real)
)
precond_pd_interface_type = callback_type(c_int, POINTER(c_real), POINTER(c_real))
project_interface_type = callback_type(c_int, POINTER(c_real))
get_extra_trial_vectors_interface_type = callback_type(c_int, POINTER(c_real), c_int)
conv_check_interface_type = callback_type(c_int, POINTER(c_bool))
logger_interface_type = callback_type(None, c_char_p)


def adopt_collector(*callbacks: Any) -> "Dict[str, Exception]":
    """
    this function returns the exception collector to use for a run
    """
    for callback in callbacks:
        collector = getattr(callback, "_otr_exception", None)
        if collector is not None:
            collector.clear()
            return collector

    return {}


# define interface factories
class CallbackInterface:
    """
    this class is the base class of the callback function interfaces, the C interface
    hands every callback function the opaque host context as its last argument while
    Python hosts carry their data in closures instead, so it is dropped here and the
    subclasses only implement call without it
    """

    def __call__(self, *args):
        return self.call(*args[:-1])


@dataclass
class ObjFuncInterface(CallbackInterface):
    """
    this class provides the interface for the objective function
    """

    obj_func: Callable[[np.ndarray], float]
    n_param: int
    exception: Dict[str, Exception]

    def call(self, kappa_ptr, func_ptr) -> int:
        # convert matrix pointers to numpy arrays
        kappa = np.ctypeslib.as_array(kappa_ptr, shape=(self.n_param,))

        try:
            func_ptr[0] = self.obj_func(kappa)
        except Exception as e:
            self.exception["exc"] = e
            return 1

        return 0


@dataclass
class HessXInterface(CallbackInterface):
    """
    this class provides the interface for the Hessian linear transformation
    """

    hess_x: Callable[[np.ndarray, np.ndarray], None]
    n_param: int
    exception: Dict[str, Exception]

    def call(self, x_ptr, hx_ptr) -> int:
        # convert trial vector pointer to numpy array
        x = np.ctypeslib.as_array(x_ptr, shape=(self.n_param,))
        hx = np.ctypeslib.as_array(hx_ptr, shape=(self.n_param,))

        # perform linear transformation
        try:
            self.hess_x(x, hx)
        except Exception as e:
            self.exception["exc"] = e
            return 1

        return 0


@dataclass
class UpdateOrbsInterface(CallbackInterface):
    """
    this class provides the interface to the orbital updating function
    """

    update_orbs: Callable[
        [np.ndarray, np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ]
    n_param: int
    exception: Dict[str, Exception]
    hess_x_funptr: Optional[Any] = None

    def call(self, kappa_ptr, func_ptr, grad_ptr, h_diag_ptr, hess_x_funptr) -> int:
        # convert matrix pointers to numpy arrays
        kappa = np.ctypeslib.as_array(kappa_ptr, shape=(self.n_param,))
        grad = np.ctypeslib.as_array(grad_ptr, shape=(self.n_param,))
        h_diag = np.ctypeslib.as_array(h_diag_ptr, shape=(self.n_param,))

        # update orbitals and retrieve objective function, gradient, Hessian diagonal
        # and Hessian linear transformation function
        try:
            func_ptr[0], hess_x = self.update_orbs(kappa, grad, h_diag)
        except Exception as e:
            self.exception["exc"] = e
            return 1

        # an orbital update which does not provide a Hessian linear transformation
        # leaves the pointer unset so that the solver reports the missing one, while
        # one that is not callable is reported as an error of the orbital update
        if hess_x is None:
            return 0
        if not callable(hess_x):
            self.exception["exc"] = TypeError(
                "Orbital update returned a Hessian linear transformation of type "
                f"{type(hess_x).__name__}, which is not callable."
            )
            return 1

        # attach the response interface to the object so that it persists in Python
        # to ensure that it is not garbage collected when the factory completes
        self.hess_x_funptr = hess_x_interface_type(
            HessXInterface(hess_x, self.n_param, self.exception)
        )
        hess_x_funptr[0] = self.hess_x_funptr

        return 0


@dataclass
class SettingsCallbackInterface(CallbackInterface):
    """
    this class is the base class of the interfaces to the callback functions of a
    settings object, which the settings build once per callback function and which
    call the callable currently assigned to it with the number of parameters and the
    exception collector the settings hold for the current call into the library
    """

    settings: Settings
    name: str

    @property
    def func(self) -> Callable:
        return self.settings._callables[self.name]

    @property
    def n_param(self) -> Optional[int]:
        return self.settings._n_param

    @property
    def exception(self) -> Dict[str, Exception]:
        return self.settings._exception


class PrecondInterface(SettingsCallbackInterface):
    """
    this class provides the interface to the preconditioning function
    """

    def call(self, residual_ptr, mu_ptr, precond_residual_ptr) -> int:
        try:
            # convert pointers to numpy arrays and float
            residual = np.ctypeslib.as_array(residual_ptr, shape=(self.n_param,))
            mu = mu_ptr[0]
            precond_residual = np.ctypeslib.as_array(
                precond_residual_ptr, shape=(self.n_param,)
            )

            # call preconditioner
            self.func(residual, mu, precond_residual)
        except Exception as e:
            self.exception["exc"] = e
            return 1

        return 0


class PrecondPDInterface(SettingsCallbackInterface):
    """
    this class provides the interface to the positive-definite preconditioning function
    used to define the ellipsoidal trust-region metric for the TCG and GLTR subsystem
    solvers
    """

    def call(self, residual_ptr, precond_residual_ptr) -> int:
        try:
            # convert pointers to numpy arrays
            residual = np.ctypeslib.as_array(residual_ptr, shape=(self.n_param,))
            precond_residual = np.ctypeslib.as_array(
                precond_residual_ptr, shape=(self.n_param,)
            )

            # call positive-definite preconditioner
            self.func(residual, precond_residual)
        except Exception as e:
            self.exception["exc"] = e
            return 1

        return 0


class ProjectInterface(SettingsCallbackInterface):
    """
    this class provides the interface to the projection function
    """

    def call(self, vector_ptr) -> int:
        try:
            # convert matrix pointers to numpy arrays
            vector = np.ctypeslib.as_array(vector_ptr, shape=(self.n_param,))

            # call projection function
            self.func(vector)
        except Exception as e:
            self.exception["exc"] = e
            return 1

        return 0


class ConvCheckInterface(SettingsCallbackInterface):
    """
    this class provides the interface to the convergence check function
    """

    def call(self, conv_ptr) -> int:
        # call convergence check
        try:
            conv_ptr[0] = self.func()
        except Exception as e:
            self.exception["exc"] = e
            return 1

        return 0


class GetExtraTrialVectorsInterface(SettingsCallbackInterface):
    """
    this class provides the interface to the extra trial vector function
    """

    def call(self, trial_vectors_ptr, n_extra_trial_vectors) -> int:
        try:
            # convert matrix pointers to numpy arrays
            trial_vectors = np.ctypeslib.as_array(
                trial_vectors_ptr, shape=(n_extra_trial_vectors, self.n_param)
            )

            # call extra trial vector function
            self.func(trial_vectors)
        except Exception as e:
            self.exception["exc"] = e
            return 1

        return 0


class SettingsHessXInterface(SettingsCallbackInterface):
    """
    this class provides the interface to a Hessian linear transformation function of
    the settings
    """

    def call(self, x_ptr, hx_ptr) -> int:
        try:
            # convert trial vector pointer to numpy array
            x = np.ctypeslib.as_array(x_ptr, shape=(self.n_param,))
            hx = np.ctypeslib.as_array(hx_ptr, shape=(self.n_param,))

            # perform linear transformation
            self.func(x, hx)
        except Exception as e:
            self.exception["exc"] = e
            return 1

        return 0


class LoggerInterface(SettingsCallbackInterface):
    """
    this class provides the interface to the logging function
    """

    def call(self, message):
        # call logger
        try:
            self.func(string_at(message).decode("utf-8"))
        except Exception as e:
            # the logger has no error channel since the result does not depend on it,
            # so the first failure is kept and reported once the call into the library
            # has returned
            self.exception.setdefault("logger", e)


def raise_callback_error(error: int, exception: Optional[Dict[str, Exception]], msg):
    """
    this function raises if a C callback function called from Python failed, with the
    exception raised by a failing Python callback function it reached as its cause
    """
    if error != 0:
        if exception is not None and "exc" in exception:
            raise RuntimeError(msg) from exception["exc"]
        raise RuntimeError(msg)


def c_buffer(array: np.ndarray) -> np.ndarray:
    """
    this function returns an array a C function can read and write in place, the given
    array itself if it is a writeable, C-contiguous float64 array, and a C-contiguous
    float64 copy otherwise, into which the values of an output have to be written back
    """
    if array.flags.writeable:
        return np.ascontiguousarray(array, dtype=np.float64)
    return np.array(array, dtype=np.float64, order="C")


def write_back(array: np.ndarray, buffer: np.ndarray):
    """
    this function writes the values a C function wrote into the buffer of an output
    back into the given array, if the buffer is a copy
    """
    if buffer is not array:
        array[...] = buffer


@dataclass
class CFunctionInterface:
    """
    this class is the base class of the Python interfaces to C callback functions,
    such as those a factory returns or wires into settings, which are called with the
    context they belong to and keep the object they rely on alive
    """

    funptr: Any
    _otr_exception: Optional[Dict[str, Exception]] = None
    context: Optional[int] = None
    _owner: Any = None


class HessXPyInterface(CFunctionInterface):
    """
    this class provides the Python interface to a Hessian linear transformation C
    function
    """

    def __call__(self, x: np.ndarray, hess_x: np.ndarray):
        x_buffer, hess_x_buffer = c_buffer(x), c_buffer(hess_x)
        error = self.funptr(
            x_buffer.ctypes.data_as(POINTER(c_real)),
            hess_x_buffer.ctypes.data_as(POINTER(c_real)),
            self.context,
        )
        write_back(hess_x, hess_x_buffer)
        raise_callback_error(
            error,
            self._otr_exception,
            "Hessian linear transformation function raised error.",
        )


class PrecondPyInterface(CFunctionInterface):
    """
    this class provides the Python interface to a level-shifted preconditioner C
    function
    """

    def __call__(self, residual: np.ndarray, mu: float, precond_residual: np.ndarray):
        residual_buffer = c_buffer(residual)
        precond_residual_buffer = c_buffer(precond_residual)
        error = self.funptr(
            residual_buffer.ctypes.data_as(POINTER(c_real)),
            c_real(mu),
            precond_residual_buffer.ctypes.data_as(POINTER(c_real)),
            self.context,
        )
        write_back(precond_residual, precond_residual_buffer)
        raise_callback_error(
            error, self._otr_exception, "Preconditioner function raised error."
        )


class PrecondPDPyInterface(CFunctionInterface):
    """
    this class provides the Python interface to a positive-definite preconditioner C
    function
    """

    def __call__(self, residual: np.ndarray, precond_residual: np.ndarray):
        residual_buffer = c_buffer(residual)
        precond_residual_buffer = c_buffer(precond_residual)
        error = self.funptr(
            residual_buffer.ctypes.data_as(POINTER(c_real)),
            precond_residual_buffer.ctypes.data_as(POINTER(c_real)),
            self.context,
        )
        write_back(precond_residual, precond_residual_buffer)
        raise_callback_error(
            error,
            self._otr_exception,
            "Positive-definite preconditioner function raised error.",
        )


class ProjectPyInterface(CFunctionInterface):
    """
    this class provides the Python interface to a projection C function
    """

    def __call__(self, vector: np.ndarray):
        vector_buffer = c_buffer(vector)
        error = self.funptr(vector_buffer.ctypes.data_as(POINTER(c_real)), self.context)
        write_back(vector, vector_buffer)
        raise_callback_error(
            error, self._otr_exception, "Projection function raised error."
        )


class GetExtraTrialVectorsPyInterface(CFunctionInterface):
    """
    this class provides the Python interface to an extra trial vector C function
    """

    def __call__(self, trial_vectors: np.ndarray):
        trial_vectors_buffer = c_buffer(trial_vectors)
        error = self.funptr(
            trial_vectors_buffer.ctypes.data_as(POINTER(c_real)),
            trial_vectors.shape[0],
            self.context,
        )
        write_back(trial_vectors, trial_vectors_buffer)
        raise_callback_error(
            error, self._otr_exception, "Extra trial vector function raised error."
        )


class ConvCheckPyInterface(CFunctionInterface):
    """
    this class provides the Python interface to a convergence check C function
    """

    def __call__(self) -> bool:
        converged = c_bool()
        error = self.funptr(byref(converged), self.context)
        raise_callback_error(
            error, self._otr_exception, "Convergence check function raised error."
        )
        return converged.value


class LoggerPyInterface(CFunctionInterface):
    """
    this class provides the Python interface to a logging C function
    """

    def __call__(self, message: str):
        self.funptr(message.encode("utf-8"), self.context)


# define classes corresponding to C structs for settings
class StabilitySettingsC(Structure):
    _fields_ = [
        ("precond", c_void_p),
        ("project", c_void_p),
        ("approx_hess_x", c_void_p),
        ("get_extra_trial_vectors", c_void_p),
        ("conv_check", c_void_p),
        ("logger", c_void_p),
        ("hess_symm", c_bool),
        ("stop_on_instability", c_bool),
        ("initialized", c_bool),
        ("conv_tol", c_real),
        ("stability_thresh", c_real),
        ("n_random_trial_vectors", c_int),
        ("n_extra_trial_vectors", c_int),
        ("n_iter", c_int),
        ("jacobi_davidson_start", c_int),
        ("seed", c_int),
        ("verbose", c_int),
        ("diag_solver", c_char * (kw_len + 1)),
        ("context", c_void_p),
    ]


class SolverSettingsC(Structure):
    _fields_ = [
        ("precond", c_void_p),
        ("precond_pd", c_void_p),
        ("project", c_void_p),
        ("get_extra_trial_vectors", c_void_p),
        ("conv_check", c_void_p),
        ("stability_hess_x", c_void_p),
        ("logger", c_void_p),
        ("stability", c_bool),
        ("line_search", c_bool),
        ("refresh_hess", c_bool),
        ("grad_noise", c_bool),
        ("hess_symm", c_bool),
        ("initialized", c_bool),
        ("conv_tol", c_real),
        ("start_trust_radius", c_real),
        ("global_red_factor", c_real),
        ("local_red_factor", c_real),
        ("n_random_trial_vectors", c_int),
        ("n_extra_trial_vectors", c_int),
        ("n_macro", c_int),
        ("n_micro", c_int),
        ("jacobi_davidson_start", c_int),
        ("seed", c_int),
        ("verbose", c_int),
        ("subsystem_solver", c_char * (kw_len + 1)),
        ("trust_region_shape", c_char * (kw_len + 1)),
        ("stability_settings", StabilitySettingsC),
        ("context", c_void_p),
    ]


# define classes corresponding to C structs for states
class StabilityStateC(Structure):
    _fields_ = [
        ("min_eigval", c_real),
        ("residual_rms", c_real),
        ("iteration", c_int),
        ("n_hess_x", c_int),
    ]


class SolverStateC(Structure):
    _fields_ = [
        ("max_precision_reached", c_bool),
        ("func", c_real),
        ("grad_rms", c_real),
        ("trust_radius", c_real),
        ("relative_residual", c_real),
        ("level_shift", c_real),
        ("step_size", c_real),
        ("macro_iteration", c_int),
        ("micro_iteration", c_int),
        ("jacobi_davidson_iteration", c_int),
        ("n_update_orbs", c_int),
        ("n_obj_func", c_int),
        ("n_hess_x", c_int),
        ("stability_state", StabilityStateC),
    ]


# define classes for settings
class Settings:

    c_struct: type[Structure]
    init_c_struct: Any

    # the interface, the C function type and the Python interface to a C function of
    # every callback function of the settings
    callback_fields: Dict[str, Tuple[type, Any, type]] = {}

    def __init__(self, settings_c: Optional[Structure] = None):
        """
        this function initializes the settings class
        """
        # specify C interface signature
        self.init_c_struct.argtypes = [POINTER(self.c_struct)]
        self.init_c_struct.restype = None

        # call C-side initialization to populate defaults, use passed settings if
        # provided
        if settings_c is None:
            self.settings_c = self.c_struct()
            self.init_c_struct(byref(self.settings_c))
        else:
            self.settings_c = settings_c

        # the callables assigned to the callback functions, and the number of
        # parameters and the exception collector of the current call into the library,
        # which the interfaces of the callback functions read when they are called
        self._callables: Dict[str, Optional[Callable]] = {}
        self._n_param: Optional[int] = None
        self._exception: Dict[str, Exception] = {}

        # build the C callback function of every callback function once, which calls
        # the callable currently assigned to it, so that the C settings only change
        # when a callable is assigned and a factory can keep them like those of a C
        # host
        self._callbacks = {
            name: interface_type(interface(self, name))
            for name, (interface, interface_type, _) in self.callback_fields.items()
        }
        self._callback_addresses = {
            name: cast(callback, c_void_p).value
            for name, callback in self._callbacks.items()
        }

    def prepare_call(self, n_param: Optional[int], exception: Dict[str, Exception]):
        """
        this function hands the callback functions the number of parameters and the
        exception collector of a call into the library
        """
        self._n_param = n_param
        self._exception = exception


class SolverSettings(Settings):

    c_struct = SolverSettingsC
    init_c_struct = lib.init_solver_settings
    callback_fields = {
        "precond": (PrecondInterface, precond_interface_type, PrecondPyInterface),
        "precond_pd": (
            PrecondPDInterface,
            precond_pd_interface_type,
            PrecondPDPyInterface,
        ),
        "project": (ProjectInterface, project_interface_type, ProjectPyInterface),
        "get_extra_trial_vectors": (
            GetExtraTrialVectorsInterface,
            get_extra_trial_vectors_interface_type,
            GetExtraTrialVectorsPyInterface,
        ),
        "conv_check": (
            ConvCheckInterface,
            conv_check_interface_type,
            ConvCheckPyInterface,
        ),
        "stability_hess_x": (
            SettingsHessXInterface,
            hess_x_interface_type,
            HessXPyInterface,
        ),
        "logger": (LoggerInterface, logger_interface_type, LoggerPyInterface),
    }

    precond: Optional[Callable[[np.ndarray, float, np.ndarray], None]]
    precond_pd: Optional[Callable[[np.ndarray, np.ndarray], None]]
    project: Optional[Callable[[np.ndarray], None]]
    get_extra_trial_vectors: Optional[Callable[[np.ndarray], None]]
    conv_check: Optional[Callable[[], bool]]
    stability_hess_x: Optional[Callable[[np.ndarray, np.ndarray], None]]
    logger: Optional[Callable[[str], None]]

    def __init__(self):
        super().__init__()

        # the nested stability check settings keep these settings alive, which own the
        # extension objects wired into both and the callback functions the nested
        # settings inherit, so that the nested settings and the callables read from
        # them stay valid once nothing else refers to these settings
        self._stability_settings = StabilitySettings(
            settings_c=self.settings_c.stability_settings
        )
        self._stability_settings._parent = self

    @property
    def stability_settings(self) -> StabilitySettings:
        return self._stability_settings

    def prepare_call(self, n_param: Optional[int], exception: Dict[str, Exception]):
        """
        this function hands the callback functions of the settings and of their
        nested stability check settings the number of parameters and the exception
        collector of a call into the library
        """
        super().prepare_call(n_param, exception)
        self.stability_settings.prepare_call(n_param, exception)


class StabilitySettings(Settings):

    c_struct = StabilitySettingsC
    init_c_struct = lib.init_stability_settings
    callback_fields = {
        "precond": (PrecondInterface, precond_interface_type, PrecondPyInterface),
        "project": (ProjectInterface, project_interface_type, ProjectPyInterface),
        "approx_hess_x": (
            SettingsHessXInterface,
            hess_x_interface_type,
            HessXPyInterface,
        ),
        "get_extra_trial_vectors": (
            GetExtraTrialVectorsInterface,
            get_extra_trial_vectors_interface_type,
            GetExtraTrialVectorsPyInterface,
        ),
        "conv_check": (
            ConvCheckInterface,
            conv_check_interface_type,
            ConvCheckPyInterface,
        ),
        "logger": (LoggerInterface, logger_interface_type, LoggerPyInterface),
    }

    precond: Optional[Callable[[np.ndarray, float, np.ndarray], None]]
    project: Optional[Callable[[np.ndarray], None]]
    approx_hess_x: Optional[Callable[[np.ndarray, np.ndarray], None]]
    get_extra_trial_vectors: Optional[Callable[[np.ndarray], None]]
    conv_check: Optional[Callable[[], bool]]
    logger: Optional[Callable[[str], None]]


def bind_callback_fields(cls: type[Settings]):
    """
    this function generates properties for the callback functions of a settings class:
    assigning a callable writes the C callback function the settings built for it into
    the C settings, or a null pointer for None, and reading returns the callable while
    the C settings hold that C callback function, a Python interface to the C function
    a factory wired into them otherwise, called with the context of the settings and
    keeping them alive, and None for an unset callback function
    """
    for name, (_, interface_type, py_interface) in cls.callback_fields.items():

        def make_property(name, interface_type, py_interface):

            def getter(self):
                funptr = getattr(self.settings_c, name)
                if not funptr:
                    return None
                if funptr == self._callback_addresses[name]:
                    return self._callables.get(name)
                return py_interface(
                    interface_type(funptr),
                    _otr_exception=self._exception,
                    context=self.settings_c.context,
                    _owner=self,
                )

            def setter(self, func):
                self._callables[name] = func
                setattr(
                    self.settings_c,
                    name,
                    None if func is None else self._callback_addresses[name],
                )

            return property(getter, setter)

        setattr(cls, name, make_property(name, interface_type, py_interface))


def auto_bind_fields(cls: type[Settings]):
    """
    this function automatically generates properties for non-pointer fields in a
    settings_c object that is an attribute of cls
    """
    for field_info in cls.c_struct._fields_:
        field_name, field_type = field_info[:2]

        # skip if function pointer (these will be initialized separately) or nested
        # structures
        if field_type is c_void_p or (
            isinstance(field_type, type) and issubclass(field_type, Structure)
        ):
            continue

        # character arrays need to be handled separately
        elif field_type is c_char * (kw_len + 1):

            def make_property(name):

                def getter(self):
                    # get bytes from ctypes Structure
                    raw = getattr(self.settings_c, name)

                    # only get elements before null bytes and convert to strings
                    return raw.split(b"\0", 1)[0].decode("utf-8")

                def setter(self, value):
                    # check for type (will be a string if supplied by a user and a
                    # bytes when extracted from a ctypes Structure)
                    if isinstance(value, str):
                        b = value.encode("utf-8")[:kw_len]
                    elif isinstance(value, bytes):
                        b = value[:kw_len]

                    # pad with null bytes up to kw_len
                    b = b.ljust(kw_len, b"\0")

                    # copy to memory in ctypes Structure
                    offset = getattr(type(self.settings_c), name).offset
                    memmove(addressof(self.settings_c) + offset, b, kw_len)

                return property(getter, setter)

        else:

            def make_property(name):

                def getter(self):
                    return getattr(self.settings_c, name)

                def setter(self, value):
                    setattr(self.settings_c, name, value)

                return property(getter, setter)

        # capture field_name in current loop scope
        setattr(cls, field_name, make_property(field_name))


# ensure that appropriate fields are automatically set in settings_c object
auto_bind_fields(SolverSettings)
auto_bind_fields(StabilitySettings)
bind_callback_fields(SolverSettings)
bind_callback_fields(StabilitySettings)


# define classes for states, which the library keeps current during a call and which
# can therefore be read by the callback functions
class State:

    c_struct: type[Structure]

    def __init__(self, state_c: Optional[Structure] = None):
        """
        this function initializes the state class on a new C struct, which holds the
        default values, or on a given one
        """
        self.state_c = self.c_struct() if state_c is None else state_c


class SolverState(State):

    c_struct = SolverStateC

    def __init__(self):
        super().__init__()
        self._stability_state = StabilityState(state_c=self.state_c.stability_state)

    @property
    def stability_state(self) -> StabilityState:
        return self._stability_state


class StabilityState(State):

    c_struct = StabilityStateC


def bind_state_fields(cls: type[State]):
    """
    this function generates read-only properties for the fields of the C struct of a
    state class, apart from the nested state
    """
    for field_info in cls.c_struct._fields_:
        field_name, field_type = field_info[:2]

        # skip nested structures
        if isinstance(field_type, type) and issubclass(field_type, Structure):
            continue

        def make_property(name):

            def getter(self):
                return getattr(self.state_c, name)

            return property(getter)

        # capture field_name in current loop scope
        setattr(cls, field_name, make_property(field_name))


# ensure that the fields of the state_c object can be read
bind_state_fields(SolverState)
bind_state_fields(StabilityState)


def raise_on_failure(error: int, exception: Dict[str, Exception], name: str):
    """
    this function raises if a call into the library failed, with the exception raised
    by a failing callback function as its cause, a failed logger cannot affect the
    result and is therefore only reported as a warning, both before an error and
    instead of one
    """
    if "logger" in exception:
        warning = RuntimeWarning(
            f"OpenTrustRegion {name} called a logging function that raised "
            f"{exception['logger']!r}, so log messages from that point on may be "
            "missing."
        )
        warning.__cause__ = exception["logger"]
        warnings.warn(warning, stacklevel=3)
    if error:
        raise RuntimeError(
            f"OpenTrustRegion {name} produced error (code {error})."
        ) from exception.get("exc")


def solver(
    update_orbs: Callable[
        [np.ndarray, np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ],
    obj_func: Callable[[np.ndarray], float],
    n_param: int,
    settings: Optional[SolverSettings] = None,
    state: Optional[SolverState] = None,
):
    # collector for exceptions raised inside the callback interfaces, adopted from a
    # factory-produced callback when there is one so its exceptions are chained too
    exception = adopt_collector(update_orbs, obj_func)

    # define interfaces for callback functions
    obj_func_interface = obj_func_interface_type(
        ObjFuncInterface(obj_func, n_param, exception)
    )
    update_orbs_interface = update_orbs_interface_type(
        UpdateOrbsInterface(update_orbs, n_param, exception)
    )

    # hand the callback functions of the settings the number of parameters and the
    # exception collector of the call
    if settings is not None:
        settings.prepare_call(n_param, exception)

    # define result and argument types
    lib.solver.restype = c_int
    lib.solver.argtypes = [
        update_orbs_interface_type,
        obj_func_interface_type,
        c_int,
        POINTER(SolverSettingsC),
        POINTER(SolverStateC),
    ]

    # call Fortran function
    error = lib.solver(
        update_orbs_interface,
        obj_func_interface,
        n_param,
        byref(settings.settings_c) if settings is not None else None,
        byref(state.state_c) if state is not None else None,
    )

    raise_on_failure(error, exception, "solver")


def stability_check(
    h_diag: np.ndarray,
    hess_x: Callable[[np.ndarray, np.ndarray], None],
    n_param: int,
    settings: Optional[StabilitySettings] = None,
    state: Optional[StabilityState] = None,
    kappa: Optional[np.ndarray] = None,
) -> bool:
    # collector for exceptions raised inside the callback interfaces, adopted from a
    # factory-produced callback when there is one so its exceptions are chained too
    exception = adopt_collector(hess_x)

    # define interfaces for callback functions
    hess_x_interface = hess_x_interface_type(HessXInterface(hess_x, n_param, exception))

    # hand the callback functions of the settings the number of parameters and the
    # exception collector of the call
    if settings is not None:
        settings.prepare_call(n_param, exception)

    # define result and argument types
    lib.stability_check.restype = c_int
    lib.stability_check.argtypes = [
        POINTER(c_real),
        hess_x_interface_type,
        c_int,
        POINTER(c_bool),
        POINTER(StabilitySettingsC),
        POINTER(StabilityStateC),
        POINTER(c_real),
    ]

    # initialize return variable
    stable = c_bool(False)

    # call Fortran function
    error = lib.stability_check(
        h_diag.ctypes.data_as(POINTER(c_real)),
        hess_x_interface,
        n_param,
        byref(stable),
        byref(settings.settings_c) if settings is not None else None,
        byref(state.state_c) if state is not None else None,
        kappa.ctypes.data_as(POINTER(c_real)) if kappa is not None else None,
    )

    raise_on_failure(error, exception, "stability check")

    return stable.value
