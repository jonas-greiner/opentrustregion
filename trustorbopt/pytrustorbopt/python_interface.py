from __future__ import annotations


import sys
import numpy as np
from importlib import resources
from ctypes import CDLL, CFUNCTYPE, POINTER, byref, c_double, c_long, c_bool
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Tuple, Callable, Optional


# import dynamic library and define result and argument types
ext = "dylib" if sys.platform == "darwin" else "so"
with resources.path("pytrustorbopt", "libtrustorbopt." + ext) as lib_path:
    ff = CDLL(str(lib_path))


def solver_interface(
    func: Callable[[np.ndarray], float],
    update_orbs: Callable[
        [np.ndarray],
        Tuple[float, np.ndarray, np.ndarray, Callable[[np.ndarray], np.ndarray]],
    ],
    n_param: int,
    stability: Optional[bool] = None,
    line_search: Optional[bool] = None,
    conv_tol: Optional[float] = None,
    n_random_trial_vectors: Optional[int] = None,
    start_trust_radius: Optional[float] = None,
    n_macro: Optional[int] = None,
    n_micro: Optional[int] = None,
    global_red_factor: Optional[float] = None,
    local_red_factor: Optional[float] = None,
    verbose: Optional[int] = None,
    seed: Optional[int] = None,
):
    # callback function ctypes specifications, ctypes can only deal with simple return
    # types so we interface to Fortran subroutines by creating data to the relevant data
    hess_x_interface_type = CFUNCTYPE(
        None, POINTER(POINTER(c_double)), POINTER(POINTER(c_double))
    )
    update_orbs_interface_type = CFUNCTYPE(
        None,
        POINTER(POINTER(c_double)),
        POINTER(c_double),
        POINTER(POINTER(c_double)),
        POINTER(POINTER(c_double)),
        POINTER(hess_x_interface_type),
    )
    obj_func_interface_type = CFUNCTYPE(c_double, POINTER(POINTER(c_double)))

    @update_orbs_interface_type
    def update_orbs_interface(kappa_ptr, func_ptr, grad_ptr, h_diag_ptr, hess_x_ptr):
        """
        this function updates the orbitals and writes pointers to the function value,
        gradient, Hessian diagonal and returns a procedure function to the Hessian
        linear transformation
        """
        # function needs to be defined as a global variable to ensure that it is not
        # garbage collected when the current function completes
        global hess_x_interface

        # convert orbital rotation matrix pointer to numpy array
        kappa = np.ctypeslib.as_array(kappa_ptr[0], shape=(n_param,))

        # update orbitals and retrieve objective function, gradient, Hessian diagonal
        # and Hessian linear transformation function
        func_ptr[0], grad, h_diag, hess_x = update_orbs(kappa)

        # convert numpy arrays to pointers
        grad_ptr[0] = grad.ctypes.data_as(POINTER(c_double))
        h_diag_ptr[0] = h_diag.ctypes.data_as(POINTER(c_double))

        @hess_x_interface_type
        def hess_x_interface(x_ptr, hx_ptr):
            # convert trial vector pointer to numpy array
            x = np.ctypeslib.as_array(x_ptr[0], shape=(n_param,))

            # perform linear transformation and convert numpy array to pointer
            hx_ptr[0] = hess_x(x).ctypes.data_as(POINTER(c_double))

        # store the function pointer in hess_x_ptr
        hess_x_ptr[0] = hess_x_interface

    @obj_func_interface_type
    def obj_func_interface(kappa_ptr):
        """
        this function returns the function value
        """
        # convert orbital rotation matrix pointer to numpy array
        kappa = np.ctypeslib.as_array(kappa_ptr[0], shape=(n_param,))

        return func(kappa)

    # define result and argument types
    ff.solver.restype = None
    ff.solver.argtypes = [
        POINTER(update_orbs_interface_type),
        POINTER(obj_func_interface_type),
        c_long,
        POINTER(c_bool),
        POINTER(c_bool),
        POINTER(c_double),
        POINTER(c_long),
        POINTER(c_double),
        POINTER(c_long),
        POINTER(c_long),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_long),
        POINTER(c_long),
    ]

    # call Fortran function
    ff.solver(
        update_orbs_interface,
        obj_func_interface,
        n_param,
        None if stability is None else byref(c_bool(stability)),
        None if line_search is None else byref(c_bool(line_search)),
        None if conv_tol is None else byref(c_double(conv_tol)),
        (
            None
            if n_random_trial_vectors is None
            else byref(c_long(n_random_trial_vectors))
        ),
        None if start_trust_radius is None else byref(c_double(start_trust_radius)),
        None if n_macro is None else byref(c_long(n_macro)),
        None if n_micro is None else byref(c_long(n_micro)),
        None if global_red_factor is None else byref(c_double(global_red_factor)),
        None if local_red_factor is None else byref(c_double(local_red_factor)),
        None if verbose is None else byref(c_long(verbose)),
        None if seed is None else byref(c_long(seed)),
    )


def stability_check_interface(
    grad: np.ndarray,
    h_diag: np.ndarray,
    hess_x: Callable[[np.ndarray], np.ndarray],
    n_param: int,
    conv_tol: Optional[float] = None,
    n_random_trial_vectors: Optional[int] = None,
    n_iter: Optional[int] = None,
) -> Tuple[bool, np.ndarray]:
    # callback function ctypes specifications, ctypes can only deal with simple return
    # types so we interface to Fortran subroutines by creating data to the relevant data
    hess_x_interface_type = CFUNCTYPE(
        None, POINTER(POINTER(c_double)), POINTER(POINTER(c_double))
    )

    @hess_x_interface_type
    def hess_x_interface(x_ptr, hx_ptr):
        # convert trial vector pointer to numpy array
        x = np.ctypeslib.as_array(x_ptr[0], shape=(n_param,))

        # perform linear transformation and convert numpy array to pointer
        hx_ptr[0] = hess_x(x).ctypes.data_as(POINTER(c_double))

    # define result and argument types
    ff.stability_check.restype = None
    ff.stability_check.argtypes = [
        POINTER(c_double),
        POINTER(c_double),
        POINTER(hess_x_interface_type),
        c_long,
        POINTER(c_bool),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_long),
        POINTER(c_long),
    ]

    # initialize return variables
    stable = c_bool(False)
    kappa = np.empty(n_param, dtype=np.float64)

    # call Fortran function
    ff.stability_check(
        grad.ctypes.data_as(POINTER(c_double)),
        h_diag.ctypes.data_as(POINTER(c_double)),
        hess_x_interface,
        n_param,
        byref(stable),
        kappa.ctypes.data_as(POINTER(c_double)),
        None if conv_tol is None else byref(c_double(conv_tol)),
        (
            None
            if n_random_trial_vectors is None
            else byref(c_long(n_random_trial_vectors))
        ),
        None if n_iter is None else byref(c_long(n_iter)),
    )

    return bool(stable), kappa
