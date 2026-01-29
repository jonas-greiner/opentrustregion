# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import annotations

import numpy as np
from ctypes import CFUNCTYPE, POINTER, c_bool, c_void_p, Structure
from typing import TYPE_CHECKING
from pyopentrustregion.python_interface import (
    lib,
    c_int,
    c_real,
    obj_func_interface_type,
    update_orbs_interface_type,
    hess_x_interface_type,
    precond_interface_type,
    logger_interface_factory,
    Settings,
    auto_bind_fields,
)

if TYPE_CHECKING:
    from typing import Tuple, Callable, Optional, Any, Union


# callback function ctypes specifications, ctypes can only deal with simple return
# types so we interface to Fortran subroutines by creating pointers to the relevant
# data
get_energy_interface_type = CFUNCTYPE(c_int, POINTER(c_real), POINTER(c_real))
get_fock_interface_type = CFUNCTYPE(
    c_int, POINTER(c_real), POINTER(c_real), POINTER(c_real)
)
get_fock_jk_interface_type = CFUNCTYPE(
    c_int,
    POINTER(c_real),
    POINTER(c_real),
    POINTER(c_real),
    POINTER(c_real),
    POINTER(c_real),
)


# define classes corresponding to C structs for settings
class ARHSettingsC(Structure):
    _fields_ = [
        ("logger", c_void_p),
        ("initialized", c_bool),
        ("restricted", c_bool),
        ("verbose", c_int),
    ]


class ARHSettings(Settings):

    c_struct = ARHSettingsC
    try:
        init_c_struct = lib.init_arh_settings
    except AttributeError:
        raise AttributeError(
            "Please reinstall the package with: "
            "CMAKE_FLAGS='-DENABLE_ARH=ON' pip install ."
        )

    logger: Optional[Callable[[str], None]]
    logger_interface: Any


# ensure that appropriate fields are automatically set in settings_c object
auto_bind_fields(ARHSettings)


def arh_factory(
    dm_ao: np.ndarray,
    ao_overlap: np.ndarray,
    n_particle: int,
    n_ao: int,
    get_energy: Callable[[np.ndarray], float],
    get_fock: Union[
        Callable[[np.ndarray, np.ndarray], float],
        Callable[[np.ndarray, np.ndarray, np.ndarray, np.ndarray], float],
    ],
    settings: ARHSettings,
) -> Tuple[
    Callable[[np.ndarray], float],
    Callable[
        [np.ndarray, np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ],
    Callable[[np.ndarray, np.ndarray], None],
]:
    # get pointers to arrays
    dm_ao_ptr = dm_ao.ctypes.data_as(POINTER(c_real))
    ao_overlap_ptr = ao_overlap.ctypes.data_as(POINTER(c_real))

    # determine if closed-shell or open-shell formalism is used
    closed_shell = dm_ao.ndim == 2

    # define interfaces for callback functions
    @get_energy_interface_type
    def get_energy_interface(dm_ao_ptr, energy_ptr):
        """
        this function provides the interface to get the energy for a given density
        matrix
        """
        # convert matrix pointers to numpy arrays
        dm_ao = np.ctypeslib.as_array(
            dm_ao_ptr,
            shape=2 * (n_ao,) if closed_shell else (n_particle, n_ao, n_ao),
        )

        # get energy
        try:
            energy_ptr[0] = get_energy_interface.get_energy(dm_ao)
        except RuntimeError:
            return 1

        return 0

    @get_fock_interface_type
    def get_fock_interface(dm_ao_ptr, energy_ptr, fock_ptr):
        """
        this function provides the interface to get the energy and write the Fock
        matrix to the memory provided for a given density matrix
        """
        # convert matrix pointers to numpy arrays
        dm_ao = np.ctypeslib.as_array(dm_ao_ptr, shape=2 * (n_ao,))
        fock = np.ctypeslib.as_array(fock_ptr, shape=2 * (n_ao,))

        # get energy and Fock matrix
        try:
            energy_ptr[0] = get_fock_interface.get_fock(dm_ao, fock)
        except RuntimeError:
            return 1

        return 0

    @get_fock_jk_interface_type
    def get_fock_jk_interface(
        dm_ao_ptr, energy_ptr, fock_ptr, coulomb_ptr, exchange_ptr
    ):
        """
        this function provides the interface to get the energy and write the Fock,
        Coulomb, and Exchange matrices to the memory provided for a given density matrix
        """
        # convert matrix pointers to numpy arrays
        dm_ao = np.ctypeslib.as_array(dm_ao_ptr, shape=(n_particle, n_ao, n_ao))
        fock = np.ctypeslib.as_array(fock_ptr, shape=(n_particle, n_ao, n_ao))
        coulomb = np.ctypeslib.as_array(coulomb_ptr, shape=(n_particle, n_ao, n_ao))
        exchange = np.ctypeslib.as_array(exchange_ptr, shape=(n_particle, n_ao, n_ao))

        # get energy and Fock matrix
        try:
            energy_ptr[0] = get_fock_jk_interface.get_fock(
                dm_ao, fock, coulomb, exchange
            )
        except RuntimeError:
            return 1

        return 0

    # attach the energy and Fock matrix functions to the returned interface so that it
    # persists in Python to ensure that it is not garbage collected when the factory
    # completes
    get_energy_interface.get_energy = get_energy
    if closed_shell:
        get_fock_interface.get_fock = get_fock
    else:
        get_fock_jk_interface.get_fock = get_fock

    # set interfaces for optional callback functions, these need to be set here since
    # the interface might need parameters that are not known when the attribute to
    # settings is set (e.g. n_param)
    settings.set_optional_callback("logger", logger_interface_factory(settings.logger))

    if not hasattr(lib, "arh_factory"):
        raise RuntimeError(
            "Please reinstall the package with: "
            "CMAKE_FLAGS='-DENABLE_ARH=ON' pip install ."
        )

    # define result and argument types
    lib.arh_factory.restype = c_int
    lib.arh_factory.argtypes = [
        POINTER(c_real),
        POINTER(c_real),
        c_int,
        c_int,
        get_energy_interface_type,
        get_fock_interface_type if closed_shell else get_fock_jk_interface_type,
        POINTER(obj_func_interface_type),
        POINTER(update_orbs_interface_type),
        POINTER(precond_interface_type),
        ARHSettingsC,
    ]

    # call Fortran function
    obj_func_arh_funptr = obj_func_interface_type()
    update_orbs_arh_funptr = update_orbs_interface_type()
    precond_arh_funptr = precond_interface_type()
    error = lib.arh_factory(
        dm_ao_ptr,
        ao_overlap_ptr,
        n_particle,
        n_ao,
        get_energy_interface,
        get_fock_interface if closed_shell else get_fock_jk_interface,
        obj_func_arh_funptr,
        update_orbs_arh_funptr,
        precond_arh_funptr,
        settings.settings_c,
    )

    if error:
        raise RuntimeError("OpenTrustRegion ARH factory produced error.")

    # python wrappers which converts numpy arrays to ctypes and calls the C functions
    def obj_func_arh_interface(kappa: np.ndarray) -> float:
        # initialize real
        func = c_real()

        # get pointers to arrays
        kappa_ptr = kappa.ctypes.data_as(POINTER(c_real))

        # update orbital function
        error = obj_func_arh_funptr(kappa_ptr, func)
        if error != 0:
            raise RuntimeError("ARH orbital updating function raised error.")

        return func.value

    def update_orbs_arh_interface(
        kappa: np.ndarray, grad: np.ndarray, h_diag: np.ndarray
    ) -> Tuple[float, Callable[[np.ndarray, np.ndarray], None]]:
        # initialize real
        func = c_real()

        # get pointers to arrays
        kappa_ptr = kappa.ctypes.data_as(POINTER(c_real))
        grad_ptr = grad.ctypes.data_as(POINTER(c_real))
        h_diag_ptr = h_diag.ctypes.data_as(POINTER(c_real))

        # initialize Hessian linear transformation function pointer
        hess_x_arh_funptr = hess_x_interface_type()

        # update orbital function
        error = update_orbs_arh_funptr(
            kappa_ptr, func, grad_ptr, h_diag_ptr, hess_x_arh_funptr
        )
        if error != 0:
            raise RuntimeError("ARH orbital updating function raised error.")

        def hess_x_arh_interface(x: np.ndarray, hess_x: np.ndarray):
            # get pointers to arrays
            x_ptr = x.ctypes.data_as(POINTER(c_real))
            hess_x_ptr = hess_x.ctypes.data_as(POINTER(c_real))

            # update orbital function
            error = hess_x_arh_funptr(x_ptr, hess_x_ptr)
            if error != 0:
                raise RuntimeError(
                    "ARH Hessian linear transformation function raised error."
                )

            return hess_x

        return func.value, hess_x_arh_interface

    def precond_arh_interface(
        residual: np.ndarray, mu: float, precond_residual: np.ndarray
    ):
        # get pointers to arrays
        residual_ptr = residual.ctypes.data_as(POINTER(c_real))
        precond_residual_ptr = precond_residual.ctypes.data_as(POINTER(c_real))

        # update orbital function
        error = precond_arh_funptr(residual_ptr, c_real(mu), precond_residual_ptr)
        if error != 0:
            raise RuntimeError("ARH preconditioning function raised error.")

        return

    # keep all functions alive that are involved in the ARH orbital update
    obj_func_arh_interface.get_energy_interface = get_energy_interface
    obj_func_arh_interface.obj_func_arh_funptr = obj_func_arh_funptr
    update_orbs_arh_interface.get_energy_interface = get_energy_interface
    update_orbs_arh_interface.get_fock_interface = (
        get_fock_interface if closed_shell else get_fock_jk_interface
    )
    update_orbs_arh_interface.update_orbs_arh_funptr = update_orbs_arh_funptr
    precond_arh_interface.precond_arh_funptr = precond_arh_funptr

    return obj_func_arh_interface, update_orbs_arh_interface, precond_arh_interface


def arh_deconstructor(dm_ao: np.ndarray):
    if not hasattr(lib, "arh_deconstructor"):
        raise RuntimeError(
            "Please reinstall the package with: "
            "CMAKE_FLAGS='-DENABLE_ARH=ON' pip install ."
        )

    # define result and argument types
    lib.arh_deconstructor.restype = c_int
    lib.arh_deconstructor.argtypes = [POINTER((c_real))]

    # call Fortran function
    error = lib.arh_deconstructor(dm_ao.ctypes.data_as(POINTER(c_real)))

    if error:
        raise RuntimeError("OpenTrustRegion ARH deconstructor produced error.")

    return
