# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import annotations

import sys
import traceback
import numpy as np
from ctypes import CFUNCTYPE, POINTER, byref, c_bool, c_void_p, Structure
from dataclasses import dataclass
from typing import TYPE_CHECKING
from pyopentrustregion.python_interface import (
    lib,
    c_int,
    c_real,
    obj_func_interface_type,
    update_orbs_interface_type,
    project_interface_type,
    logger_interface_type,
    LoggerInterface,
    Settings,
    auto_bind_fields,
)
from pyopentrustregion.extensions.common.python_interface import UpdateOrbsPyInterface

if TYPE_CHECKING:
    from typing import Tuple, Callable, Optional, Any


# callback function ctypes specifications, ctypes can only deal with simple return
# types so we interface to Fortran subroutines by creating pointers to the relevant
# data
get_energy_interface_type = CFUNCTYPE(c_int, POINTER(c_real), POINTER(c_real))
get_response_interface_type = CFUNCTYPE(c_int, POINTER(c_real), POINTER(c_real))
update_dm_interface_type = CFUNCTYPE(
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
        ("restricted", c_bool),
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
class GetEnergyInterface:
    """
    this class provides the interface to get the energy for a given density matrix
    """

    get_energy: Callable[[np.ndarray], float]
    n_ao: int
    n_particle: int
    closed_shell: bool

    def __call__(self, dm_ao_ptr, energy_ptr) -> int:
        # convert matrix pointers to numpy arrays
        dm_ao = np.ctypeslib.as_array(
            dm_ao_ptr,
            shape=(
                2 * (self.n_ao,)
                if self.closed_shell
                else (self.n_particle, self.n_ao, self.n_ao)
            ),
        )

        # get energy
        try:
            energy_ptr[0] = self.get_energy(dm_ao)
        except RuntimeError:
            traceback.print_exc(file=sys.stderr)
            return 1

        return 0


@dataclass
class GetResponseInterface:
    """
    this class provides the interface to write the response to the memory provided for
    a given density matrix
    """

    get_response: Callable[[np.ndarray, np.ndarray], None]
    n_ao: int
    n_particle: int
    closed_shell: bool

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
        except RuntimeError:
            traceback.print_exc(file=sys.stderr)
            return 1

        return 0


@dataclass
class UpdateDMInterface:
    """
    this class provides the interface to the density matrix updating function
    """

    update_dm: Callable[
        [np.ndarray, np.ndarray], Tuple[float, Callable[[np.ndarray, np.ndarray], None]]
    ]
    n_ao: int
    n_particle: int
    closed_shell: bool
    get_response_funptr: Optional[Any] = None

    def __call__(self, dm_ao_ptr, energy_ptr, fock_ptr, get_response_funptr) -> int:
        # convert matrix pointers to numpy arrays
        if self.closed_shell:
            dm_ao = np.ctypeslib.as_array(dm_ao_ptr, shape=2 * (self.n_ao,))
            fock = np.ctypeslib.as_array(fock_ptr, shape=2 * (self.n_ao,))
        else:
            dm_ao = np.ctypeslib.as_array(
                dm_ao_ptr, shape=(self.n_particle, self.n_ao, self.n_ao)
            )
            fock = np.ctypeslib.as_array(
                fock_ptr, shape=(self.n_particle, self.n_ao, self.n_ao)
            )

        # get energy, Fock matrix, and response function
        try:
            energy_ptr[0], get_response = self.update_dm(dm_ao, fock)
        except RuntimeError:
            traceback.print_exc(file=sys.stderr)
            return 1

        # attach the response interface to the object so that it persists in Python
        # to ensure that it is not garbage collected when the factory completes
        self.get_response_funptr = get_response_interface_type(
            GetResponseInterface(
                get_response, self.n_ao, self.n_particle, self.closed_shell
            )
        )
        get_response_funptr[0] = self.get_response_funptr

        return 0


@dataclass
class ObjFuncPyInterface:
    """
    this class provides the Python interface to the objective function,
    get_energy_interface is stored to ensure that it is not garbage collected when the
    factory completes
    """

    obj_func_funptr: Any
    get_energy_interface: Any

    def __call__(self, kappa: np.ndarray) -> float:
        # initialize real
        func = c_real()

        # get pointers to arrays
        kappa_ptr = kappa.ctypes.data_as(POINTER(c_real))

        # update orbital function
        error = self.obj_func_funptr(kappa_ptr, byref(func))
        if error != 0:
            raise RuntimeError("Objective function raised error.")

        return func.value


@dataclass
class ProjectPyInterface:
    """
    this class provides the Python interface to the projection function
    """

    project_funptr: Any

    def __call__(self, vector: np.ndarray):
        # get pointers to arrays
        vector_ptr = vector.ctypes.data_as(POINTER(c_real))

        # update orbital function
        error = self.project_funptr(vector_ptr)
        if error != 0:
            raise RuntimeError("Projection function raised error.")

        return


def oao_factory(
    dm_ao: np.ndarray,
    ao_overlap: np.ndarray,
    n_particle: int,
    n_ao: int,
    get_energy: Callable[[np.ndarray], float],
    update_dm: Callable[
        [np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ],
    settings: OAOSettings,
) -> Tuple[
    Callable[[np.ndarray], float],
    Callable[
        [np.ndarray, np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ],
    Callable[[np.ndarray], None],
]:
    # get pointers to arrays
    dm_ao_ptr = dm_ao.ctypes.data_as(POINTER(c_real))
    ao_overlap_ptr = ao_overlap.ctypes.data_as(POINTER(c_real))

    # determine if closed-shell or open-shell formalism is used
    closed_shell = dm_ao.ndim == 2

    # define interfaces for callback functions
    get_energy_interface = get_energy_interface_type(
        GetEnergyInterface(
            get_energy=get_energy,
            n_ao=n_ao,
            n_particle=n_particle,
            closed_shell=closed_shell,
        )
    )
    update_dm_interface = update_dm_interface_type(
        UpdateDMInterface(update_dm, n_ao, n_particle, closed_shell)
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
        get_energy_interface_type,
        update_dm_interface_type,
        POINTER(obj_func_interface_type),
        POINTER(update_orbs_interface_type),
        POINTER(project_interface_type),
        POINTER(OAOSettingsC),
    ]

    # call Fortran function
    obj_func_oao_funptr = obj_func_interface_type()
    update_orbs_oao_funptr = update_orbs_interface_type()
    project_oao_funptr = project_interface_type()
    error = lib.oao_factory(
        dm_ao_ptr,
        ao_overlap_ptr,
        n_particle,
        n_ao,
        get_energy_interface,
        update_dm_interface,
        byref(obj_func_oao_funptr),
        byref(update_orbs_oao_funptr),
        byref(project_oao_funptr),
        byref(settings.settings_c),
    )

    if error:
        raise RuntimeError("OpenTrustRegion OAO factory produced error.")

    return (
        ObjFuncPyInterface(
            obj_func_funptr=obj_func_oao_funptr,
            get_energy_interface=get_energy_interface,
        ),
        UpdateOrbsPyInterface(
            update_orbs_funptr=update_orbs_oao_funptr,
            saved_objects={
                "get_energy_interface": get_energy_interface,
                "update_dm_interface": update_dm_interface,
            },
        ),
        ProjectPyInterface(project_funptr=project_oao_funptr),
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
