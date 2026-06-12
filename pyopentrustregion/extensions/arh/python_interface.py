# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import annotations

import numpy as np
from ctypes import CFUNCTYPE, POINTER, c_bool, c_void_p, c_char, Structure, byref
from dataclasses import dataclass
from inspect import signature
from typing import TYPE_CHECKING
from pyopentrustregion.python_interface import (
    lib,
    c_int,
    c_real,
    kw_len,
    obj_func_interface_type,
    update_orbs_interface_type,
    project_interface_type,
    logger_interface_type,
    LoggerInterface,
    Settings,
    auto_bind_fields,
)
from pyopentrustregion.extensions.common.python_interface import UpdateOrbsPyInterface
from pyopentrustregion.extensions.oao.python_interface import (
    get_energy_interface_type,
    get_response_interface_type,
    update_dm_interface_type,
    GetEnergyInterface,
    UpdateDMInterface,
    GetResponseInterface,
    ObjFuncPyInterface,
    ProjectPyInterface,
)

if TYPE_CHECKING:
    from typing import Tuple, Callable, Optional, Any, Union, Callable, TypeGuard

    # Python type specifications
    UpdateDMType = Callable[
        [np.ndarray, np.ndarray], Tuple[float, Callable[[np.ndarray, np.ndarray], None]]
    ]
    UpdateDMJKType = Callable[
        [np.ndarray, np.ndarray, np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ]


# type guards
def is_update_dm(func: Any) -> TypeGuard[UpdateDMType]:
    try:
        sig = signature(func)
        return len(sig.parameters) == 2
    except (ValueError, TypeError):
        return False


def is_update_dm_jk(func: Any) -> TypeGuard[UpdateDMJKType]:
    try:
        sig = signature(func)
        return len(sig.parameters) == 4
    except (ValueError, TypeError):
        return False


# callback function ctypes specifications, ctypes can only deal with simple return
# types so we interface to Fortran subroutines by creating pointers to the relevant
# data
update_dm_jk_interface_type = CFUNCTYPE(
    c_int,
    POINTER(c_real),
    POINTER(c_real),
    POINTER(c_real),
    POINTER(c_real),
    POINTER(c_real),
    POINTER(get_response_interface_type),
)


# define classes corresponding to C structs for settings
class ARHSettingsC(Structure):
    _fields_ = [
        ("logger", c_void_p),
        ("initialized", c_bool),
        ("restricted", c_bool),
        ("verbose", c_int),
        ("arh_type", c_char * (kw_len + 1)),
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


# define interface factories
@dataclass
class UpdateDMJKInterface:
    """
    this class provides the interface density matrix updating function with Coulomb and
    exchange contributions
    """

    update_dm_jk: Callable[
        [np.ndarray, np.ndarray, np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ]
    n_ao: int
    n_particle: int
    closed_shell: bool
    get_response_funptr: Optional[Any] = None

    def __call__(
        self,
        dm_ao_ptr,
        energy_ptr,
        fock_ptr,
        coulomb_ptr,
        exchange_ptr,
        get_response_funptr,
    ) -> int:
        # convert matrix pointers to numpy arrays
        if self.closed_shell:
            dm_ao = np.ctypeslib.as_array(dm_ao_ptr, shape=2 * (self.n_ao,))
            fock = np.ctypeslib.as_array(fock_ptr, shape=2 * (self.n_ao,))
            coulomb = np.ctypeslib.as_array(coulomb_ptr, shape=2 * (self.n_ao,))
            exchange = np.ctypeslib.as_array(exchange_ptr, shape=2 * (self.n_ao,))
        else:
            dm_ao = np.ctypeslib.as_array(
                dm_ao_ptr, shape=(self.n_particle, self.n_ao, self.n_ao)
            )
            fock = np.ctypeslib.as_array(
                fock_ptr, shape=(self.n_particle, self.n_ao, self.n_ao)
            )
            coulomb = np.ctypeslib.as_array(
                coulomb_ptr, shape=(self.n_particle, self.n_ao, self.n_ao)
            )
            exchange = np.ctypeslib.as_array(
                exchange_ptr, shape=(self.n_particle, self.n_ao, self.n_ao)
            )

        # get energy, Fock matrix, and response function
        try:
            energy_ptr[0], get_response = self.update_dm_jk(
                dm_ao, fock, coulomb, exchange
            )
        except RuntimeError:
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


def arh_factory(
    dm_ao: np.ndarray,
    ao_overlap: np.ndarray,
    n_particle: int,
    n_ao: int,
    get_energy: Callable[[np.ndarray], float],
    update_dm: Union[
        Callable[
            [np.ndarray, np.ndarray],
            Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
        ],
        Callable[
            [np.ndarray, np.ndarray, np.ndarray, np.ndarray],
            Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
        ],
    ],
    settings: ARHSettings,
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
    if is_update_dm(update_dm):
        update_dm_interface = update_dm_interface_type(
            UpdateDMInterface(update_dm, n_ao, n_particle, closed_shell)
        )
    elif is_update_dm_jk(update_dm):
        update_dm_jk_interface = update_dm_jk_interface_type(
            UpdateDMJKInterface(update_dm, n_ao, n_particle, closed_shell)
        )

    # set interfaces for optional callback functions, these need to be set here since
    # the interface might need parameters that are not known when the attribute to
    # settings is set (e.g. n_param)
    settings.set_optional_callback(
        "logger", settings.logger, LoggerInterface, logger_interface_type
    )

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
        update_dm_interface_type if closed_shell else update_dm_jk_interface_type,
        POINTER(obj_func_interface_type),
        POINTER(update_orbs_interface_type),
        POINTER(project_interface_type),
        POINTER(ARHSettingsC),
    ]

    # call Fortran function
    obj_func_arh_funptr = obj_func_interface_type()
    update_orbs_arh_funptr = update_orbs_interface_type()
    project_arh_funptr = project_interface_type()
    error = lib.arh_factory(
        dm_ao_ptr,
        ao_overlap_ptr,
        n_particle,
        n_ao,
        get_energy_interface,
        update_dm_interface if closed_shell else update_dm_jk_interface,
        byref(obj_func_arh_funptr),
        byref(update_orbs_arh_funptr),
        byref(project_arh_funptr),
        byref(settings.settings_c),
    )

    if error:
        raise RuntimeError("OpenTrustRegion ARH factory produced error.")

    return (
        ObjFuncPyInterface(
            obj_func_funptr=obj_func_arh_funptr,
            get_energy_interface=get_energy_interface,
        ),
        UpdateOrbsPyInterface(
            update_orbs_funptr=update_orbs_arh_funptr,
            saved_objects={
                "get_energy_interface": get_energy_interface,
                "update_dm_interface": (
                    update_dm_interface if closed_shell else update_dm_jk_interface
                ),
            },
        ),
        ProjectPyInterface(project_funptr=project_arh_funptr),
    )
