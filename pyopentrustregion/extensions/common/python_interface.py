# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import annotations

import numpy as np
from ctypes import POINTER, byref
from dataclasses import dataclass
from typing import TYPE_CHECKING
from pyopentrustregion.python_interface import c_real, hess_x_interface_type

if TYPE_CHECKING:
    from typing import Tuple, Callable, Optional, Any, Dict


# define interface factories
@dataclass
class HessXPyInterface:
    """
    this class provides the Python interface to the Hessian linear transformation
    function
    """

    hess_x_funptr: Any

    def __call__(self, x: np.ndarray, hess_x: np.ndarray):
        # get pointers to arrays
        x_ptr = x.ctypes.data_as(POINTER(c_real))
        hess_x_ptr = hess_x.ctypes.data_as(POINTER(c_real))

        # update orbital function
        error = self.hess_x_funptr(x_ptr, hess_x_ptr)
        if error != 0:
            raise RuntimeError("Hessian linear transformation function raised error.")


@dataclass
class UpdateOrbsPyInterface:
    """
    this class provides the Python interface to the orbital updating function,
    get_energy_interface and get_fock_interface are stored to ensure that they are not
    garbage collected when the factory completes
    """

    update_orbs_funptr: Any
    saved_objects: Dict[str, Any]
    hess_x: Optional[Callable[[np.ndarray, np.ndarray], None]] = None

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
            raise RuntimeError("OAO orbital updating function raised error.")

        # attach the Hessian linear transformation interface to the object so that it
        # persists in Python to ensure that it is not garbage collected when the
        # factory completes
        self.hess_x = HessXPyInterface(hess_x_funptr)

        return func.value, self.hess_x
