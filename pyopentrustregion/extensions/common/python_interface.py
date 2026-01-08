# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import annotations

import numpy as np
from ctypes import CFUNCTYPE, POINTER
from typing import TYPE_CHECKING
from pyopentrustregion.python_interface import c_int, c_real

if TYPE_CHECKING:
    from typing import Callable, Any


# callback function ctypes specifications, ctypes can only deal with simple return
# types so we interface to Fortran subroutines by creating pointers to the relevant
# data
change_reference_interface_type = CFUNCTYPE(
    c_int, POINTER(c_real), c_int, POINTER(c_real), POINTER(c_real), POINTER(c_real)
)


# define interface factories
def change_reference_interface_factory(
    change_reference: Callable[[np.ndarray, np.ndarray, np.ndarray, np.ndarray], None],
    n_param: int,
) -> Any:
    """
    this function is a factory for the change reference interface
    """

    @change_reference_interface_type
    def change_reference_interface(
        new_ref_ptr, n_points, kappa_list_ptr, local_grad_list_ptr, grad_list_ptr
    ):
        """
        this function provides the interface to change the reference and to write the
        parameter, local gradient, and gradient lists to the memory provided through
        pointers
        """
        # convert pointers to numpy arrays and deal with row- vs column-major order
        new_ref = np.ctypeslib.as_array(new_ref_ptr, shape=(n_param,))
        kappa_list = np.ctypeslib.as_array(kappa_list_ptr, shape=(n_points, n_param))
        local_grad_list = np.ctypeslib.as_array(
            local_grad_list_ptr, shape=(n_points, n_param)
        )
        grad_list = np.ctypeslib.as_array(grad_list_ptr, shape=(n_points, n_param))

        # change reference and retrieve parameter, local gradient, and gradient lists
        try:
            change_reference(new_ref, kappa_list, local_grad_list, grad_list)
        except RuntimeError:
            return 1

        return 0

    return change_reference_interface
