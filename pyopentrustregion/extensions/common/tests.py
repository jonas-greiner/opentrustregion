# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import unittest
from ctypes import c_void_p, POINTER

from pyopentrustregion.tests import lib, NUMPY_AVAILABLE, add_tests, print_separator
from pyopentrustregion.python_interface import c_real, c_int

if NUMPY_AVAILABLE:
    import numpy as np
    from pyopentrustregion.extensions.common.python_interface import (
        EvaluateDMInterface,
    )


# define all tests in alphabetical order
fortran_tests = {
    "common_tests": [
        "channel_rows",
        "compute_sqrt_and_inv_sqrt",
        "fill_extra_trial_vectors_orbital_basis",
        "init_orbital_settings",
        "level_shifted_divisors",
        "matrix_exponential",
        "positive_definite_divisors",
        "precond_orbital_basis",
        "precond_pd_orbital_basis",
        "refresh_response_orbital_basis",
    ],
}

# the shared extension routines are only built together with an extension which uses
# them, so raise an AttributeError on import otherwise
getattr(lib, "test_" + fortran_tests["common_tests"][0])

# number of AOs
n_ao = c_int.in_dll(lib, "test_n_ao").value


@add_tests
class CommonTests(unittest.TestCase):
    """
    this class contains unit tests for the routines shared between extensions
    """

    tests = fortran_tests["common_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for shared extension routines...")
        return super().setUpClass()


@unittest.skipUnless(NUMPY_AVAILABLE, "NumPy not available.")
class CommonPyInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the Python interface shared between extensions
    """

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for shared extension Python interface...")
        return super().setUpClass()

    def test_evaluate_dm_py_interface(self):
        """
        this function tests the Python interface to the density matrix evaluating
        function
        """
        # initialize test flag
        test_passed = True

        # check if a density matrix evaluating function which returns no response
        # function although one is requested produces an error
        def evaluate_dm_without_response(dm_ao, fock, get_response):
            return np.sum(dm_ao), None

        exception = {}
        evaluate_dm = EvaluateDMInterface(
            evaluate_dm_without_response, n_ao, 1, True, exception
        )
        dm_ao = np.full(2 * (n_ao,), 1.0, dtype=np.float64)
        energy = (c_real * 1)()
        get_response_funptr = (c_void_p * 1)()
        error = evaluate_dm(
            dm_ao.ctypes.data_as(POINTER(c_real)), energy, None, get_response_funptr
        )
        if error == 0 or not isinstance(exception.get("exc"), RuntimeError):
            print(
                " test_evaluate_dm_py_interface failed: Missing response function "
                "requested from density matrix evaluating function does not produce "
                "an error."
            )
            test_passed = False

        self.assertTrue(test_passed, "test_evaluate_dm_py_interface failed")
        print(" test_evaluate_dm_py_interface PASSED")
