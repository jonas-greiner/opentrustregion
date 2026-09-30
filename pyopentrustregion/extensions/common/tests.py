# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import unittest

from pyopentrustregion.tests import lib, add_tests, print_separator


# define all tests in alphabetical order
fortran_tests = {
    "common_tests": [
        "compute_sqrt_and_inv_sqrt",
        "level_shifted_divisors",
        "matrix_exponential",
        "positive_definite_divisors",
    ],
}

# the shared extension routines are only built together with an extension which uses
# them, so raise an AttributeError on import otherwise
getattr(lib, "test_" + fortran_tests["common_tests"][0])


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
