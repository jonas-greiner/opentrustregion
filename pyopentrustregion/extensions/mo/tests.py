# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import unittest

from pyopentrustregion.tests import lib, add_tests, print_separator


# define all tests in alphabetical order
fortran_tests = {
    "mo_tests": [
        "calculate_grad_h_diag_mo",
        "finalize_mo",
        "get_extra_trial_vectors_mo",
        "get_extra_trial_vectors_mo_callback",
        "get_hess_eigval_pairs_mo",
        "hess_x_mo_callback",
        "hess_x_static_mo",
        "init_mo_settings",
        "mo_deconstructor",
        "mo_factory_common",
        "mo_factory_cs",
        "mo_factory_os",
        "mo_sanity_check",
        "mo_set_solver_settings",
        "mo_transform",
        "obj_func_mo_callback",
        "precond_mo_callback",
        "precond_pd_mo_callback",
        "refresh_hess_eigen_mo",
        "rotate_from_hess_eigenbasis_mo",
        "rotate_mo_coeff",
        "rotate_orbitals_mo",
        "rotate_to_hess_eigenbasis_mo",
        "update_orbs_mo_callback",
    ],
}

# the MO routines are only built together with the MO extension, so raise an
# AttributeError on import otherwise
getattr(lib, "test_" + fortran_tests["mo_tests"][0])


@add_tests
class MOTests(unittest.TestCase):
    """
    this class contains unit tests for the MO extension
    """

    tests = fortran_tests["mo_tests"]

    @classmethod
    def setUpClass(cls):
        print_separator("Running unit tests for MO...")
        return super().setUpClass()
