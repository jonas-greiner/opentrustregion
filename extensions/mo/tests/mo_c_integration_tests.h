// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/* Code the C integration tests of MO share with those of ARH: the setup of a solve of
 * the toy model of common_c_integration_tests.h through the MO factory. */

#ifndef MO_C_INTEGRATION_TESTS_H
#define MO_C_INTEGRATION_TESTS_H

#include "../../common/tests/common_c_integration_tests.h"

/* set up a solve of the toy model in the MO basis and the given shell through the MO
 * factory, with MO settings whose host context is the host of the orbital basis, and
 * return whether the factory succeeded */
bool setup_mo_toy_run(toy_run *run, c_int n_particle, const char *test_name,
                      const char *case_name);

#endif
