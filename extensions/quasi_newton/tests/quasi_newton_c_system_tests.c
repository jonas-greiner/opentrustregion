// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Pure-C system test for the public C interface declared in
// opentrustregion_quasi_newton.h.
//
// The Fortran-side quasi_newton_c_interface_unit_tests cover the bind(C) wrappers but
// never compile against the C header itself. This test does, so any drift between
// qn_settings_type_c (Fortran) and qn_settings_type (C) is caught here.
//

#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "opentrustregion.h"
#include "opentrustregion_quasi_newton.h"

// ---------------------------------------------------------------------------
// Compile-time layout checks for the C struct.
//
// These only verify that the C header is self-consistent: each field sits where the
// field order claims it does, with no surprise padding before the pointer block.
// Cross-language drift (Fortran vs. C) is caught at runtime by the default-value test
// below.
// ---------------------------------------------------------------------------

_Static_assert(offsetof(qn_settings_type, logger) == 0,
               "qn_settings_type: logger must be the first field");
_Static_assert(offsetof(qn_settings_type, initialized) == 1 * sizeof(void *),
               "qn_settings_type: initialized must follow logger");

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

bool test_qn_settings_init(void) {
  bool is_default_qn_settings(const qn_settings_type *settings);
  bool ok = true;

  // call function and compare with the default settings without logging function and
  // host context
  qn_settings_type s = qn_settings_init();
  if (!is_default_qn_settings(&s)) {
    fprintf(stderr, "test_qn_settings_init failed: Settings not initialized to the "
                    "default values without logging function and host context.\n");
    ok = false;
  }

  return ok;
}
