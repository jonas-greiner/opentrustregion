// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/* Pure-C integration tests for the public C interface declared in
 * opentrustregion_quasi_newton.h
 *
 * The Fortran-side qn_c_interface_unit_tests cover the bind(C) wrappers but never
 * compile against the C header itself. These tests do, and test_qn_settings_layout
 * checks that the settings struct has the size of the bind(C) type and that every
 * field is read back under its own name, so any drift between the bind(C) type
 * (Fortran) and the struct (C) is caught here, except a field left in the struct that
 * fits into the padding of the bind(C) type, which changes neither the size nor the
 * offset of a named field. The header's inline initialization, which
 * test_qn_settings_init checks, takes both sides from the same Fortran assignment, so
 * it cannot detect such drift. */

#include <stdint.h>
#include <stdio.h>

#include "opentrustregion.h"
#include "opentrustregion_quasi_newton.h"

#include "../../common/tests/common_c_integration_tests.h"

/* ------------------------------------------------------------------
 * Reference settings provided by qn_test_reference.f90: the reference value of a field
 * by its name, C settings filled with these values field by field, the size of the
 * bind(C) settings type, and a comparison with the default settings
 * ------------------------------------------------------------------ */

void reference_qn_field(const char *name, c_real *value, char *keyword);
void get_reference_qn_values(qn_settings_type *settings);
size_t qn_settings_size(void);
bool is_default_qn_settings(const qn_settings_type *settings);

/* ------------------------------------------------------------------
 * Tests
 * ------------------------------------------------------------------ */

bool test_qn_settings_layout(void) {
  /* check that the settings struct has the size of the bind(C) settings type, so that
   * a field missing from either is detected, which the field checks below would not
   * read, and stop otherwise, since filling smaller settings would write past their
   * end */
  if (sizeof(qn_settings_type) != qn_settings_size()) {
    fprintf(stderr,
            "test_qn_settings_layout failed: Size of quasi-Newton settings differs "
            "from the library's.\n");
    return false;
  }

  /* check that every field is read back under its own name */
  bool ok = true;
  qn_settings_type s = {0};
  get_reference_qn_values(&s);
  reference_field_fn *lookup = reference_qn_field;
  const char *test_name = "test_qn_settings_layout";
  ok &= check_settings_field(lookup, (uintptr_t)s.logger, "logger", test_name);
  ok &= check_settings_field(lookup, s.initialized, "initialized", test_name);
  ok &= check_settings_field(lookup, s.verbose, "verbose", test_name);
  ok &= check_settings_field(lookup, s.max_points, "max_points", test_name);
  ok &= check_settings_keyword(lookup, s.hess_update_scheme, "hess_update_scheme",
                               test_name);
  ok &= check_settings_field(lookup, (uintptr_t)s.context, "context", test_name);

  return ok;
}

bool test_qn_settings_init(void) {
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
