// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Pure-C system test for the public C interface declared in opentrustregion_arh.h.
//
// The Fortran-side arh_c_interface_unit_tests cover the bind(C) wrappers but never 
// compile against the C header itself. This test does, so any drift between
// arh_settings_type_c (Fortran) and arh_settings_type (C) is caught here.
//

#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "opentrustregion.h"
#include "opentrustregion_arh.h"

// ---------------------------------------------------------------------------
// Compile-time layout checks for the C struct.
//
// These only verify that the C header is self-consistent: each field sits where the 
// field order claims it does, with no surprise padding before the pointer block. 
// Cross-language drift (Fortran vs. C) is caught at runtime by the default-value test 
// below.
// ---------------------------------------------------------------------------

_Static_assert(offsetof(arh_settings_type, logger) == 0,
               "arh_settings_type: logger must be the first field");
_Static_assert(offsetof(arh_settings_type, initialized) == 1 * sizeof(void*),
               "arh_settings_type: initialized must follow logger");

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

bool test_arh_settings_init(void)
{
    // get defaults
    void get_default_arh_values(arh_settings_type *settings);
    arh_settings_type defaults = {0};
    get_default_arh_values(&defaults);

    // call function
    arh_settings_type s = arh_settings_init();

    // compare values
    bool ok = true;
    if (!s.initialized) {
        fprintf(stderr,
                "test_arh_settings_init failed: Settings not initialized.\n");
        ok = false;
    }
    if (s.verbose != defaults.verbose) {
        fprintf(stderr,
                "test_arh_settings_init failed: Verbosity parameter wrong.\n");
        ok = false;
    }
    if (strcmp(s.arh_type, defaults.arh_type) != 0) {
        fprintf(stderr,
                "test_arh_settings_init failed: ARH type parameter wrong.\n");
        ok = false;
    }
    if (s.logger) {
        fprintf(stderr,
                "test_arh_settings_init failed: Callback pointers should be "
                "NULL.\n");
        ok = false;
    }

    return ok;
}
