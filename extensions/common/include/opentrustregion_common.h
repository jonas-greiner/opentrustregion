// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENTRUSTREGION_COMMON_H
#define OPENTRUSTREGION_COMMON_H

#include "opentrustregion.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------
 * Function pointer types shared by the extensions parameterizing
 * the orbitals in an orbital basis
 * ------------------------------------------------------------------ */

/* Response callback */
typedef c_int get_response_fn(const c_real *dm_ao_c, c_real *response_c, void *context);
typedef get_response_fn *get_response_fp;

/* Density matrix evaluating callback */
typedef c_int evaluate_dm_fn(const c_real *dm_ao_c, c_real *energy_c, c_real *fock_c,
                             get_response_fp *get_response_ptr, void *context);
typedef evaluate_dm_fn *evaluate_dm_fp;

/* ------------------------------------------------------------------
 * Functions shared by every extension
 * ------------------------------------------------------------------ */

/**
 * Frees the extension object a factory returned, together with everything the object
 * holds. Release every extension exactly once, never while a solve uses it, and only
 * once no extension wired after it into the same solver settings is used any more.
 *
 * @param extension                               Extension a factory returned
 */
void release_extension(void *extension);

#ifdef __cplusplus
}
#endif

#endif
