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
typedef c_int get_response_fn(const c_real *dm_ao_c, c_real *response_c);
typedef get_response_fn *get_response_fp;

/* Density matrix evaluating callback */
typedef c_int evaluate_dm_fn(const c_real *dm_ao_c, c_real *energy_c, c_real *fock_c,
                             get_response_fp *get_response_ptr);
typedef evaluate_dm_fn *evaluate_dm_fp;

#ifdef __cplusplus
}
#endif

#endif
