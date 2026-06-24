// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENTRUSTREGION_ARH_H
#define OPENTRUSTREGION_ARH_H

#include "opentrustregion.h"
#include "opentrustregion_oao.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------
 * Declarations for ARH functions and function pointer types
 * ------------------------------------------------------------------ */

/* Density matrix updating callback with same-spin and opposite-spin contributions */
typedef c_int update_dm_spin_fn(
    c_real* v_same_spin_c,
    c_real* v_opposite_spin_c,
typedef update_dm_spin_fn *update_dm_spin_fp;

/* ------------------------------------------------------------------
 * Struct corresponding to Fortran type(arh_settings_type_c)
 * ------------------------------------------------------------------ */
typedef struct {
    logger_fp logger;
    c_bool initialized;
    c_bool restricted;
    c_int verbose;
    char arh_type[OTR_KW_LEN + 1];
} arh_settings_type;

// Fortran-callable init routine for ARH settings
void init_arh_settings(arh_settings_type* settings);

/* ------------------------------------------------------------------
 * Fortran wrappers
 * ------------------------------------------------------------------ */

* @param update_dm_c_funptr         C pointer to update_dm or update_dm_spin callback
    update_dm_spin_fp update_dm_c_funptr,

/**
 * Fortran-callable ARH deconstructor.
 */
void arh_deconstructor();

#ifdef __cplusplus
}
#endif

/* ------------------------------------------------------------------
 * Small C helper functions to mimic Fortran settings%init()
 * ------------------------------------------------------------------ */

static inline arh_settings_type arh_settings_init(void) {
    arh_settings_type s = {0};
    init_arh_settings(&s);
    return s;
}

#endif
