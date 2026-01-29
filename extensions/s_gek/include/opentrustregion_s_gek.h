// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENTRUSTREGION_S_GEK_H
#define OPENTRUSTREGION_S_GEK_H

#include "opentrustregion.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------
 * Struct corresponding to Fortran type(s_gek_settings_type_c)
 * ------------------------------------------------------------------ */
typedef struct {
    logger_fp logger;
    c_bool initialized;
    c_bool use_subspace;
    c_int verbose;
    c_int max_points;
} s_gek_settings_type;

// Fortran-callable init routine for S-GEK settings
void init_s_gek_settings(s_gek_settings_type* settings);

/* ------------------------------------------------------------------
 * Fortran wrappers
 * ------------------------------------------------------------------ */

/**
 * Fortran-callable update_orbs_s_gek_factory interface
 *
 * @param update_orbs_orig_c_funptr    Pointer to original update_orbs callback
 * @param n_param                      Number of parameters
 * @param settings                     S-GEK settings
 * @param update_orbs_s_gek_c_funptr   Pointer to new update_orbs callback
 * @return                             Integer error code from Fortran
 */
c_int update_orbs_s_gek_factory(
    update_orbs_fp update_orbs_orig_c_funptr,
    c_int n_param,
    s_gek_settings_type settings,
    update_orbs_fp* update_orbs_s_gek_c_funptr
);

/**
 * Fortran-callable update_orbs_s_gek_deconstructor interface
 */
void update_orbs_s_gek_deconstructor(void);

#ifdef __cplusplus
}
#endif

/* ------------------------------------------------------------------
 * Small C helper functions to mimic Fortran settings%init()
 * ------------------------------------------------------------------ */

static inline s_gek_settings_type s_gek_settings_init(void) {
    s_gek_settings_type s = {0};
    init_s_gek_settings(&s);
    return s;
}

#endif
