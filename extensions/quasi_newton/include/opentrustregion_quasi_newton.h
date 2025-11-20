// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENTRUSTREGION_QUASI_NEWTON_H
#define OPENTRUSTREGION_QUASI_NEWTON_H

#include "opentrustregion.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------
 * Struct corresponding to Fortran type(qn_settings_type_c)
 * ------------------------------------------------------------------ */
typedef struct {
    logger_fp logger;
    c_bool initialized;
    c_int verbose;
    char hess_update_scheme[OTR_KW_LEN + 1];
} qn_settings_type;

// Fortran-callable init routine for quasi-Newton settings
void init_qn_settings(qn_settings_type* settings);

/* ------------------------------------------------------------------
 * Fortran wrappers
 * ------------------------------------------------------------------ */

/**
 * Fortran-callable update_orbs_qn_factory interface
 *
 * @param update_orbs_orig_c_funptr    Pointer to original update_orbs callback
 * @param n_param                      Number of parameters
 * @param settings                     Quasi-Newton settings
 * @param update_orbs_qn_c_funptr      Pointer to new update_orbs callback
 * @return                             Integer error code from Fortran
 */
c_int update_orbs_qn_factory(
    update_orbs_fp update_orbs_orig_c_funptr,
    c_int n_param,
    qn_settings_type settings,
    update_orbs_fp* update_orbs_qn_c_funptr
);

/**
 * Fortran-callable update_orbs_qn_deconstructor interface
 */
void update_orbs_qn_deconstructor(void);

#ifdef __cplusplus
}
#endif

/* ------------------------------------------------------------------
 * Small C helper functions to mimic Fortran settings%init()
 * ------------------------------------------------------------------ */

static inline qn_settings_type qn_settings_init(void) {
    qn_settings_type s = {0};
    init_qn_settings(&s);
    return s;
}

#endif
