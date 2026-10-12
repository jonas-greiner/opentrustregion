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
 * Declarations for S-GEK functions and function pointer types
 * ------------------------------------------------------------------ */

/* Change of reference callback */
typedef c_int change_reference_fn(const c_real *new_ref_c, c_int n_points_c,
                                  c_real *kappa_list_c, c_real *local_grad_list_c,
                                  c_real *grad_list_c, void *context);
typedef change_reference_fn *change_reference_fp;

/* ------------------------------------------------------------------
 * Struct corresponding to Fortran type(s_gek_settings_type_c)
 * ------------------------------------------------------------------ */
typedef struct {
  logger_fp logger;
  c_bool initialized;
  c_bool use_subspace;
  c_int verbose;
  c_int max_points;
  void *context;
} s_gek_settings_type;

// Fortran-callable init routine for S-GEK settings
void init_s_gek_settings(s_gek_settings_type *settings);

/* ------------------------------------------------------------------
 * Fortran wrappers
 * ------------------------------------------------------------------ */

/**
 * Fortran-callable S-GEK factory interface.
 *
 * @param update_orbs_orig_c_funptr               C pointer to original update_orbs
 *                                                callback
 * @param change_reference_c_funptr               C pointer to change of reference
 *                                                callback
 * @param n_param_c                               Number of parameters
 * @param settings_c                              S-GEK settings
 * @param update_orbs_s_gek_c_funptr              Output: wrapped update_orbs function
 *                                                pointer
 *
 * @return                                        Integer error code from Fortran
 */
c_int s_gek_factory(update_orbs_fp update_orbs_orig_c_funptr,
                    change_reference_fp change_reference_c_funptr, c_int n_param_c,
                    const s_gek_settings_type *settings_c,
                    update_orbs_fp *update_orbs_s_gek_c_funptr);

/**
 * Fortran-callable S-GEK deconstructor.
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
