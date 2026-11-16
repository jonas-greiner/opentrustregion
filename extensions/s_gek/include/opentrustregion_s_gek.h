// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENTRUSTREGION_S_GEK_H
#define OPENTRUSTREGION_S_GEK_H

#include "opentrustregion.h"
#include "opentrustregion_common.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------
 * Declarations for S-GEK functions and function pointer types
 * ------------------------------------------------------------------ */

/* Chart transformation callback */
typedef c_int change_reference_fn(const c_real *new_ref_c, c_int n_points_c,
                                  c_real *kappa_list_c, c_real *local_grad_list_c,
                                  c_real *grad_list_c, void *context);
typedef change_reference_fn *change_reference_fp;

/* ------------------------------------------------------------------
 * Struct corresponding to Fortran type(s_gek_settings_type_c)
 * ------------------------------------------------------------------ */
typedef struct {
  obj_func_grad_fp obj_func_grad;
  hess_x_fp init_hess;
  logger_fp logger;
  c_bool initialized;
  c_bool excess_dispersion;
  c_int verbose;
  char inner_solver[OTR_KW_LEN + 1];
  char prior_update[OTR_KW_LEN + 1];
  char acceptance[OTR_KW_LEN + 1];
  void *context;
} s_gek_settings_type;

// Fortran-callable init routine for S-GEK settings
void init_s_gek_settings(s_gek_settings_type *settings);

/* ------------------------------------------------------------------
 * Fortran wrappers
 * ------------------------------------------------------------------ */

/**
 * C interface to the S-GEK solver
 *
 * @param update_orbs_ptr      Pointer to update_orbs callback
 * @param obj_func_ptr         Pointer to objective function callback
 * @param change_reference_ptr Pointer to the chart transformation callback, which
 *                             receives the host context of the S-GEK settings
 * @param n_param              Number of parameters
 * @param settings             Pointer to struct of solver settings, may be NULL to use
 *                             the default settings
 * @param state                Pointer to struct that receives the state of the
 *                             optimization, may be NULL to skip returning it
 * @param s_gek_settings       Pointer to struct of S-GEK settings, may be NULL to use
 *                             the default settings
 * @return                     Integer error code from Fortran
 */
c_int s_gek_solver(update_orbs_fp update_orbs_ptr, obj_func_fp obj_func_ptr,
                   change_reference_fp change_reference_ptr, c_int n_param,
                   const solver_settings_type *settings, solver_state_type *state,
                   const s_gek_settings_type *s_gek_settings);

#ifdef __cplusplus
}
#endif

/* ------------------------------------------------------------------
 * Small C helper functions returning the default settings, which a declared Fortran
 * settings object holds
 * ------------------------------------------------------------------ */

static inline s_gek_settings_type s_gek_settings_init(void) {
  s_gek_settings_type s;
  init_s_gek_settings(&s);
  return s;
}

#endif
