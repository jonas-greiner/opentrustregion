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
 * @param n_param_c                               Number of parameters
 * @param update_orbs_c_funptr                    C pointer to update_orbs callback
 * @param obj_func_c_funptr                       C pointer to obj_func callback
 * @param change_reference_c_funptr               C pointer to change of reference
 *                                                callback
 * @param update_orbs_s_gek_c_funptr              Output: wrapped update_orbs function
 *                                                pointer
 * @param obj_func_s_gek_c_funptr                 Output: wrapped objective function
 *                                                pointer
 * @param extension_c                             Output: S-GEK extension, released
 *                                                with release_extension
 * @param solver_settings_c                       Input/output: solver settings, whose
 *                                                context the returned functions have
 *                                                to be called with
 * @param solver_state_c                          Solver state the host also passes to
 *                                                the solver, which has to outlive the
 *                                                extension
 * @param settings_c                              S-GEK settings, may be NULL
 *                                                to use the default settings
 *
 * @return                                        Integer error code from Fortran
 */
c_int s_gek_factory(c_int n_param_c, update_orbs_fp update_orbs_c_funptr,
                    obj_func_fp obj_func_c_funptr,
                    change_reference_fp change_reference_c_funptr,
                    update_orbs_fp *update_orbs_s_gek_c_funptr,
                    obj_func_fp *obj_func_s_gek_c_funptr, void **extension_c,
                    solver_settings_type *solver_settings_c,
                    solver_state_type *solver_state_c,
                    const s_gek_settings_type *settings_c);

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
