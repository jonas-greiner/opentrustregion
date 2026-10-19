// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENTRUSTREGION_ROTATION_QN_H
#define OPENTRUSTREGION_ROTATION_QN_H

#include "opentrustregion.h"
#include "opentrustregion_common.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------
 * Declarations for rotation QN functions and function pointer types
 * ------------------------------------------------------------------ */

/* Transport callback */
typedef c_int transport_fn(const c_real *geodesic_c, c_real *tangent_vector_c,
                           void *context);
typedef transport_fn *transport_fp;

/* Hessian initialization callback */
typedef c_int init_hess_fn(c_real *vector_c, void *context);
typedef init_hess_fn *init_hess_fp;

/* ------------------------------------------------------------------
 * Struct corresponding to Fortran type(rotation_qn_settings_type_c)
 * ------------------------------------------------------------------ */
typedef struct {
  logger_fp logger;
  c_bool initialized;
  c_int verbose;
  c_int max_points;
  char qn_method[OTR_KW_LEN + 1];
  void *context;
} rotation_qn_settings_type;

// Fortran-callable init routine for rotation QN settings
void init_rotation_qn_settings(rotation_qn_settings_type *settings);

/* ------------------------------------------------------------------
 * Fortran wrappers
 * ------------------------------------------------------------------ */

/**
 * Fortran-callable rotation QN factory interface.
 *
 * @param n_param_c                               Number of parameters
 * @param update_orbs_c_funptr                    C pointer to update_orbs callback
 * @param obj_func_c_funptr                       C pointer to obj_func callback
 * @param transport_c_funptr                      C pointer to transport callback
 * @param init_hess_c_funptr                      C pointer to init_hess callback
 * @param update_orbs_rotation_qn_c_funptr        Output: wrapped update_orbs function
 *                                                pointer
 * @param obj_func_rotation_qn_c_funptr           Output: wrapped objective function
 *                                                pointer
 * @param extension_c                             Output: rotation QN extension,
 *                                                released with release_extension
 * @param solver_settings_c                       Input/output: solver settings, whose
 *                                                context the returned functions have
 *                                                to be called with
 * @param solver_state_c                          Solver state the host also passes to
 *                                                the solver, which has to outlive the
 *                                                extension
 * @param settings_c                              Rotation QN settings, may be NULL
 *                                                to use the default settings
 *
 * @return                                        Integer error code from Fortran
 */
c_int rotation_qn_factory(c_int n_param_c, update_orbs_fp update_orbs_c_funptr,
                          obj_func_fp obj_func_c_funptr,
                          transport_fp transport_c_funptr,
                          init_hess_fp init_hess_c_funptr,
                          update_orbs_fp *update_orbs_rotation_qn_c_funptr,
                          obj_func_fp *obj_func_rotation_qn_c_funptr,
                          void **extension_c, solver_settings_type *solver_settings_c,
                          solver_state_type *solver_state_c,
                          const rotation_qn_settings_type *settings_c);

#ifdef __cplusplus
}
#endif

/* ------------------------------------------------------------------
 * Small C helper functions to mimic Fortran settings%init()
 * ------------------------------------------------------------------ */

static inline rotation_qn_settings_type rotation_qn_settings_init(void) {
  rotation_qn_settings_type s = {0};
  init_rotation_qn_settings(&s);
  return s;
}

#endif
