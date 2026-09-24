// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENTRUSTREGION_OAO_H
#define OPENTRUSTREGION_OAO_H

#include "opentrustregion.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------
 * Declarations for OAO functions and function pointer types
 * ------------------------------------------------------------------ */

/* Response callback */
typedef c_int get_response_fn(const c_real *dm_ao_c, c_real *response_c);
typedef get_response_fn *get_response_fp;

/* Density matrix evaluating callback */
typedef c_int evaluate_dm_fn(const c_real *dm_ao_c, c_real *energy_c, c_real *fock_c,
                             get_response_fp *get_response_ptr);
typedef evaluate_dm_fn *evaluate_dm_fp;

/* ------------------------------------------------------------------
 * Struct corresponding to Fortran type(oao_settings_type_c)
 * ------------------------------------------------------------------ */
typedef struct {
  logger_fp logger;
  c_bool initialized;
  c_int verbose;
} oao_settings_type;

// Fortran-callable init routine for OAO settings
void init_oao_settings(oao_settings_type *settings);

/* ------------------------------------------------------------------
 * Fortran wrappers
 * ------------------------------------------------------------------ */

/**
 * Fortran-callable OAO factory interface.
 *
 * @param dm_ao_c                                 Flattened AO density matrix (size
 *                                                n_ao^2)
 * @param ao_overlap_c                            Flattened AO overlap matrix (size
 *                                                n_ao^2)
 * @param n_particle_c                            Number of particles
 * @param n_ao_c                                  Number of AO basis functions
 * @param evaluate_dm_c_funptr                    C pointer to evaluate_dm callback
 * @param settings_c                              OAO settings
 * @param obj_func_oao_c_funptr                   Output: wrapped objective function
 *                                                pointer
 * @param update_orbs_oao_c_funptr                Output: wrapped update_orbs function
 *                                                pointer
 * @param precond_oao_c_funptr                    Output: wrapped level-shifted
 *                                                preconditioner function pointer
 * @param precond_pd_oao_c_funptr                 Output: wrapped positive-definite
 *                                                preconditioner function pointer
 * @param project_oao_c_funptr                    Output: wrapped projection function
 *                                                pointer
 * @param get_extra_trial_vectors_oao_c_funptr    Output: wrapped extra trial vector
 *                                                function pointer
 *
 * @return                                        Integer error code from Fortran
 */
c_int oao_factory(const c_real *dm_ao_c, const c_real *ao_overlap_c, c_int n_particle_c,
                  c_int n_ao_c, evaluate_dm_fp evaluate_dm_c_funptr,
                  obj_func_fp *obj_func_oao_c_funptr,
                  update_orbs_fp *update_orbs_oao_c_funptr,
                  precond_fp *precond_oao_c_funptr,
                  precond_pd_fp *precond_pd_oao_c_funptr,
                  project_fp *project_oao_c_funptr,
                  get_extra_trial_vectors_fp *get_extra_trial_vectors_oao_c_funptr,
                  oao_settings_type *settings_c);

/**
 * Fortran-callable OAO deconstructor.
 */
void oao_deconstructor();

#ifdef __cplusplus
}
#endif

/* ------------------------------------------------------------------
 * Small C helper functions to mimic Fortran settings%init()
 * ------------------------------------------------------------------ */

static inline oao_settings_type oao_settings_init(void) {
  oao_settings_type s = {0};
  init_oao_settings(&s);
  return s;
}

#endif
