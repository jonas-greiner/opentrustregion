// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENTRUSTREGION_MO_H
#define OPENTRUSTREGION_MO_H

#include "opentrustregion.h"
#include "opentrustregion_common.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------
 * Struct corresponding to Fortran type(mo_settings_type_c)
 * ------------------------------------------------------------------ */
typedef struct {
  logger_fp logger;
  c_bool initialized;
  c_int verbose;
} mo_settings_type;

// Fortran-callable init routine for MO settings
void init_mo_settings(mo_settings_type *settings);

/* ------------------------------------------------------------------
 * Fortran wrappers
 * ------------------------------------------------------------------ */

/**
 * Fortran-callable MO factory interface, whose parameters are the occupied-virtual
 * rotations of every particle channel.
 *
 * @param mo_coeff_c                              Flattened MO coefficients, occupied
 *                                                first, column-major per particle
 *                                                channel (size n_ao * n_mo *
 *                                                n_particle); rotated in place
 * @param ao_overlap_c                            Flattened AO overlap matrix (size
 *                                                n_ao^2)
 * @param n_occ_c                                 Number of occupied orbitals of every
 *                                                particle channel (size n_particle)
 * @param n_particle_c                            Number of particles
 * @param n_ao_c                                  Number of AO basis functions
 * @param n_mo_c                                  Number of MOs
 * @param evaluate_dm_c_funptr                    C pointer to evaluate_dm callback
 * @param obj_func_mo_c_funptr                    Output: wrapped objective function
 *                                                pointer
 * @param update_orbs_mo_c_funptr                 Output: wrapped update_orbs function
 *                                                pointer
 * @param solver_settings_c                       Input/output: solver settings
 * @param settings_c                              MO settings
 *
 * @return                                        Integer error code from Fortran
 */
c_int mo_factory(c_real *mo_coeff_c, const c_real *ao_overlap_c, const c_int *n_occ_c,
                 c_int n_particle_c, c_int n_ao_c, c_int n_mo_c,
                 evaluate_dm_fp evaluate_dm_c_funptr, obj_func_fp *obj_func_mo_c_funptr,
                 update_orbs_fp *update_orbs_mo_c_funptr,
                 solver_settings_type *solver_settings_c, mo_settings_type *settings_c);

/**
 * Fortran-callable MO deconstructor.
 */
void mo_deconstructor();

#ifdef __cplusplus
}
#endif

/* ------------------------------------------------------------------
 * Small C helper functions to mimic Fortran settings%init()
 * ------------------------------------------------------------------ */

static inline mo_settings_type mo_settings_init(void) {
  mo_settings_type s = {0};
  init_mo_settings(&s);
  return s;
}

#endif
