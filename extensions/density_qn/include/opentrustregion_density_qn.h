// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENTRUSTREGION_DENSITY_QN_H
#define OPENTRUSTREGION_DENSITY_QN_H

#include "opentrustregion.h"
#include "opentrustregion_mo.h"
#include "opentrustregion_oao.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------
 * Declarations for density QN functions and function pointer types
 * ------------------------------------------------------------------ */

/* Density matrix evaluating callback with Coulomb, exact-exchange and non-linear
 * potential contributions */
typedef c_int density_qn_evaluate_dm_fn(const c_real *dm_ao_c, c_real *energy_c,
                                        c_real *fock_c, c_real *v_coulomb_c,
                                        c_real *v_exchange_c, c_real *v_nonlinear_c,
                                        void *context);
typedef density_qn_evaluate_dm_fn *density_qn_evaluate_dm_fp;

/* ------------------------------------------------------------------
 * Struct corresponding to Fortran type(density_qn_settings_type_c)
 * ------------------------------------------------------------------ */
typedef struct {
  logger_fp logger;
  c_bool initialized;
  c_int verbose;
  char qn_method[OTR_KW_LEN + 1];
  void *context;
} density_qn_settings_type;

// Fortran-callable init routine for density QN settings
void init_density_qn_settings(density_qn_settings_type *settings);

/* ------------------------------------------------------------------
 * Fortran wrappers
 * ------------------------------------------------------------------ */

/**
 * Fortran-callable density QN factory interface, which builds on the MO or OAO
 * extension whose factory wired the solver settings last and moves its orbitals.
 *
 * @param evaluate_dm_c_funptr            C pointer to evaluate_dm callback of the
 *                                        shell of the orbitals
 * @param update_orbs_density_qn_c_funptr Output: wrapped update_orbs function pointer
 * @param obj_func_density_qn_c_funptr    Output: wrapped objective function pointer
 * @param extension_c                     Output: density QN extension, released with
 *                                        release_extension before the MO or OAO
 *                                        extension
 * @param solver_settings_c               Input/output: solver settings the MO or OAO
 *                                        factory wired last, whose context the
 *                                        returned functions have to be called with
 * @param solver_state_c                  Solver state the host also passes to the
 *                                        solver, which has to outlive the extension
 * @param settings_c                      Density QN settings, may be NULL to use the
 *                                        default settings
 *
 * @return                           Integer error code from Fortran
 */
c_int density_qn_factory(density_qn_evaluate_dm_fp evaluate_dm_c_funptr,
                         update_orbs_fp *update_orbs_density_qn_c_funptr,
                         obj_func_fp *obj_func_density_qn_c_funptr, void **extension_c,
                         solver_settings_type *solver_settings_c,
                         solver_state_type *solver_state_c,
                         const density_qn_settings_type *settings_c);

#ifdef __cplusplus
}
#endif

/* ------------------------------------------------------------------
 * Small C helper functions returning the default settings, which a declared Fortran
 * settings object holds
 * ------------------------------------------------------------------ */

static inline density_qn_settings_type density_qn_settings_init(void) {
  density_qn_settings_type s;
  init_density_qn_settings(&s);
  return s;
}

#endif
