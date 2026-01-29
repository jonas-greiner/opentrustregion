// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENTRUSTREGION_ARH_H
#define OPENTRUSTREGION_ARH_H

#include "opentrustregion.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------
 * Declarations for ARH functions and function pointer types
 * ------------------------------------------------------------------ */

/* Energy callback */
typedef c_int get_energy_fn(const c_real* dm_ao_c, c_real* energy_c);
typedef get_energy_fn* get_energy_fp;

/* Fock matrix callback */
typedef c_int get_fock_fn(const c_real* dm_ao_c, c_real* energy_c, c_real* fock_c);
typedef get_fock_fn* get_fock_fp;

/* ------------------------------------------------------------------
 * Struct corresponding to Fortran type(arh_settings_type_c)
 * ------------------------------------------------------------------ */
typedef struct {
    logger_fp logger;
    c_bool initialized;
    c_bool restricted;
    c_int verbose;
} arh_settings_type;

// Fortran-callable init routine for ARH settings
void init_arh_settings(arh_settings_type* settings);

/* ------------------------------------------------------------------
 * Fortran wrappers
 * ------------------------------------------------------------------ */

/**
 * Fortran-callable ARH factory interface.
 *
 * @param dm_ao_c                    Flattened AO density matrix (size n_ao^2)
 * @param ao_overlap_c               Flattened AO overlap matrix (size n_ao^2)
 * @param n_particle_c               Number of particles
 * @param n_ao_c                     Number of AO basis functions
 * @param get_energy_c_funptr        C pointer to get_energy callback
 * @param get_fock_c_funptr          C pointer to get_fock callback
 * @param settings_c                 ARH settings
 * @param obj_func_arh_c_funptr      Output: wrapped objective function pointer
 * @param update_orbs_arh_c_funptr   Output: wrapped update_orbs function pointer
 * @param precond_arh_c_funptr       Output: wrapped preconditioner function pointer
 *
 * @return                           Integer error code from Fortran
 */
c_int arh_factory(
    const c_real* dm_ao_c,
    const c_real* ao_overlap_c,
    c_int n_particle_c,
    c_int n_ao_c,
    get_energy_fp get_energy_c_funptr,
    get_fock_fp get_fock_c_funptr,
    obj_func_fp* obj_func_arh_c_funptr,
    update_orbs_fp* update_orbs_arh_c_funptr,
    precond_fp* precond_arh_c_funptr,
    arh_settings_type settings_c
);

/**
 * Fortran-callable ARH deconstructor.
 *
 * @param dm_ao_c    Flattened AO density matrix (size n_ao^2).
 * @return           Integer error code.
 */
c_int arh_deconstructor(c_real* dm_ao_c);

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
