// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENTRUSTREGION_ARH_H
#define OPENTRUSTREGION_ARH_H

#include "opentrustregion.h"
#include "opentrustregion_oao.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------
 * Declarations for ARH functions and function pointer types
 * ------------------------------------------------------------------ */

/* Density matrix updating callback with Coulomb and exchange contributions */
typedef c_int update_dm_jk_fn(
    const c_real* dm_ao_c, 
    c_real* energy_c, 
    c_real* fock_c, 
    c_real* coulomb_c,
    c_real* exchange_c,
    get_response_fp* get_response_ptr
);
typedef update_dm_jk_fn* update_dm_jk_fp;
typedef void (*generic_fp)(void);

/* ------------------------------------------------------------------
 * Struct corresponding to Fortran type(arh_settings_type_c)
 * ------------------------------------------------------------------ */
typedef struct {
    logger_fp logger;
    c_bool initialized;
    c_bool restricted;
    c_int verbose;
    char arh_type[OTR_KW_LEN + 1];
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
 * @param update_dm_c_funptr         C pointer to update_dm or update_dm_jk callback
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
    generic_fp update_dm_c_funptr,
    obj_func_fp* obj_func_arh_c_funptr,
    update_orbs_fp* update_orbs_arh_c_funptr,
    precond_fp* precond_arh_c_funptr,
    arh_settings_type *settings_c
);

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
