// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENTRUSTREGION_OAO_H
#define OPENTRUSTREGION_OAO_H

#include "opentrustregion.h"

#ifdef __cplusplus
extern "C"
{
#endif

    /* ------------------------------------------------------------------
     * Declarations for OAO functions and function pointer types
     * ------------------------------------------------------------------ */

    /* Energy callback */
    typedef c_int get_energy_fn(const c_real *dm_ao_c, c_real *energy_c);
    typedef get_energy_fn *get_energy_fp;

    /* Response callback */
    typedef c_int get_response_fn(const c_real *dm_ao_c, c_real *response_c);
    typedef get_response_fn *get_response_fp;

    /* Density matrix updating callback */
    typedef c_int update_dm_fn(
        const c_real *dm_ao_c,
        c_real *energy_c,
        c_real *fock_c,
        get_response_fp *get_response_ptr);
    typedef update_dm_fn *update_dm_fp;

    /* ------------------------------------------------------------------
     * Struct corresponding to Fortran type(oao_settings_type_c)
     * ------------------------------------------------------------------ */
    typedef struct
    {
        logger_fp logger;
        c_bool initialized;
        c_bool restricted;
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
     * @param dm_ao_c                    Flattened AO density matrix (size n_ao^2)
     * @param ao_overlap_c               Flattened AO overlap matrix (size n_ao^2)
     * @param n_particle_c               Number of particles
     * @param n_ao_c                     Number of AO basis functions
     * @param get_energy_c_funptr        C pointer to get_energy callback
     * @param update_dm_c_funptr         C pointer to update_dm callback
     * @param settings_c                 OAO settings
     * @param obj_func_oao_c_funptr      Output: wrapped objective function pointer
     * @param update_orbs_oao_c_funptr   Output: wrapped update_orbs function pointer
     * @param precond_oao_c_funptr       Output: wrapped preconditioner function pointer
     *
     * @return                           Integer error code from Fortran
     */
    c_int oao_factory(
        const c_real *dm_ao_c,
        const c_real *ao_overlap_c,
        c_int n_particle_c,
        c_int n_ao_c,
        get_energy_fp get_energy_c_funptr,
        update_dm_fp update_dm_c_funptr,
        obj_func_fp *obj_func_oao_c_funptr,
        update_orbs_fp *update_orbs_oao_c_funptr,
        precond_fp *precond_oao_c_funptr,
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

static inline oao_settings_type oao_settings_init(void)
{
    oao_settings_type s = {0};
    init_oao_settings(&s);
    return s;
}

#endif
