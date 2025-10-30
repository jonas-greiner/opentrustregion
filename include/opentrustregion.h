// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENTRUSTREGION_H
#define OPENTRUSTREGION_H

#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------
 * Type aliases matching Fortran kinds
 * ------------------------------------------------------------------ */
#ifdef USE_ILP64
typedef int64_t c_int; /* corresponds to integer(c_ip) */
#else
typedef int32_t c_int; /* corresponds to integer(c_ip) */
#endif
typedef double c_real;  /* corresponds to real(c_rp) */
typedef bool c_bool;  /* corresponds to logical(c_bool) */

/* ------------------------------------------------------------------
 * Maximum keyword length
 * ------------------------------------------------------------------ */
#define OTR_KW_LEN 64

/* ------------------------------------------------------------------
 * Forward declarations for function pointer types
 * ------------------------------------------------------------------ */

/* Hessian-vector product callback */
typedef c_int (*hess_x_type)(const c_real* x_c, c_real* hess_x);

/* Orbital update callback */
typedef c_int (*update_orbs_type)(
    const c_real* kappa,
    c_real* func,
    c_real* grad,
    c_real* h_diag,
    void** hess_x_ptr  /* corresponds to type(c_funptr) intent(out) */
);

/* Objective function callback */
typedef c_int (*obj_func_type)(const c_real* kappa, c_real* func);

/* Preconditioner callback */
typedef c_int (*precond_type)(
    const c_real* residual, const c_real* mu, c_real* precond_residual
);

/* Convergence check callback */
typedef c_int (*conv_check_type)(c_bool* converged);

/* Logger callback */
typedef void (*logger_type)(const char* message);

/* ------------------------------------------------------------------
 * Struct corresponding to Fortran type(solver_settings_type_c)
 * ------------------------------------------------------------------ */
typedef struct {
    void* precond;
    void* conv_check;
    void* logger;

    c_bool stability;
    c_bool line_search;
    c_bool initialized;

    c_real conv_tol;
    c_real start_trust_radius;
    c_real global_red_factor;
    c_real local_red_factor;

    c_int n_random_trial_vectors;
    c_int n_macro;
    c_int n_micro;
    c_int jacobi_davidson_start;
    c_int seed;
    c_int verbose;

    char subsystem_solver[OTR_KW_LEN + 1];
} solver_settings_type;


/* ------------------------------------------------------------------
 * Struct corresponding to Fortran type(stability_settings_type_c)
 * ------------------------------------------------------------------ */
typedef struct {
    void* precond;
    void* logger;

    c_bool initialized;

    c_real conv_tol;

    c_int n_random_trial_vectors;
    c_int n_iter;
    c_int jacobi_davidson_start;
    c_int seed;
    c_int verbose;

    char diag_solver[OTR_KW_LEN + 1];
} stability_settings_type;

/* ------------------------------------------------------------------
 * Fortran solver wrapper
 * ------------------------------------------------------------------ */

/**
 * Fortran-callable solver interface
 *
 * @param update_orbs_ptr   Pointer to update_orbs callback
 * @param obj_func_ptr      Pointer to objective function callback
 * @param n_param           Number of parameters
 * @param settings          Struct of solver settings
 * @return                  Integer error code from Fortran
 */
c_int solver(
    void* update_orbs_ptr, 
    void* obj_func_ptr, 
    c_int n_param, 
    solver_settings_type settings
);

/**
 * Fortran-callable stability check interface
 *
 * @param h_diag_ptr        Pointer to diagonal Hessian array
 * @param hess_x_ptr        Pointer to Hessian–vector product callback
 * @param n_param           Number of parameters
 * @param stable            Pointer to bool that receives stability result
 * @param settings          Struct of stability solver settings
 * @param kappa_ptr         Pointer to orbital rotation vector
 * @return                  Integer error code from Fortran
 */
c_int stability_check(
    const void* h_diag_ptr,
    const void* hess_x_ptr,
    c_int n_param,
    c_bool* stable,
    stability_settings_type settings,
    const void* kappa_ptr
);

// Fortran-callable init subroutine for solver settings
void init_solver_settings(solver_settings_type* settings);

// Fortran-callable init subroutine for stability check settings
void init_stability_settings(stability_settings_type* settings);

#ifdef __cplusplus
}
#endif

// Small helper functions to mimic Fortran settings%init()
static inline solver_settings_type solver_settings_init(void) {
    solver_settings_type s = {0};
    init_solver_settings(&s);
    return s;
}

static inline stability_settings_type stability_settings_init(void) {
    stability_settings_type s = {0};
    init_stability_settings(&s);
    return s;
}

#endif
