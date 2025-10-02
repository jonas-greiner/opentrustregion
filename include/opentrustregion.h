// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

int64_t solver(
    const int64_t (*update_orbs)(const double*, double*, double**, double**, int64_t (**)(const double*, double**)),
    const int64_t (*obj_func)(const double*, double*),
    const int64_t n_param,
    const int64_t (*precond)(const double*, const double, double**),
    const int64_t (*conv_check)(bool*),
    const bool* stability_ptr,
    const bool* line_search_ptr,
    const bool* davidson_ptr,
    const bool* jacobi_davidson_ptr,
    const bool* prefer_jacobi_davidson,
    const double* conv_tol_ptr,
    const int64_t* n_random_trial_vectors_ptr,
    const double* start_trust_radius_ptr,
    const int64_t* n_macro_ptr,
    const int64_t* n_micro_ptr,
    const bool* global_red_factor_ptr,
    const bool* local_red_factor_ptr,
    const int64_t* seed_ptr,
    const int64_t* verbose_ptr,
    const void (*logger)(const char*)
);

int64_t stability_check(
    const double* h_diag,
    const int64_t (*hess_x)(const double*, double**),
    const int64_t n_param,
    bool stable,
    double* kappa,
    const int64_t (*precond)(const double*, const double, double**),
    const bool* jacobi_davidson_ptr,
    const double* conv_tol_ptr,
    const int64_t* n_random_trial_vectors_ptr,
    const int64_t* n_iter_ptr,
    const int64_t* verbose_ptr,
    const void (*logger)(const char*)
);

#ifdef __cplusplus
}
#endif
