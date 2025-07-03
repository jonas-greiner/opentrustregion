// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifdef __cplusplus
extern "C" {
#endif

void solver(
    const void (*update_orbs)(const double*, double*, void**, void**, void (**)(const double*, void**)),
    const double (*obj_func)(const double*),
    const int n_param,
    const bool* error,
    const void* precond_ptr,
    const void (*conv_check)(),
    const void* stability_ptr,
    const void* line_search_ptr,
    const void* davidson_c_ptr,
    const void* jacobi_davidson_c_ptr,
    const void* prefer_jacobi_davidson,
    const void* conv_tol_ptr,
    const void* n_random_trial_vectors_ptr,
    const void* start_trust_radius_ptr,
    const void* n_macro_ptr,
    const void* n_micro_ptr,
    const void* global_red_factor_ptr,
    const void* local_red_factor_ptr,
    const void* seed_ptr,
    const void* verbose_ptr,
    const void (*logger)(const char*)
);

void stability_check(
    const double* h_diag,
    const void (*hess_x)(const double*, void**),
    const int n_param,
    bool stable,
    double* kappa,
    const bool* error,
    const void* precond_ptr,
    const void* jacobi_davidson_ptr,
    const void* conv_tol_ptr,
    const void* n_random_trial_vectors_ptr,
    const void* n_iter_ptr,
    const void* verbose_ptr,
    const void (*logger)(const char*)
);

#ifdef __cplusplus
}
#endif
