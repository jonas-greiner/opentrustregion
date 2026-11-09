// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/* Pure-C integration tests for the public C interface declared in
 * opentrustregion_rotation_qn.h
 *
 * test_rotation_qn_factory_c solves the Hartmann 6D problem with the real factory and
 * solver through the header, with either updating scheme.
 *
 * The Fortran-side rotation_qn_c_interface_unit_tests cover the bind(C) wrappers but
 * never compile against the C header itself. These tests do, and
 * test_rotation_qn_settings_layout checks that the settings struct has the size of the
 * bind(C) type and that every field is read back under its own name, so any drift
 * between the bind(C) type (Fortran) and the struct (C) is caught here, except a field
 * left in the struct that fits into the padding of the bind(C) type, which changes
 * neither the size nor the offset of a named field. The header's inline
 * initialization, which test_rotation_qn_settings_init checks, takes both sides from
 * the same Fortran assignment, so it cannot detect such drift. */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "opentrustregion.h"
#include "opentrustregion_rotation_qn.h"

#include "../../common/tests/common_c_integration_tests.h"

/* ------------------------------------------------------------------
 * Reference settings provided by rotation_qn_test_reference.f90: the reference value
 * of a field by its name, C settings filled with these values field by field, the size
 * of the bind(C) settings type, and a comparison with the default settings
 * ------------------------------------------------------------------ */

void reference_rotation_qn_field(const char *name, c_real *value, char *keyword);
void get_reference_rotation_qn_values(rotation_qn_settings_type *settings);
size_t rotation_qn_settings_size(void);
bool is_default_rotation_qn_settings(const rotation_qn_settings_type *settings);

/* ------------------------------------------------------------------
 * Host of the Hartmann 6D problem the factory test solves: the context of the solver
 * settings, which the orbital update, the objective function and the exact Hessian
 * the orbital update returns receive, holds the current point, its Hessian and the
 * Hessian diagonal there, and that of the rotation QN settings, which the initial
 * Hessian, its preconditioner, the objective function and the transport function of
 * the settings receive, refers to it; both count the calls of the callback functions
 * ------------------------------------------------------------------ */

#define N_PARAM HARTMANN6D_N_PARAM

enum { SOLVER_HOST = 1, SETTINGS_HOST = 2 };

typedef struct {
  int kind;
  c_real curr_vars[N_PARAM];
  c_real hess[N_PARAM][N_PARAM];
  int n_update_orbs, n_obj_func, n_exact_hess_x;
} solver_host_type;

typedef struct {
  int kind;
  solver_host_type *solver_host;
  int n_init_hess, n_init_precond, n_obj_func_grad, n_transport;
} settings_host_type;

static solver_host_type *solver_host(void *context) {
  solver_host_type *host = context;
  return host && host->kind == SOLVER_HOST ? host : NULL;
}

static settings_host_type *settings_host(void *context) {
  settings_host_type *host = context;
  return host && host->kind == SETTINGS_HOST ? host : NULL;
}

static c_int exact_hess_x(const c_real *x, c_real *hx, void *context) {
  solver_host_type *host = solver_host(context);
  if (!host)
    return 1;
  host->n_exact_hess_x++;
  for (int i = 0; i < N_PARAM; i++) {
    hx[i] = 0.0;
    for (int j = 0; j < N_PARAM; j++)
      hx[i] += host->hess[i][j] * x[j];
  }
  return 0;
}

static c_int update_orbs(const c_real *kappa, c_real *func, c_real *grad,
                         c_real *h_diag, hess_x_fp *hess_x_ptr, void *context) {
  solver_host_type *host = solver_host(context);
  if (!host)
    return 1;
  host->n_update_orbs++;
  for (int i = 0; i < N_PARAM; i++)
    host->curr_vars[i] += kappa[i];
  *func = hartmann_func(host->curr_vars);
  hartmann_grad(host->curr_vars, grad);
  hartmann_hess(host->curr_vars, host->hess);
  for (int i = 0; i < N_PARAM; i++)
    h_diag[i] = host->hess[i][i];
  *hess_x_ptr = exact_hess_x;
  return 0;
}

static c_int obj_func(const c_real *kappa, c_real *func, void *context) {
  solver_host_type *host = solver_host(context);
  if (!host)
    return 1;
  host->n_obj_func++;
  c_real x[N_PARAM];
  for (int i = 0; i < N_PARAM; i++)
    x[i] = host->curr_vars[i] + kappa[i];
  *func = hartmann_func(x);
  return 0;
}

static c_int init_hess(const c_real *x, c_real *hx, void *context) {
  settings_host_type *host = settings_host(context);
  if (!host)
    return 1;
  host->n_init_hess++;
  for (int i = 0; i < N_PARAM; i++)
    hx[i] = host->solver_host->hess[i][i] * x[i];
  return 0;
}

static c_int init_precond(const c_real *residual, const c_real *mu,
                          c_real *precond_residual, void *context) {
  settings_host_type *host = settings_host(context);
  if (!host)
    return 1;
  host->n_init_precond++;
  for (int i = 0; i < N_PARAM; i++) {
    c_real divisor = host->solver_host->hess[i][i] - *mu;
    if (fabs(divisor) < 1e-10)
      divisor = 1e-10;
    precond_residual[i] = residual[i] / divisor;
  }
  return 0;
}

static c_int obj_func_grad(const c_real *kappa, c_real *grad, c_real *func,
                           void *context) {
  settings_host_type *host = settings_host(context);
  if (!host)
    return 1;
  host->n_obj_func_grad++;
  c_real x[N_PARAM];
  for (int i = 0; i < N_PARAM; i++)
    x[i] = host->solver_host->curr_vars[i] + kappa[i];
  *func = hartmann_func(x);
  hartmann_grad(x, grad);
  return 0;
}

static c_int transport(const c_real *geodesic, c_real *tangent_vector, void *context) {
  settings_host_type *host = settings_host(context);
  if (!host)
    return 1;
  host->n_transport++;
  return 0;
}

/* ------------------------------------------------------------------
 * Tests
 * ------------------------------------------------------------------ */

bool test_rotation_qn_factory_c(void) {
  bool ok = true;
  const char *qn_methods[2] = {"sr1", "bfgs"};

  /* solve the Hartmann 6D problem from near its minimum with either updating scheme,
   * SR1 with every callback function of the settings, the Hessian diagonal as the
   * initial Hessian, its level-shifted inverse as its preconditioner, the objective
   * function of the settings, which has to replace the host's, and an identity as
   * the transport function, and BFGS with none, whose defaults are the same apart
   * from the objective function, so that the host's is called, and check that the
   * solve converges to the minimum, that every callback function is called with its
   * context, while the exact Hessian the orbital update returns is not, that the
   * factory lets the solver refresh the model after a rejected step only with the
   * objective function of the settings and that the object is released */
  for (int i_method = 0; i_method < 2; i_method++) {
    solver_host_type host = {SOLVER_HOST};
    memcpy(host.curr_vars, hartmann6d_near_minimum, sizeof(host.curr_vars));
    settings_host_type settings_ctx = {SETTINGS_HOST, &host};

    solver_settings_type settings = solver_settings_init();
    settings.context = &host;
    rotation_qn_settings_type rotation_qn_settings = rotation_qn_settings_init();
    strcpy(rotation_qn_settings.qn_method, qn_methods[i_method]);
    bool callbacks = i_method == 0;
    if (callbacks) {
      rotation_qn_settings.init_hess = init_hess;
      rotation_qn_settings.init_precond = init_precond;
      rotation_qn_settings.obj_func_grad = obj_func_grad;
      rotation_qn_settings.transport = transport;
    }
    rotation_qn_settings.context = &settings_ctx;
    solver_state_type state;

    update_orbs_fp update_orbs_rotation_qn;
    obj_func_fp obj_func_rotation_qn;
    void *extension = NULL;
    c_int error = rotation_qn_factory(
        N_PARAM, update_orbs, obj_func, &update_orbs_rotation_qn, &obj_func_rotation_qn,
        &extension, &settings, &state, &rotation_qn_settings);
    if (error != 0) {
      fprintf(stderr,
              "test_rotation_qn_factory_c failed: Factory produced error for "
              "%s.\n",
              qn_methods[i_method]);
      ok = false;
      continue;
    }
    if (settings.refresh_hess != callbacks) {
      fprintf(stderr,
              "test_rotation_qn_factory_c failed: Refresh of the model after "
              "a rejected step not set with exactly the objective function of the "
              "settings for %s.\n",
              qn_methods[i_method]);
      ok = false;
    }

    error = solver(update_orbs_rotation_qn, obj_func_rotation_qn, N_PARAM, &settings,
                   &state);
    /* keep a deviation that is not a number, which fmax would drop */
    c_real distance = 0.0;
    for (int i = 0; i < N_PARAM; i++) {
      c_real deviation = fabs(host.curr_vars[i] - hartmann6d_minimum1[i]);
      if (!(deviation <= distance))
        distance = deviation;
    }
    if (error != 0 || !(distance < 1e-4)) {
      fprintf(stderr,
              "test_rotation_qn_factory_c failed: Solve did not converge to "
              "the minimum for %s.\n",
              qn_methods[i_method]);
      ok = false;
    }
    if (host.n_update_orbs == 0 || (settings_ctx.n_init_hess > 0) != callbacks ||
        (settings_ctx.n_init_precond > 0) != callbacks ||
        (settings_ctx.n_obj_func_grad > 0) != callbacks ||
        (settings_ctx.n_transport > 0) != callbacks ||
        (host.n_obj_func > 0) == callbacks || host.n_exact_hess_x != 0) {
      fprintf(stderr,
              "test_rotation_qn_factory_c failed: Callback functions not "
              "called with their contexts, host's objective function called "
              "together with the one of the settings or exact Hessian called "
              "for %s.\n",
              qn_methods[i_method]);
      ok = false;
    }

    release_extension(extension);
  }

  return ok;
}

bool test_rotation_qn_settings_layout(void) {
  /* check that the settings struct has the size of the bind(C) settings type, so that
   * a field missing from either is detected, which the field checks below would not
   * read, and stop otherwise, since filling smaller settings would write past their
   * end */
  if (sizeof(rotation_qn_settings_type) != rotation_qn_settings_size()) {
    fprintf(stderr, "test_rotation_qn_settings_layout failed: Size of rotation QN "
                    "settings differs "
                    "from the library's.\n");
    return false;
  }

  /* check that every field is read back under its own name */
  bool ok = true;
  rotation_qn_settings_type s = {0};
  get_reference_rotation_qn_values(&s);
  reference_field_fn *lookup = reference_rotation_qn_field;
  const char *test_name = "test_rotation_qn_settings_layout";
  ok &= check_settings_field(lookup, (uintptr_t)s.init_hess, "init_hess", test_name);
  ok &= check_settings_field(lookup, (uintptr_t)s.init_precond, "init_precond",
                             test_name);
  ok &= check_settings_field(lookup, (uintptr_t)s.obj_func_grad, "obj_func_grad",
                             test_name);
  ok &= check_settings_field(lookup, (uintptr_t)s.transport, "transport", test_name);
  ok &= check_settings_field(lookup, (uintptr_t)s.logger, "logger", test_name);
  ok &= check_settings_field(lookup, s.initialized, "initialized", test_name);
  ok &= check_settings_field(lookup, s.verbose, "verbose", test_name);
  ok &= check_settings_field(lookup, s.max_points, "max_points", test_name);
  ok &= check_settings_keyword(lookup, s.qn_method, "qn_method", test_name);
  ok &= check_settings_field(lookup, (uintptr_t)s.context, "context", test_name);

  return ok;
}

bool test_rotation_qn_settings_init(void) {
  bool ok = true;

  // call function and compare with the default settings without callback functions and
  // host context
  rotation_qn_settings_type s = rotation_qn_settings_init();
  if (!is_default_rotation_qn_settings(&s)) {
    fprintf(stderr,
            "test_rotation_qn_settings_init failed: Settings not initialized to the "
            "default values without callback functions and host context.\n");
    ok = false;
  }

  return ok;
}
