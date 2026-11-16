// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/* Pure-C integration tests for the public C interface declared in
 * opentrustregion_s_gek.h
 *
 * test_s_gek_solver_c solves the Hartmann 6D problem with the real S-GEK solver
 * through the header.
 *
 * The Fortran-side s_gek_c_interface_unit_tests cover the bind(C) wrappers but never
 * compile against the C header itself. These tests do, and test_s_gek_settings_layout
 * checks that the settings struct has the size of the bind(C) type and that every
 * field is read back under its own name, so any drift between the bind(C) type
 * (Fortran) and the struct (C) is caught here, except a field left in the struct that
 * fits into the padding of the bind(C) type, which changes neither the size nor the
 * offset of a named field. The header's inline initialization, which
 * test_s_gek_settings_init checks, takes both sides from the same Fortran assignment,
 * so it cannot detect such drift. */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "opentrustregion.h"
#include "opentrustregion_s_gek.h"

#include "../../common/tests/common_c_integration_tests.h"

/* ------------------------------------------------------------------
 * Reference settings provided by s_gek_test_reference.f90: the reference value of a
 * field by its name, C settings filled with these values field by field, the size of
 * the bind(C) settings type, a comparison with the default settings and the warning
 * for settings that were not initialized
 * ------------------------------------------------------------------ */

void reference_s_gek_field(const char *name, c_real *value, char *keyword);
void get_reference_s_gek_values(s_gek_settings_type *settings);
size_t s_gek_settings_size(void);
bool is_default_s_gek_settings(const s_gek_settings_type *settings);
extern const char test_s_gek_settings_uninitialized_warning[];

/* ------------------------------------------------------------------
 * Host of the Hartmann 6D problem the solver test solves: the context of the solver
 * settings, which the orbital update, the objective function and the exact Hessian
 * the orbital update returns receive, holds the current point, its Hessian, the record
 * of the state and whether the record was not current in a callback function, and
 * that of the S-GEK settings, which the chart transformation, the objective function,
 * the initial Hessian and the logging function of the settings receive, refers to it;
 * both count the calls of the callback functions
 * ------------------------------------------------------------------ */

#define N_PARAM HARTMANN6D_N_PARAM

enum { SOLVER_HOST = 1, SETTINGS_HOST = 2 };

typedef struct {
  int kind;
  c_real curr_vars[N_PARAM];
  c_real hess[N_PARAM][N_PARAM];
  const solver_state_type *state;
  bool state_not_current;
  int n_update_orbs, n_obj_func;
} solver_host_type;

typedef struct {
  int kind;
  solver_host_type *solver_host;
  int n_change_reference, n_obj_func_grad, n_init_hess, n_logger;
} settings_host_type;

static solver_host_type *solver_host(void *context) {
  solver_host_type *host = context;
  return host && host->kind == SOLVER_HOST ? host : NULL;
}

static settings_host_type *settings_host(void *context) {
  settings_host_type *host = context;
  return host && host->kind == SETTINGS_HOST ? host : NULL;
}

/* record whether the record of the state, which is copied before every callback
 * function, holds the number of orbital updates so far */
static void check_state_current(solver_host_type *host) {
  if (host->state && host->state->n_update_orbs != host->n_update_orbs)
    host->state_not_current = true;
}

static c_int exact_hess_x(const c_real *x, c_real *hx, void *context) {
  solver_host_type *host = solver_host(context);
  if (!host)
    return 1;
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
  check_state_current(host);
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

/* the transformation of the flat parameter space, which moves the points to the new
 * reference and leaves their gradients local */
static void flat_chart(const c_real *new_ref, c_int n_points, c_real *kappa_list,
                       const c_real *local_grad_list, c_real *grad_list) {
  for (c_int p = 0; p < n_points; p++) {
    for (int i = 0; i < N_PARAM; i++) {
      kappa_list[p * N_PARAM + i] -= new_ref[i];
      grad_list[p * N_PARAM + i] = local_grad_list[p * N_PARAM + i];
    }
  }
}

/* the chart transformation of the flat parameter space, which receives the host
 * context of the S-GEK settings */
static c_int change_reference(const c_real *new_ref, c_int n_points, c_real *kappa_list,
                              c_real *local_grad_list, c_real *grad_list,
                              void *context) {
  settings_host_type *host = settings_host(context);
  if (!host)
    return 1;
  check_state_current(host->solver_host);
  host->n_change_reference++;
  flat_chart(new_ref, n_points, kappa_list, local_grad_list, grad_list);
  return 0;
}

/* the chart transformation of the flat parameter space without S-GEK settings, which
 * then receives no host context */
static c_int change_reference_without_context(const c_real *new_ref, c_int n_points,
                                              c_real *kappa_list,
                                              c_real *local_grad_list,
                                              c_real *grad_list, void *context) {
  if (context)
    return 1;
  flat_chart(new_ref, n_points, kappa_list, local_grad_list, grad_list);
  return 0;
}

static c_int obj_func_grad(const c_real *kappa, c_real *grad, c_real *func,
                           void *context) {
  settings_host_type *host = settings_host(context);
  if (!host)
    return 1;
  /* the record counts the evaluations at trial points before this one */
  check_state_current(host->solver_host);
  const solver_state_type *state = host->solver_host->state;
  if (state && state->n_obj_func != host->n_obj_func_grad)
    host->solver_host->state_not_current = true;
  host->n_obj_func_grad++;
  c_real x[N_PARAM];
  for (int i = 0; i < N_PARAM; i++)
    x[i] = host->solver_host->curr_vars[i] + kappa[i];
  *func = hartmann_func(x);
  hartmann_grad(x, grad);
  return 0;
}

static c_int init_hess(const c_real *x, c_real *hx, void *context) {
  settings_host_type *host = settings_host(context);
  if (!host)
    return 1;
  check_state_current(host->solver_host);
  host->n_init_hess++;
  return exact_hess_x(x, hx, host->solver_host);
}

static void logger(const char *message, void *context) {
  settings_host_type *host = settings_host(context);
  if (host)
    host->n_logger++;
}

/* ------------------------------------------------------------------
 * Tests
 * ------------------------------------------------------------------ */

bool test_s_gek_solver_c(void) {
  bool ok = true;

  /* solve the Hartmann 6D problem from near its minimum in its flat chart, whose chart
   * transformation receives the context of the S-GEK settings, with every callback
   * function of the S-GEK settings, the exact Hessian as the initial Hessian, and
   * check that the solve converges to the minimum, that every callback function is
   * called with its context and finds the record of the state current and that the
   * state is returned */
  solver_host_type host = {SOLVER_HOST};
  memcpy(host.curr_vars, hartmann6d_near_minimum, sizeof(host.curr_vars));
  settings_host_type settings_ctx = {SETTINGS_HOST, &host};
  solver_state_type state = {0};
  host.state = &state;

  solver_settings_type settings = solver_settings_init();
  settings.context = &host;
  s_gek_settings_type s_gek_settings = s_gek_settings_init();
  s_gek_settings.obj_func_grad = obj_func_grad;
  s_gek_settings.init_hess = init_hess;
  s_gek_settings.logger = logger;
  s_gek_settings.context = &settings_ctx;

  c_int error = s_gek_solver(update_orbs, obj_func, change_reference, N_PARAM,
                             &settings, &state, &s_gek_settings);
  /* keep a deviation that is not a number, which fmax would drop */
  c_real distance = 0.0;
  for (int i = 0; i < N_PARAM; i++) {
    c_real deviation = fabs(host.curr_vars[i] - hartmann6d_minimum1[i]);
    if (!(deviation <= distance))
      distance = deviation;
  }
  if (error != 0 || !(distance < 1e-4)) {
    fprintf(stderr, "test_s_gek_solver_c failed: Solve did not converge to the "
                    "minimum.\n");
    ok = false;
  }
  if (host.n_update_orbs == 0 || settings_ctx.n_change_reference == 0 ||
      settings_ctx.n_obj_func_grad == 0 || settings_ctx.n_init_hess == 0) {
    fprintf(stderr, "test_s_gek_solver_c failed: Callback functions not called with "
                    "their contexts.\n");
    ok = false;
  }
  if (host.state_not_current || state.n_update_orbs != host.n_update_orbs ||
      state.n_obj_func != settings_ctx.n_obj_func_grad) {
    fprintf(stderr, "test_s_gek_solver_c failed: State was not current in the "
                    "callbacks or not returned.\n");
    ok = false;
  }

  /* the logging function of the S-GEK settings receives the error messages of S-GEK,
   * here for S-GEK settings with an option S-GEK does not know */
  strcpy(s_gek_settings.inner_solver, "unknown");
  s_gek_settings.verbose = 1; /* error messages only */
  settings_ctx.n_logger = 0;
  error = s_gek_solver(update_orbs, obj_func, change_reference, N_PARAM, &settings,
                       NULL, &s_gek_settings);
  if (error == 0 || settings_ctx.n_logger == 0) {
    fprintf(stderr, "test_s_gek_solver_c failed: Error of S-GEK settings with an "
                    "unknown option not produced or not logged with their "
                    "context.\n");
    ok = false;
  }

  /* S-GEK settings that were not initialized are replaced by the default settings
   * with a warning on the standard error, since they provide no logging function,
   * while no S-GEK settings stand for the default settings without a warning, with
   * both of which the solve converges, its trial acceptance evaluating the energy at
   * trial points with the objective function of the solver and its chart
   * transformation receiving no host context */
  s_gek_settings_type uninitialized_settings = {0};
  const s_gek_settings_type *case_settings[2] = {NULL, &uninitialized_settings};
  const char *case_warnings[2] = {NULL, test_s_gek_settings_uninitialized_warning};
  const char *case_names[2] = {"without S-GEK settings",
                               "with S-GEK settings that were not initialized"};
  for (int i_case = 0; i_case < 2; i_case++) {
    memcpy(host.curr_vars, hartmann6d_near_minimum, sizeof(host.curr_vars));
    host.state = NULL;
    host.n_obj_func = 0;
    start_stderr_capture();
    error = s_gek_solver(update_orbs, obj_func, change_reference_without_context,
                         N_PARAM, &settings, NULL, case_settings[i_case]);
    if (!stop_stderr_capture(case_warnings[i_case])) {
      fprintf(stderr,
              "test_s_gek_solver_c failed: Warning on the standard error wrong %s.\n",
              case_names[i_case]);
      ok = false;
    }
    distance = 0.0;
    for (int i = 0; i < N_PARAM; i++) {
      c_real deviation = fabs(host.curr_vars[i] - hartmann6d_minimum1[i]);
      if (!(deviation <= distance))
        distance = deviation;
    }
    if (error != 0 || !(distance < 1e-4)) {
      fprintf(stderr,
              "test_s_gek_solver_c failed: Solve did not converge to the minimum "
              "%s.\n",
              case_names[i_case]);
      ok = false;
    }
    if (host.n_obj_func == 0) {
      fprintf(stderr,
              "test_s_gek_solver_c failed: Trial points not evaluated with the "
              "objective function of the solver %s.\n",
              case_names[i_case]);
      ok = false;
    }
  }

  return ok;
}

bool test_s_gek_settings_layout(void) {
  /* check that the settings struct has the size of the bind(C) settings type, so that
   * a field missing from either is detected, which the field checks below would not
   * read, and stop otherwise, since filling smaller settings would write past their
   * end */
  if (sizeof(s_gek_settings_type) != s_gek_settings_size()) {
    fprintf(stderr, "test_s_gek_settings_layout failed: Size of S-GEK settings differs "
                    "from the library's.\n");
    return false;
  }

  /* check that every field is read back under its own name, from settings filled
   * with zero and then with one bytes before, so that no reference value can match
   * bytes the library did not write */
  bool ok = true;
  for (int fill = 0; fill < 2; fill++) {
    s_gek_settings_type s;
    memset(&s, fill, sizeof(s));
    get_reference_s_gek_values(&s);
    reference_field_fn *lookup = reference_s_gek_field;
    const char *test_name = "test_s_gek_settings_layout";
    ok &= check_settings_field(lookup, (uintptr_t)s.obj_func_grad, "obj_func_grad",
                               test_name);
    ok &= check_settings_field(lookup, (uintptr_t)s.init_hess, "init_hess", test_name);
    ok &= check_settings_field(lookup, (uintptr_t)s.logger, "logger", test_name);
    ok &= check_settings_field(lookup, s.initialized, "initialized", test_name);
    ok &= check_settings_field(lookup, s.excess_dispersion, "excess_dispersion",
                               test_name);
    ok &= check_settings_field(lookup, s.verbose, "verbose", test_name);
    ok &= check_settings_keyword(lookup, s.inner_solver, "inner_solver", test_name);
    ok &= check_settings_keyword(lookup, s.prior_update, "prior_update", test_name);
    ok &= check_settings_keyword(lookup, s.acceptance, "acceptance", test_name);
    ok &= check_settings_field(lookup, (uintptr_t)s.context, "context", test_name);
  }

  return ok;
}

bool test_s_gek_settings_init(void) {
  bool ok = true;

  // call function and compare with the default settings without callback functions and
  // host context
  s_gek_settings_type s = s_gek_settings_init();
  if (!is_default_s_gek_settings(&s)) {
    fprintf(stderr, "test_s_gek_settings_init failed: Settings not initialized to the "
                    "default values without callback functions and host context.\n");
    ok = false;
  }

  return ok;
}
