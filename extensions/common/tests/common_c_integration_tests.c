// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/* Code the C integration tests of the extensions share, see
 * common_c_integration_tests.h */

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "common_c_integration_tests.h"

/* ------------------------------------------------------------------
 * Layout tests
 * ------------------------------------------------------------------ */

bool check_settings_field(reference_field_fn *lookup, c_real value, const char *name,
                          const char *test_name) {
  c_real ref_value;
  char ref_keyword[OTR_KW_LEN + 1];
  lookup(name, &ref_value, ref_keyword);
  if (value != ref_value) {
    fprintf(stderr, "%s failed: Field %s misplaced.\n", test_name, name);
    return false;
  }
  return true;
}

bool check_settings_keyword(reference_field_fn *lookup, const char *keyword,
                            const char *name, const char *test_name) {
  c_real ref_value;
  char ref_keyword[OTR_KW_LEN + 1];
  lookup(name, &ref_value, ref_keyword);
  if (strcmp(keyword, ref_keyword) != 0) {
    fprintf(stderr, "%s failed: Field %s misplaced.\n", test_name, name);
    return false;
  }
  return true;
}

/* ------------------------------------------------------------------
 * Toy model
 * ------------------------------------------------------------------ */

/* the toy model provided by common_test_reference.f90 */
extern const c_int test_n_ao;
void get_toy_model(c_real *h, c_real *s, c_real *g, c_real *nonlinear_factor,
                   c_real *mo_coeff, c_int *n_occ);

toy_model toy;

bool toy_model_load(const char *test_name) {
  if (test_n_ao > TOY_MAX_N_AO) {
    fprintf(stderr, "%s failed: Toy model too large for its arrays.\n", test_name);
    return false;
  }
  toy.n_ao = test_n_ao;
  get_toy_model(toy.h, toy.s, toy.g, &toy.c, toy.mo_coeff, toy.n_occ);
  return true;
}

/* the hosts of the current run, the only contexts the callback functions accept, and
 * the number of calls that received another context */
static toy_host *current_hosts[2];
static int n_wrong_context = 0;

static toy_host *host_of(void *context) {
  for (int i = 0; i < 2; i++)
    if (context == current_hosts[i])
      return current_hosts[i];
  n_wrong_context++;
  return NULL;
}

/* record whether the number of orbital updates and objective function evaluations in
 * the record of the state never decreased between the calls of the callback functions
 * of a host, which it does once the record is not kept current */
static void note_state(toy_host *host) {
  c_int n_updates = host->state->n_update_orbs + host->state->n_obj_func;
  if (n_updates < host->n_updates_seen)
    host->state_not_current = true;
  if (n_updates > host->n_updates_seen)
    host->n_updates_seen = n_updates;
}

/* the Coulomb and the exchange matrix of a density matrix of one channel */
static void coulomb(const c_real *dm, c_real *j) {
  c_int n = toy.n_ao;
  for (c_int p = 0; p < n; p++)
    for (c_int q = 0; q < n; q++) {
      c_real sum = 0.0;
      for (c_int r = 0; r < n; r++)
        for (c_int s = 0; s < n; s++)
          sum += toy.g[p + n * (q + n * (r + n * s))] * dm[r + n * s];
      j[p + n * q] = sum;
    }
}

static void exchange(const c_real *dm, c_real *k) {
  c_int n = toy.n_ao;
  for (c_int p = 0; p < n; p++)
    for (c_int q = 0; q < n; q++) {
      c_real sum = 0.0;
      for (c_int r = 0; r < n; r++)
        for (c_int s = 0; s < n; s++)
          sum += toy.g[p + n * (r + n * (q + n * s))] * dm[r + n * s];
      k[p + n * q] = sum;
    }
}

/* the total diagonal density of the closed-shell density matrix of one spin channel,
 * counted twice, or of the open-shell density matrices of both channels */
static void total_diagonal(c_int n_particle, const c_real *dm, c_real *n_diag) {
  c_int n = toy.n_ao;
  for (c_int p = 0; p < n; p++)
    n_diag[p] =
        n_particle == 1 ? 2.0 * dm[p + n * p] : dm[p + n * p] + dm[n * n + p + n * p];
}

/* the energy and, where wanted, the Fock matrix and the Coulomb, exact-exchange and
 * non-linear potentials of the closed-shell density matrix of one spin channel or the
 * open-shell density matrices of both channels */
static c_real toy_outputs(c_int n_particle, const c_real *dm, c_real *fock,
                          c_real *v_coulomb, c_real *v_exchange, c_real *v_nonlinear) {
  c_int n = toy.n_ao, nn = n * n;
  c_real n_diag[TOY_MAX_N_AO], j[TOY_MAX_SIZE], k[TOY_MAX_SIZE],
      j_total[TOY_MAX_N_AO * TOY_MAX_N_AO], dm_total[TOY_MAX_N_AO * TOY_MAX_N_AO];
  total_diagonal(n_particle, dm, n_diag);
  for (c_int i = 0; i < n_particle; i++) {
    coulomb(dm + i * nn, j + i * nn);
    exchange(dm + i * nn, k + i * nn);
  }
  for (c_int pq = 0; pq < nn; pq++)
    dm_total[pq] = n_particle == 1 ? 2.0 * dm[pq] : dm[pq] + dm[nn + pq];
  coulomb(dm_total, j_total);

  /* closed shell: E = 2 tr(h D) + tr(D J[2 D]) - tr(D K[D]) + c sum_p n_p^3, open
   * shell: E = sum_s tr(h D_s) + tr(D_s J[D_a + D_b]) / 2 - tr(D_s K[D_s]) / 2 + c
   * sum_p n_p^3, whose Fock matrices are h + J[D_a + D_b] - K[D_s] + 3 c n_p^2 */
  c_real weight = n_particle == 1 ? 2.0 : 1.0, energy = 0.0;
  for (c_int i = 0; i < n_particle; i++)
    for (c_int pq = 0; pq < nn; pq++)
      energy += weight * dm[i * nn + pq] *
                (toy.h[pq] + 0.5 * j_total[pq] - 0.5 * k[i * nn + pq]);
  for (c_int p = 0; p < n; p++)
    energy += toy.c * n_diag[p] * n_diag[p] * n_diag[p];
  for (c_int i = 0; i < n_particle; i++)
    for (c_int q = 0; q < n; q++)
      for (c_int p = 0; p < n; p++) {
        c_int pq = p + n * q, ipq = i * nn + pq;
        c_real nonlinear = p == q ? 3.0 * toy.c * n_diag[p] * n_diag[p] : 0.0;
        if (fock)
          fock[ipq] = toy.h[pq] + j_total[pq] - k[ipq] + nonlinear;
        if (v_coulomb)
          v_coulomb[ipq] = n_particle == 1 ? j_total[pq] : j[ipq];
        if (v_exchange)
          v_exchange[ipq] = k[ipq];
        if (v_nonlinear)
          v_nonlinear[ipq] = nonlinear;
      }
  return energy;
}

/* the response of the Fock matrices to a density matrix displacement at the density
 * matrices the density matrix evaluating function returned it for */
static c_int toy_get_response(const c_real *ddm, c_real *response, void *context) {
  toy_host *host = host_of(context);
  if (!host)
    return 1;
  host->n_get_response++;
  c_int n = toy.n_ao, nn = n * n, n_particle = host->n_particle;
  c_real n_diag[TOY_MAX_N_AO], dn_diag[TOY_MAX_N_AO], k[TOY_MAX_SIZE],
      j_total[TOY_MAX_N_AO * TOY_MAX_N_AO], ddm_total[TOY_MAX_N_AO * TOY_MAX_N_AO];
  total_diagonal(n_particle, host->response_dm, n_diag);
  total_diagonal(n_particle, ddm, dn_diag);
  for (c_int pq = 0; pq < nn; pq++)
    ddm_total[pq] = n_particle == 1 ? 2.0 * ddm[pq] : ddm[pq] + ddm[nn + pq];
  coulomb(ddm_total, j_total);
  for (c_int i = 0; i < n_particle; i++) {
    exchange(ddm + i * nn, k + i * nn);
    for (c_int q = 0; q < n; q++)
      for (c_int p = 0; p < n; p++) {
        c_int pq = p + n * q, ipq = i * nn + pq;
        response[ipq] = j_total[pq] - k[ipq] +
                        (p == q ? 6.0 * toy.c * n_diag[p] * dn_diag[p] : 0.0);
      }
  }
  return 0;
}

c_int toy_evaluate_dm(const c_real *dm, c_real *energy, c_real *fock,
                      get_response_fp *get_response, void *context) {
  toy_host *host = host_of(context);
  if (!host)
    return 1;
  host->n_evaluate_dm++;
  note_state(host);
  *energy = toy_outputs(host->n_particle, dm, fock, NULL, NULL, NULL);

  /* the response belongs to the density matrices it was requested for, not to those of
   * a later evaluation without one */
  if (get_response) {
    memcpy(host->response_dm, dm,
           (size_t)(host->n_particle * toy.n_ao * toy.n_ao) * sizeof(c_real));
    *get_response = toy_get_response;
  }
  return 0;
}

c_int toy_density_qn_evaluate_dm(const c_real *dm, c_real *energy, c_real *fock,
                                 c_real *v_coulomb, c_real *v_exchange,
                                 c_real *v_nonlinear, void *context) {
  toy_host *host = host_of(context);
  if (!host)
    return 1;
  host->n_evaluate_dm++;
  note_state(host);
  *energy = toy_outputs(host->n_particle, dm, fock, v_coulomb, v_exchange, v_nonlinear);
  return 0;
}

void toy_logger(const char *message, void *context) {
  toy_host *host = host_of(context);
  if (host)
    strncat(host->log, message, sizeof(host->log) - strlen(host->log) - 1);
}

/* ------------------------------------------------------------------
 * Runs
 * ------------------------------------------------------------------ */

static void init_host(toy_host *host, c_int n_particle,
                      const solver_state_type *state) {
  memset(host, 0, sizeof(*host));
  host->n_particle = n_particle;
  host->state = state;
}

void toy_run_init(toy_run *run, bool mo_basis, c_int n_particle) {
  c_int n = toy.n_ao;
  memset(run, 0, sizeof(*run));
  run->mo_basis = mo_basis;
  run->n_particle = n_particle;
  run->settings = solver_settings_init();
  run->settings.stability = true;
  init_host(&run->host, n_particle, &run->state);
  init_host(&run->density_qn_host, n_particle, &run->state);
  current_hosts[0] = &run->host;
  current_hosts[1] = &run->density_qn_host;
  n_wrong_context = 0;

  /* MO coefficients of fewer MOs than AOs for every channel, so that the numbers of
   * AOs and MOs cannot be confused, or the density matrices of their occupied orbitals
   * for every channel */
  if (mo_basis) {
    run->n_mo = n - 1;
    for (c_int i = 0; i < n_particle; i++) {
      for (c_int pq = 0; pq < n * run->n_mo; pq++)
        run->orbitals[i * n * run->n_mo + pq] = toy.mo_coeff[pq];
      run->n_param += toy.n_occ[i] * (run->n_mo - toy.n_occ[i]);
    }
  } else {
    for (c_int i = 0; i < n_particle; i++)
      for (c_int q = 0; q < n; q++)
        for (c_int p = 0; p < n; p++) {
          c_real sum = 0.0;
          for (c_int o = 0; o < toy.n_occ[i]; o++)
            sum += toy.mo_coeff[p + n * o] * toy.mo_coeff[q + n * o];
          run->orbitals[i * n * n + p + n * q] = sum;
        }
    run->n_param = n_particle * n * (n - 1) / 2;
  }
  memcpy(run->start_orbitals, run->orbitals, sizeof(run->orbitals));
}

/* check that the callback functions of a host were called with the record of the
 * state, which never fell behind and had moved on by the solve */
static bool check_host(const toy_host *host, const toy_run *run, const char *name,
                       const char *test_name, const char *case_name) {
  bool ok = true;
  if (host->n_evaluate_dm == 0) {
    fprintf(stderr,
            "%s failed: Density matrix evaluating function of the %s not "
            "called for the %s case.\n",
            test_name, name, case_name);
    ok = false;
  }
  if (host->state_not_current || host->n_updates_seen == 0 ||
      host->n_updates_seen > run->state.n_update_orbs + run->state.n_obj_func) {
    fprintf(stderr,
            "%s failed: Record of the state not current in the callback "
            "functions of the %s for the %s case.\n",
            test_name, name, case_name);
    ok = false;
  }
  return ok;
}

/* the largest deviation of the overlap of the orbitals, or of the density matrices
 * multiplied by the overlap matrix from both sides, from what orthonormal orbitals
 * give */
static c_real orthonormality_error(const toy_run *run) {
  c_int n = toy.n_ao;
  c_real error = 0.0;
  for (c_int i = 0; i < run->n_particle; i++) {
    c_int m = run->mo_basis ? run->n_mo : n;
    const c_real *x = run->orbitals + i * n * m;
    for (c_int a = 0; a < m; a++)
      for (c_int b = 0; b < m; b++) {
        /* MO basis: C^T S C = I, OAO basis: D S D = D */
        c_real sum = 0.0;
        for (c_int p = 0; p < n; p++)
          for (c_int q = 0; q < n; q++)
            sum += x[p + n * a] * toy.s[p + n * q] * x[q + n * b];
        c_real target = run->mo_basis ? (a == b ? 1.0 : 0.0) : x[a + n * b];
        error = fmax(error, fabs(sum - target));
      }
  }
  return error;
}

bool toy_run_solve(toy_run *run, const char *test_name, const char *case_name) {
  bool ok = true;
  c_int error = solver(run->update_orbs, run->obj_func, run->n_param, &run->settings,
                       &run->state);
  if (error != 0) {
    fprintf(stderr, "%s failed: Solver produced error %d for the %s case.\n", test_name,
            (int)error, case_name);
    return false;
  }
  if (!(run->state.grad_rms < run->settings.conv_tol ||
        run->state.max_precision_reached)) {
    fprintf(stderr, "%s failed: Solver did not converge for the %s case.\n", test_name,
            case_name);
    ok = false;
  }
  if (n_wrong_context != 0) {
    fprintf(stderr,
            "%s failed: Callback functions received a wrong context for the "
            "%s case.\n",
            test_name, case_name);
    ok = false;
  }

  /* the exact Hessian of the orbital basis is used for the solve or, with density QN,
   * for the stability check, so the response function is called */
  ok &= check_host(&run->host, run, "orbital basis", test_name, case_name);
  if (run->host.n_get_response == 0) {
    fprintf(stderr,
            "%s failed: Response function of the orbital basis not called for "
            "the %s case.\n",
            test_name, case_name);
    ok = false;
  }
  if (run->density_qn_extension)
    ok &= check_host(&run->density_qn_host, run, "density QN object", test_name,
                     case_name);

  /* the orbitals of the caller moved and stayed orthonormal */
  c_real moved = 0.0;
  for (c_int i = 0; i < TOY_MAX_SIZE; i++)
    moved = fmax(moved, fabs(run->orbitals[i] - run->start_orbitals[i]));
  if (moved < 1e-3) {
    fprintf(stderr, "%s failed: Orbitals of the caller not moved for the %s case.\n",
            test_name, case_name);
    ok = false;
  }
  if (orthonormality_error(run) > 1e-8) {
    fprintf(stderr,
            "%s failed: Orbitals of the caller not orthonormal for the %s "
            "case.\n",
            test_name, case_name);
    ok = false;
  }
  return ok;
}
