// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/* Code the C integration tests of the extensions share: the field checks of their
 * layout tests, the warning of the orbital basis factories without a density matrix
 * evaluating function and the toy model, whose solves through the factories of the
 * orbital bases and density QN show that what crosses their C interfaces arrives,
 * together with what they share with the core's integration tests through
 * c_integration_tests.h. It calls no factory, since it is built with every extension.
 */

#ifndef COMMON_C_INTEGRATION_TESTS_H
#define COMMON_C_INTEGRATION_TESTS_H

#include <stdbool.h>

#include "opentrustregion.h"
#include "opentrustregion_common.h"

#include "../../../tests/c_integration_tests.h"

/* ------------------------------------------------------------------
 * Layout tests
 * ------------------------------------------------------------------ */

/* lookup of the reference value of a field of the settings of an extension by its
 * name, which the test reference module of every extension provides */
typedef void reference_field_fn(const char *name, c_real *value, char *keyword);

/* check that a field of settings filled with the reference values holds the reference
 * value of the field of that name, which is either a number or a keyword */
bool check_settings_field(reference_field_fn *lookup, c_real value, const char *name,
                          const char *test_name);
bool check_settings_keyword(reference_field_fn *lookup, const char *keyword,
                            const char *name, const char *test_name);

/* ------------------------------------------------------------------
 * Warning the factories of the orbital bases print without a density matrix evaluating
 * function, provided by common_test_reference.f90
 * ------------------------------------------------------------------ */

extern const char test_evaluate_dm_missing_warning[];

/* ------------------------------------------------------------------
 * Toy model
 * ------------------------------------------------------------------ */

/* largest number of AOs the arrays of the toy model hold */
#define TOY_MAX_N_AO 8
#define TOY_MAX_SIZE (TOY_MAX_N_AO * TOY_MAX_N_AO * 2)

/* the toy model of common_test_reference.f90: a mean-field energy in a non-orthogonal
 * AO basis with the core Hamiltonian h, the AO overlap matrix s, the two-electron
 * kernel g, whose Coulomb and exchange matrices are J[D]_pq = sum_rs g_pqrs D_rs and
 * K[D]_pq = sum_rs g_prqs D_rs, a non-linear term c sum_p n_p^3 of the total diagonal
 * density n_p and starting orbitals orthonormal in the AO metric, all column-major */
typedef struct {
  c_int n_ao;
  c_int n_occ[2];
  c_real h[TOY_MAX_N_AO * TOY_MAX_N_AO];
  c_real s[TOY_MAX_N_AO * TOY_MAX_N_AO];
  c_real g[TOY_MAX_N_AO * TOY_MAX_N_AO * TOY_MAX_N_AO * TOY_MAX_N_AO];
  c_real c;
  c_real mo_coeff[TOY_MAX_N_AO * TOY_MAX_N_AO];
} toy_model;

/* largest length of the messages the logging function of the toy model collects */
#define TOY_MAX_LOG 1024

/* host of the callback functions of the toy model, whose address is the host context
 * of the settings of an orbital basis or density QN, which counts their calls, records
 * whether the record of the state was current in them and collects the messages of the
 * logging function */
typedef struct {
  c_int n_particle;
  const solver_state_type *state;
  c_int n_evaluate_dm;
  c_int n_get_response;
  c_int n_project;
  c_int n_updates_seen;
  bool state_not_current;
  c_real response_dm[TOY_MAX_SIZE];
  char log[TOY_MAX_LOG];
} toy_host;

/* one solve of the toy model through the factory of an orbital basis and optionally
 * density QN, with the solver settings and state, the hosts of the callback functions
 * of both objects, the returned functions and objects and the orbitals the orbital
 * basis moves, MO coefficients or density matrices, with the ones it started from */
typedef struct {
  bool mo_basis;
  c_int n_particle;
  c_int n_mo;
  c_int n_param;
  solver_settings_type settings;
  solver_state_type state;
  toy_host host;
  toy_host density_qn_host;
  update_orbs_fp update_orbs;
  obj_func_fp obj_func;
  void *extension;
  void *density_qn_extension;
  c_real orbitals[TOY_MAX_SIZE];
  c_real start_orbitals[TOY_MAX_SIZE];
} toy_run;

/* the toy model, read from common_test_reference.f90 by toy_model_load, which returns
 * whether its dimensions fit into the arrays */
extern toy_model toy;
bool toy_model_load(const char *test_name);

/* prepare a solve of the toy model in the given orbital basis and shell, with solver
 * settings requesting a stability check and hosts reached through no other run, which
 * the factory of the orbital basis then fills, keeping the run in place since the hosts
 * point at its state */
void toy_run_init(toy_run *run, bool mo_basis, c_int n_particle);

/* the closed- or open-shell density matrix evaluating functions of the orbital bases
 * and density QN, which have to receive a host of the current run as their context */
c_int toy_evaluate_dm(const c_real *dm, c_real *energy, c_real *fock,
                      get_response_fp *get_response, void *context);
c_int toy_density_qn_evaluate_dm(const c_real *dm, c_real *energy, c_real *fock,
                                 c_real *v_coulomb, c_real *v_exchange,
                                 c_real *v_nonlinear, void *context);

/* the logging function of the settings of an orbital basis or density QN, which
 * appends the messages to the log of the host it receives as its context if it is one
 * of the current run */
void toy_logger(const char *message, void *context);

/* solve the prepared run and check that the solve converged, that the callback
 * functions received the hosts of the run as their context and the current record of
 * the state, that the exact Hessian of the orbital basis and, with density QN, its
 * density matrix evaluating function were used and that the orbitals moved and stayed
 * orthonormal */
bool toy_run_solve(toy_run *run, const char *test_name, const char *case_name);

#endif
