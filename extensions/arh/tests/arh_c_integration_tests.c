// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/* Pure-C integration tests for the public C interface declared in
 * opentrustregion_arh.h
 *
 * The Fortran-side arh_c_interface_unit_tests cover the bind(C) wrappers but never
 * compile against the C header itself. These tests do, and test_arh_settings_layout
 * checks that the settings struct has the size of the bind(C) type and that every
 * field is read back under its own name, so any drift between the bind(C) type
 * (Fortran) and the struct (C) is caught here, except a field left in the struct that
 * fits into the padding of the bind(C) type, which changes neither the size nor the
 * offset of a named field. The header's inline initialization, which
 * test_arh_settings_init checks, takes both sides from the same Fortran assignment, so
 * it cannot detect such drift. test_arh_factory_c solves the toy model of
 * common_c_integration_tests.h with ARH, built through the ARH factory the header
 * declares on the objects of both orbital bases in the closed and the open shell, and
 * checks the warning the factory prints. */

#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include "opentrustregion.h"
#include "opentrustregion_arh.h"

#include "../../mo/tests/mo_c_integration_tests.h"
#include "../../oao/tests/oao_c_integration_tests.h"

/* ------------------------------------------------------------------
 * Reference settings provided by arh_test_reference.f90: the reference value of a field
 * by its name, C settings filled with these values field by field, the size of the
 * bind(C) settings type, a comparison with the default settings, and the warning the
 * factory prints for ARH settings that were not initialized
 * ------------------------------------------------------------------ */

void reference_arh_field(const char *name, c_real *value, char *keyword);
void get_reference_arh_values(arh_settings_type *settings);
size_t arh_settings_size(void);
bool is_default_arh_settings(const arh_settings_type *settings);
extern const char test_arh_settings_uninitialized_warning[];

/* ------------------------------------------------------------------
 * Tests
 * ------------------------------------------------------------------ */

bool test_arh_settings_layout(void) {
  /* check that the settings struct has the size of the bind(C) settings type, so that
   * a field missing from either is detected, which the field checks below would not
   * read, and stop otherwise, since filling smaller settings would write past their
   * end */
  if (sizeof(arh_settings_type) != arh_settings_size()) {
    fprintf(stderr, "test_arh_settings_layout failed: Size of ARH settings differs "
                    "from the library's.\n");
    return false;
  }

  /* check that every field is read back under its own name */
  bool ok = true;
  arh_settings_type s = {0};
  get_reference_arh_values(&s);
  reference_field_fn *lookup = reference_arh_field;
  const char *test_name = "test_arh_settings_layout";
  ok &= check_settings_field(lookup, (uintptr_t)s.logger, "logger", test_name);
  ok &= check_settings_field(lookup, s.initialized, "initialized", test_name);
  ok &= check_settings_field(lookup, s.verbose, "verbose", test_name);
  ok &= check_settings_keyword(lookup, s.arh_type, "arh_type", test_name);
  ok &= check_settings_field(lookup, (uintptr_t)s.context, "context", test_name);

  return ok;
}

bool test_arh_settings_init(void) {
  bool ok = true;

  // call function and compare with the default settings without logging function and
  // host context
  arh_settings_type s = arh_settings_init();
  if (!is_default_arh_settings(&s)) {
    fprintf(stderr, "test_arh_settings_init failed: Settings not initialized to the "
                    "default values without logging function and host context.\n");
    ok = false;
  }

  return ok;
}

bool test_arh_factory_c(void) {
  bool ok = true;
  const char *test_name = "test_arh_factory_c";

  if (!toy_model_load(test_name))
    return false;

  // solve the toy model with ARH in both orbital bases and shells and compare with the
  // minimum the exact Hessian of the orbital basis reaches
  bool (*setup_basis_run[2])(toy_run *, c_int, const char *,
                             const char *) = {setup_mo_toy_run, setup_oao_toy_run};
  const char *basis_names[2] = {"MO-basis", "OAO-basis"};
  const char *shell_names[2] = {"closed-shell", "open-shell"};
  for (int i_basis = 0; i_basis < 2; i_basis++)
    for (c_int n_particle = 1; n_particle <= 2; n_particle++) {
      char case_name[64];
      snprintf(case_name, sizeof(case_name), "%s %s", basis_names[i_basis],
               shell_names[n_particle - 1]);

      // solve with the exact Hessian of the orbital basis
      toy_run run;
      if (!setup_basis_run[i_basis](&run, n_particle, test_name, case_name)) {
        ok = false;
        continue;
      }
      bool exact_ok = toy_run_solve(&run, test_name, case_name);
      c_real exact_energy = run.state.func;
      release_extension(run.extension);
      if (!exact_ok) {
        ok = false;
        continue;
      }

      // wire the settings with a new orbital basis object and set the exact Hessian
      // linear transformation its orbital update returns at the starting orbitals for
      // the stability check, before the ARH factory wires them with ARH
      if (!setup_basis_run[i_basis](&run, n_particle, test_name, case_name)) {
        ok = false;
        continue;
      }
      c_real kappa[TOY_MAX_SIZE] = {0}, func, grad[TOY_MAX_SIZE], h_diag[TOY_MAX_SIZE];
      hess_x_fp hess_x;
      c_int error =
          run.update_orbs(kappa, &func, grad, h_diag, &hess_x, run.settings.context);
      if (error == 0) {
        run.settings.stability_hess_x = hess_x;
        arh_settings_type arh_settings = arh_settings_init();
        arh_settings.context = &run.arh_host;
        error =
            arh_factory(toy_arh_evaluate_dm, &run.update_orbs, &run.obj_func,
                        &run.arh_extension, &run.settings, &run.state, &arh_settings);
      }
      if (error != 0) {
        fprintf(stderr, "%s failed: Setup of ARH produced error %d for the %s case.\n",
                test_name, (int)error, case_name);
        ok = false;
        release_extension(run.extension);
        continue;
      }

      // solve with ARH, which reaches the same minimum
      ok &= toy_run_solve(&run, test_name, case_name);
      if (fabs(run.state.func - exact_energy) > 1e-8) {
        fprintf(stderr,
                "%s failed: Minimum differs from the one of the exact Hessian "
                "for the %s case.\n",
                test_name, case_name);
        ok = false;
      }
      release_extension(run.arh_extension);
      release_extension(run.extension);
    }

  // check that the factory warns on the standard error about ARH settings that were
  // not initialized, which provide no logging function
  toy_run run;
  if (!setup_mo_toy_run(&run, 1, test_name, "closed-shell"))
    return false;
  arh_settings_type arh_settings = {0};
  start_stderr_capture();
  c_int error =
      arh_factory(toy_arh_evaluate_dm, &run.update_orbs, &run.obj_func,
                  &run.arh_extension, &run.settings, &run.state, &arh_settings);
  if (!stop_stderr_capture(test_arh_settings_uninitialized_warning)) {
    fprintf(stderr,
            "%s failed: Warning on the standard error wrong for ARH settings "
            "that were not initialized.\n",
            test_name);
    ok = false;
  }
  if (error != 0) {
    fprintf(stderr,
            "%s failed: Produced error for ARH settings that were not "
            "initialized.\n",
            test_name);
    ok = false;
  } else {
    release_extension(run.arh_extension);
  }
  release_extension(run.extension);

  return ok;
}
