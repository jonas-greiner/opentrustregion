// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/* Pure-C integration tests for the public C interface declared in
 * opentrustregion_oao.h
 *
 * The Fortran-side oao_c_interface_unit_tests cover the bind(C) wrappers but never
 * compile against the C header itself. These tests do, and test_oao_settings_layout
 * checks that the settings struct has the size of the bind(C) type and that every
 * field is read back under its own name, so any drift between the bind(C) type
 * (Fortran) and the struct (C) is caught here, except a field left in the struct that
 * fits into the padding of the bind(C) type, which changes neither the size nor the
 * offset of a named field. The header's inline initialization, which
 * test_oao_settings_init checks, takes both sides from the same Fortran assignment, so
 * it cannot detect such drift. test_oao_factory_c solves the toy model of
 * common_c_integration_tests.h through the OAO factory the header declares, in the
 * closed and the open shell, and checks the warnings the factory prints. */

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "opentrustregion.h"
#include "opentrustregion_oao.h"

#include "oao_c_integration_tests.h"

/* ------------------------------------------------------------------
 * Reference settings provided by oao_test_reference.f90: the reference value of a field
 * by its name, C settings filled with these values field by field, the size of the
 * bind(C) settings type, a comparison with the default settings, and the warning the
 * factory prints for OAO settings that were not initialized
 * ------------------------------------------------------------------ */

void reference_oao_field(const char *name, c_real *value, char *keyword);
void get_reference_oao_values(oao_settings_type *settings);
size_t oao_settings_size(void);
bool is_default_oao_settings(const oao_settings_type *settings);
extern const char test_oao_settings_uninitialized_warning[];

/* ------------------------------------------------------------------
 * Helpers
 * ------------------------------------------------------------------ */

/* the setup declared in oao_c_integration_tests.h, which the integration tests of ARH
 * use too */
bool setup_oao_toy_run(toy_run *run, c_int n_particle, const char *test_name,
                       const char *case_name) {
  toy_run_init(run, false, n_particle);
  oao_settings_type settings = oao_settings_init();
  settings.context = &run->host;
  c_int error = oao_factory(run->orbitals, toy.s, n_particle, toy.n_ao, toy_evaluate_dm,
                            &run->update_orbs, &run->obj_func, &run->extension,
                            &run->settings, &run->state, &settings);
  if (error != 0) {
    fprintf(stderr, "%s failed: OAO factory produced error %d for the %s case.\n",
            test_name, (int)error, case_name);
    return false;
  }
  return true;
}

/* ------------------------------------------------------------------
 * Tests
 * ------------------------------------------------------------------ */

bool test_oao_settings_layout(void) {
  /* check that the settings struct has the size of the bind(C) settings type, so that
   * a field missing from either is detected, which the field checks below would not
   * read, and stop otherwise, since filling smaller settings would write past their
   * end */
  if (sizeof(oao_settings_type) != oao_settings_size()) {
    fprintf(stderr, "test_oao_settings_layout failed: Size of OAO settings differs "
                    "from the library's.\n");
    return false;
  }

  /* check that every field is read back under its own name */
  bool ok = true;
  oao_settings_type s = {0};
  get_reference_oao_values(&s);
  reference_field_fn *lookup = reference_oao_field;
  const char *test_name = "test_oao_settings_layout";
  ok &= check_settings_field(lookup, (uintptr_t)s.logger, "logger", test_name);
  ok &= check_settings_field(lookup, s.initialized, "initialized", test_name);
  ok &= check_settings_field(lookup, s.verbose, "verbose", test_name);
  ok &= check_settings_field(lookup, (uintptr_t)s.context, "context", test_name);

  return ok;
}

bool test_oao_settings_init(void) {
  bool ok = true;

  // call function and compare with the default settings without logging function and
  // host context
  oao_settings_type s = oao_settings_init();
  if (!is_default_oao_settings(&s)) {
    fprintf(stderr, "test_oao_settings_init failed: Settings not initialized to the "
                    "default values without logging function and host context.\n");
    ok = false;
  }

  return ok;
}

bool test_oao_factory_c(void) {
  bool ok = true;

  if (!toy_model_load("test_oao_factory_c"))
    return false;

  // solve the toy model in the closed and the open shell
  const char *shell_names[2] = {"closed-shell", "open-shell"};
  for (c_int n_particle = 1; n_particle <= 2; n_particle++) {
    toy_run run;
    if (!setup_oao_toy_run(&run, n_particle, "test_oao_factory_c",
                           shell_names[n_particle - 1])) {
      ok = false;
      continue;
    }
    ok &= toy_run_solve(&run, "test_oao_factory_c", shell_names[n_particle - 1]);
    release_extension(run.extension);
  }

  // check that the factory warns on the standard error about OAO settings and solver
  // settings that were not initialized, which provide no logging function, and only
  // through the logging function of the OAO settings about a missing density matrix
  // evaluating function
  const char *case_names[3] = {"for OAO settings that were not initialized",
                               "for solver settings that were not initialized",
                               "without a density matrix evaluating function"};
  const char *stderr_warnings[3] = {test_oao_settings_uninitialized_warning,
                                    test_settings_uninitialized_warning, NULL};
  for (int i_case = 0; i_case < 3; i_case++) {
    toy_run run;
    toy_run_init(&run, false, 1);
    oao_settings_type settings = oao_settings_init();
    settings.logger = toy_logger;
    settings.verbose = 3;
    settings.context = &run.host;
    if (i_case == 0)
      settings = (oao_settings_type){0};
    if (i_case == 1)
      run.settings = (solver_settings_type){0};
    start_stderr_capture();
    c_int error = oao_factory(run.orbitals, toy.s, 1, toy.n_ao,
                              i_case == 2 ? NULL : toy_evaluate_dm, &run.update_orbs,
                              &run.obj_func, &run.extension, &run.settings, &run.state,
                              &settings);
    if (!stop_stderr_capture(stderr_warnings[i_case])) {
      fprintf(stderr,
              "test_oao_factory_c failed: Warning on the standard error wrong %s.\n",
              case_names[i_case]);
      ok = false;
    }
    if ((strstr(run.host.log, test_evaluate_dm_missing_warning) != NULL) !=
        (i_case == 2)) {
      fprintf(stderr,
              "test_oao_factory_c failed: Warning about the density matrix evaluating "
              "function wrong %s.\n",
              case_names[i_case]);
      ok = false;
    }
    if (error != 0) {
      fprintf(stderr, "test_oao_factory_c failed: Produced error %s.\n",
              case_names[i_case]);
      ok = false;
      continue;
    }
    release_extension(run.extension);
  }

  return ok;
}
