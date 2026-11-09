// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/* Code the C integration tests of the core share with those of the extensions: the
 * capture of the standard error and the Hartmann 6D problem, which
 * c_integration_tests.c provides, and the warnings the library prints for settings
 * that were not initialized. */

#ifndef C_INTEGRATION_TESTS_H
#define C_INTEGRATION_TESTS_H

#include <stdbool.h>

#include "opentrustregion.h"

/* the Hartmann 6D problem, whose data opentrustregion_unit_tests.f90 provides: a
 * minimum and a point near it, and the function, its gradient and its Hessian at a
 * point, which c_integration_tests.c provides */
#define HARTMANN6D_N_PARAM 6
extern const c_real hartmann6d_minimum1[HARTMANN6D_N_PARAM];
extern const c_real hartmann6d_near_minimum[HARTMANN6D_N_PARAM];
c_real hartmann_func(const c_real x[HARTMANN6D_N_PARAM]);
void hartmann_grad(const c_real x[HARTMANN6D_N_PARAM], c_real grad[HARTMANN6D_N_PARAM]);
void hartmann_hess(const c_real x[HARTMANN6D_N_PARAM],
                   c_real hess[HARTMANN6D_N_PARAM][HARTMANN6D_N_PARAM]);

/* the warnings the library prints for settings and nested stability check settings
 * that were not initialized, provided by test_reference.f90 */
extern const char test_settings_uninitialized_warning[];
extern const char test_stability_settings_uninitialized_warning[];

/* run the calls between these two functions with the standard error redirected to a
 * temporary file, flushing what C and the Fortran runtime buffered before and after
 * them, the second returns whether the given warning was written to it, or with NULL
 * whether nothing was written to it */
void start_stderr_capture(void);
bool stop_stderr_capture(const char *warning);

#endif
