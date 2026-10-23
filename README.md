![License](https://img.shields.io/github/license/eriksen-lab/opentrustregion)
![CI](https://github.com/eriksen-lab/opentrustregion/actions/workflows/main.yml/badge.svg)
[![codecov](https://codecov.io/github/eriksen-lab/opentrustregion/graph/badge.svg?token=NJIM6FDADD)](https://codecov.io/github/eriksen-lab/opentrustregion)
[![DOI](https://zenodo.org/badge/855894860.svg)](https://doi.org/10.5281/zenodo.17142554)

# OpenTrustRegion: A Reusable Library for Second-Order Trust Region Orbital Optimization

This library provides a robust and flexible implementation for second-order trust region orbital optimization, with extensive customization options to suit various use cases.

The following paper documents the theory and implementation of the methodology in OpenTrustRegion, and should be cited in any work using OpenTrustRegion:  

- Greiner, J.; Høyvik, I.-M.; Lehtola, S.; Eriksen, J. J. 
  A Reusable Library for Second-Order Orbital Optimization Using the Trust Region Method. Journal of Chemical Theory and Computation 2026, 22(2), 881–895. 
  DOI: [10.1021/acs.jctc.5c01576](https://doi.org/10.1021/acs.jctc.5c01576). 
  arXiv: [2509.13931](https://arxiv.org/abs/2509.13931).

## Installation

The library needs a Fortran 2008 compiler, a C compiler, CMake 3.24 or newer and BLAS and LAPACK. It is tested with gfortran 10 and newer, Intel ifx and LLVM Flang 21 and newer, and against OpenBLAS, the reference BLAS and LAPACK, Intel MKL and Apple Accelerate.

### Fortran or C Usage
To build the library for Fortran or C:

```sh
mkdir build
cd build
cmake ..
cmake --build .
```

The installation can be tested by running the ```testsuite.py``` file in the ```pyopentrustregion``` directory.

A Fortran or C program has to be rebuilt against the modules or the header of the library version it links, since the settings and state types can change from one version to the next.

### Python Usage
To install the library for use with Python:

```sh
pip install .
```

The installation can be tested by running

```sh
python3 -m pyopentrustregion.testsuite
```

### CMake Configuration Options

The build process can be customized using the following CMake options:

| Option | Type | Default | Description |
|--------|------|----------|------------|
| **BUILD_SHARED_LIBS** | `BOOL` | `OFF` | Build shared libraries (`.so`, `.dylib`) instead of static ones. |
| **OpenTrustRegion_BUILD_TESTING** | `BOOL` | `ON` when built as the top-level project, `OFF` otherwise | Build the project’s testsuite. |
| **OpenTrustRegion_INSTALL_TESTING_LIBRARY** | `BOOL` | `OFF` | Install the library the testsuite runs (`libotrtestsuite`) together with the library. |
| **OpenTrustRegion_INSTALL_CMAKEDIR** | `STRING` | (auto) | Project install directory. |
| **CMAKE_BUILD_TYPE** | `STRING` | `Release` | Choose the build type (`Debug`, `Release`, etc.). |
| **INTEGER_SIZE** | `STRING` | *(auto)* | Set the integer precision to `4` (32-bit) or `8` (64-bit). Required when providing custom BLAS/LAPACK libraries. Otherwise defaults to 32-bit integers and tries to locate compatible BLAS and LAPACK libraries. Falls back to 64-bit integers if 32-bit libraries cannot be found. The resulting library name reflects the chosen integer precision (`libopentrustregion_32.*` or `libopentrustregion_64.*`). C programs compiled against the header of a library with 64-bit integers need `OTR_ILP64` defined, which linking the CMake target `OpenTrustRegion::opentrustregion` does automatically. |
| **BLAS_LIBRARIES** | `PATH` | *(auto)* | Path(s) to BLAS libraries. If not provided, CMake attempts to locate a suitable BLAS automatically. |
| **LAPACK_LIBRARIES** | `PATH` | *(auto)* | Path(s) to LAPACK libraries. If not provided, CMake attempts to locate a suitable LAPACK automatically. |
| **OpenTrustRegion_HOST_PROVIDES_BLAS** | `BOOL` | `OFF` | When enabled, OpenTrustRegion will not attempt to detect or link BLAS/LAPACK and the testsuite is automatically disabled. The calling program must provide BLAS/LAPACK routines that expose the unsuffixed symbol names (for example `ddot`, `dsyev`) with an integer width matching `INTEGER_SIZE`; otherwise linking will fail loudly. |
| **OpenTrustRegion_ENABLE_XHOST** | `BOOL` | `ON` | Optimize the Release build for the current machine's instruction set (`-march=native` for GNU, `-xHost` for Intel). Automatically disabled when cross-compiling, regardless of this setting. Turn off when building for a different machine than the one compiling (e.g. packaging/conda builds). |

### Nested and Concurrent Calls

The library keeps no state of its own between calls, so a solver or stability check call can be nested inside a callback function of another call, or run concurrently with another call in a different thread, as long as each call is given its own state record and, where an extension factory wired the solver settings, its own settings object, since these own the extension object a call works with. The same holds for calculations with extensions, whose factories create a new object for every call, see [Optional Extensions](#optional-extensions), as long as each calculation has its own solver settings, state record and orbitals. A nested call re-enters procedures of the library that are still running, which Fortran only permits for procedures compiled for recursion, and a concurrent call needs their local variables on its own thread's stack. The library then has to be built with `-frecursive` for gfortran (implied by `-fopenmp`) or `-recursive` for Intel:

```sh
cmake .. -DCMAKE_Fortran_FLAGS=-frecursive                  # Fortran or C
CMAKE_FLAGS="-DCMAKE_Fortran_FLAGS=-frecursive" pip install .  # Python
```

Every call also reseeds the random number generator of the Fortran runtime, which the whole program shares, so a nested or concurrent call changes the random trial vectors, and with them the iterations, of the other call, though not its correctness.

## Program Interfaces

The PySCF interface is available as an extension hosted at https://github.com/eriksen-lab/pyscf_opentrustregion. To install it, simply add its path to the **`PYSCF_EXT_PATH`** environment variable:
```sh
export PYSCF_EXT_PATH=path/to/pyscf_opentrustregion
```
Usage examples can be found in the **`examples`** directory of the PySCF interface repository.

The interface supports Hartree–Fock and DFT calculations via the **`mf_to_otr`** function, which wraps PySCF **`HF`** and **`KS`** objects into their OpenTrustRegion counterparts. Similarly, localization methods are available through the **`BoysOTR`**, **`PipekMezeyOTR`**, and **`EdmistonRuedenbergOTR`** classes, and state-specific CASSCF calculations are supported via the **`casscf_to_otr`** function applied to a PySCF **`CASSCF`** object. All returned objects are fully compatible with the original PySCF classes and can be used interchangeably.

Optional settings can be adjusted by modifying object attributes directly. Orbital optimization and internal stability analysis are performed using the **`kernel`** and **`stability_check`** member functions, respectively.

## Usage

The optimization process is initiated by calling a `solver` subroutine. This routine takes the following arguments:

### Arguments

- **`update_orbs`** (subroutine):  
  Accepts and applies a variable update (e.g., orbital rotation), updates the host's internal state, and provides:
  - Objective function value (real)
  - Gradient (real array, written in-place)
  - Hessian diagonal (real array, written in-place)
  - A **`hess_x`** subroutine that performs Hessian-vector products:
    - Accepts a trial vector and writes the result of the Hessian transformation into an output array (real array, written in-place)
    - Returns an integer error code (0 for success, positive integers < 100 for errors)
    - Receives the host context as its last argument.
    - Has to be provided whenever the orbital update succeeds, otherwise the solver fails. In Fortran, the procedure pointer argument is therefore `intent(inout)` rather than `intent(out)`, since the solver disassociates it before every call to detect a missing one.
  - Returns an integer error code (0 for success, positive integers < 100 for errors)
  - Receives the host context as its last argument.
- **`obj_func`** (function):  
  Accepts and applies a variable update (e.g., orbital rotation) and returns:
  - Objective function value (real)
  - An integer error code (0 for success, positive integers < 100 for errors)
  - Receives the host context as its last argument.
- **`n_param`** (integer): Specifies the number of parameters to be optimized.
- **`error`** (integer): An integer code indicating the success or failure of the solver. The error code structure is explained below.
- **`settings`** (solver_settings_type, optional): Settings object which controls optional arguments as described below. The solver only reads it. Without it, the default settings are used.
- **`state`** (solver_state_type, optional): Record that the solver keeps up to date during the optimization and that holds the final state of the optimization once the solver returns, as described below.

---

The following Fortran snippet demonstrates how to use the `solver` interface:

```fortran
use opentrustregion, only: ip, rp, update_orbs_type, obj_func_type, &
                           solver_settings_type, solver_state_type, solver

procedure(update_orbs_type), pointer :: update_orbs_funptr
procedure(obj_func_type), pointer :: obj_func_funptr
integer(ip) :: n_param, error
type(solver_settings_type) :: settings
type(solver_state_type) :: state

! set callback function pointers to existing implementations
update_orbs_funptr => update_orbs
obj_func_funptr => obj_func

! override default settings
settings%conv_tol = 1e-6_rp
settings%n_macro = 100
settings%subsystem_solver = "tcg"

! hand the callback functions whatever the host needs to reach its own data
settings%context => host_data

! run solver
call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings, state)

! read the final state
print *, "Number of orbital updates:", state%n_update_orbs
```

- Callback function pointers (`update_orbs_funptr`, `obj_func_funptr`) point to existing implementations elsewhere in the program.
- `n_param` is also assumed to be defined elsewhere.
- A declared settings object holds the default settings, which can be overridden (here, `conv_tol`, `n_macro` and `subsystem_solver`); assigning `solver_settings_type()` resets an object to them.
- `host_data` is a variable of any type declared with the `target` attribute. Every callback function receives it as its last argument, declared as `class(*), intent(in), pointer :: context`, and recovers its own type with `select type`. The `intent(in)` applies to the pointer rather than the data, so a callback function may modify `host_data` but not point `context` elsewhere. Leaving `settings%context` unset is fine; the callback functions then receive an unassociated pointer.
- Finally, the `solver` is called with the callback functions, the settings and a state record. Both of these are optional: without settings, the default settings are used, and the state is only needed to read it.
- After the call, `state` holds the final state of the optimization (here, the number of orbital updates `n_update_orbs`), see [Output](#output). The settings are left unchanged.

---

The following C snippet demonstrates the equivalent usage through the C interface:

```c
#include <stdio.h>
#include <string.h>
#include "opentrustregion.h"

c_int n_param;

// set callback function pointers to existing implementations
update_orbs_fp update_orbs_funptr = update_orbs;
obj_func_fp obj_func_funptr = obj_func;

// initialize settings
solver_settings_type settings = solver_settings_init();
solver_state_type state;

// override default settings
settings.conv_tol = 1e-6;
settings.n_macro = 100;
strcpy(settings.subsystem_solver, "tcg");

// hand the callback functions whatever the host needs to reach its own data
settings.context = &host_data;

// run solver
c_int error = solver(update_orbs_funptr, obj_func_funptr, n_param, &settings, &state);

// read the final state
printf("Number of orbital updates: %lld\n", (long long)state.n_update_orbs);
```

- Callback function pointers (`update_orbs_funptr`, `obj_func_funptr`) point to existing implementations elsewhere in the program.
- `n_param` is also assumed to be defined elsewhere.
- Solver settings are initialized via a small helper function `solver_settings_init()`, which returns a struct with default values. Individual settings (here, `conv_tol`, `n_macro` and `subsystem_solver`) can then be overridden.
- `host_data` is any host object. Every callback function receives `&host_data` as its `void *context` argument and casts it back. Leaving `settings.context` as `NULL` is fine; the callback functions then receive `NULL`.
- Finally, the `solver` is called with the callback functions and pointers to the settings and to a state struct, and directly returns an error code in typical C fashion. Either pointer may be `NULL`: without settings, the default settings are used, and without a state struct, none is returned. Since the host context is a setting, a host that needs one passes settings.
- After the call, `state` holds the final state of the optimization (here, `n_update_orbs`), also after an error, see [Output](#output). The state struct needs no initialization, and the settings are left unchanged.

---

The following Python snippet demonstrates the equivalent usage through the Python interface:

```python
from pyopentrustregion import SolverSettings, SolverState, solver

# initialize settings and state
settings = SolverSettings()
state = SolverState()

# override default settings
settings.conv_tol = 1e-6
settings.n_macro = 100
settings.subsystem_solver = "tcg"

# run solver
solver(update_orbs, obj_func, n_param, settings, state)

# read the final state
print(f"Number of orbital updates: {state.n_update_orbs}")
```

- Callback functions (`update_orbs`, `obj_func`) are defined elsewhere in the program.
- `n_param` is also assumed to be defined elsewhere.
- Solver settings are initialized via the `SolverSettings` class, which returns an object with default values; individual settings (here, `conv_tol`, `n_macro` and `subsystem_solver`) can then be overridden.
- There is no `context` setting in Python. Any callable works as a callback function, so a closure, a bound method or a `functools.partial` already carries whatever host data the callback function needs.
- Finally, the `solver` is called with the callback functions, the settings and a state object, and errors can be caught in pythonic fashion in the form of a `RuntimeError`. The settings and the state are optional: without settings, the default settings are used, and the state is only needed to read it.
- After the call, `state` holds the final state of the optimization (here, `n_update_orbs`), see [Output](#output), whose fields are read-only attributes. The settings are left unchanged.

### Optional Settings
The optimization process can be fine-tuned using the following settings:

- **`precond`** (subroutine): Applies a preconditioner to a residual vector, given a level shift `mu`. Writes the result in-place into a provided array and returns an integer error code (0 for success, positive integers < 100 for errors). The Davidson iterations of the Davidson-family subsystem solvers use it to approximate `(H - mu * I)^-1`, also at `mu = 0`, while their Jacobi-Davidson iterations do not call it; `mu` legitimately ranges over any real value, so `precond` is not guaranteed to be positive definite (e.g. near a saddle point). The `"tcg"` and `"gltr"` subsystem solvers use `precond_pd` instead. Receives the host context as its last argument.
- **`project`** (subroutine): Applies a projection in-place to a provided vector and returns an integer error code (0 for success, positive integers < 100 for errors). Required for optimization using non-redundant parameters. When this is used, all other passed routines (`update_orbs`, `hess_x`, `precond`, `precond_pd` and `stability_hess_x`) must be self-projecting. Receives the host context as its last argument.
- **`precond_pd`** (subroutine): Applies a positive-definite preconditioner to a residual vector. Writes the result in-place into a provided array and returns an integer error code (0 for success, positive integers < 100 for errors). Used by the `"tcg"` and `"gltr"` subsystem solvers to define the ellipsoidal trust-region metric (see `trust_region_shape` below); unlike `precond`, this callback takes no level shift and must always return a positive-definite result. Receives the host context as its last argument.
- **`get_extra_trial_vectors`** (subroutine): Returns extra trial vectors that seed the leading part of the initial trial space of the Davidson-family subsystem solvers, written in-place to the provided matrix, whose column count is `n_extra_trial_vectors` and is passed to the callback itself by the C and Python interfaces, since those only receive a flat pointer. In Python, the matrix holds one trial vector per row instead. The vectors do not need to be orthonormalized or projected, and a vanishing column is dropped, which is how the callback signals that it has no vector to contribute for that slot. Additionally, outputs an integer code indicating the success or failure of the function, positive integers less than 100 represent error conditions. Receives the host context as its last argument.
- **`conv_check`** (function): Returns whether the optimization has converged due to some supplied convergence criterion, which replaces the criterion on the RMS gradient (`conv_tol`). A point the solver reaches by a step along an unstable mode does not count as converged, whatever the check returns. Additionally, outputs an integer code indicating the success or failure of the function, positive integers less than 100 represent error conditions. Receives the host context as its last argument.
- **`stability_hess_x`** (subroutine): Applies a different Hessian linear transformation to a trial vector for the stability check than the one `update_orbs` returns for the optimization, and returns an integer error code (0 for success, positive integers < 100 for errors). Intended for optimization with an approximate Hessian, where the stability of a converged solution has to be decided with the exact Hessian. Receives the host context as its last argument.
- **`stability`** (boolean): Determines whether a stability check is performed upon convergence. The solver checks the stability of a starting point that is already converged and of a point with an exactly vanishing gradient in any case. From an unstable point, it steps along the unstable mode until the gradient RMS of the quadratic model is ten times the larger of the current one and `conv_tol`, or by the trust radius if both vanish or the step would be longer, and does not take the point it reaches as converged.
- **`hess_symm`** (boolean): Determines whether the supplied Hessian is symmetric. This is sometimes not the case for approximate Hessians. A non-symmetric Hessian is only supported by the `"davidson_ls"` subsystem solver. Unless `stability_hess_x` is set, the internal stability check diagonalizes the same Hessian and inherits a non-symmetric Hessian, i.e. this setting switched off also switches off `stability_settings%hess_symm`.
- **`line_search`** (boolean): Determines whether a line search is performed after every macro iteration.
- **`refresh_hess`** (boolean): Determines whether the subsystem solver rebuilds its Hessian information after every rejected trust-region step, since the Hessian linear transformation may have changed with the objective function evaluation at the rejected point.
- **`grad_noise`** (boolean): Determines whether the gradient is perturbed by a random vector before every subsystem solve. The perturbation keeps the solution from remaining exactly within a symmetry subspace of the starting point, which a Hessian with no negative curvature outside that subspace, or no information about it, cannot otherwise leave. Its norm is the micro-iteration reduction factor (`local_red_factor` or `global_red_factor`) of the region the previous subproblem was solved in, the first counting as global, and always `local_red_factor` for `"tcg"`, which always solves to it, times the gradient norm, so that it goes to zero together with the gradient at convergence; the micro iterations are never asked to reduce the residual below this relative size, so that they do not resolve the perturbation itself, whichever region the subproblem turns out to be solved in.
- **`subsystem_solver`** (string): Specifies which subsystem solver to use. Options include:
  - `"davidson_ls"`: generalized Davidson method applied to the linear system,
  - `"jacobi-davidson_ls"`: generalized Davidson method applied to the linear system with fallback to Jacobi-Davidson if convergence is difficult, or automatically after `jacobi_davidson_start` micro iterations,
  - `"davidson_ah"`: Davidson method applied to the diagonalization of the augmented Hessian,
  - `"jacobi-davidson_ah"`: Davidson method applied to the diagonalization of the augmented Hessian with fallback to Jacobi-Davidson if convergence is difficult, or automatically after `jacobi_davidson_start` micro iterations,
  - `"tcg"`: truncated conjugate gradient method,
  - `"gltr"`: generalized Lanczos trust region method.
- **`conv_tol`** (real): Specifies the convergence criterion for the RMS gradient, which a custom `conv_check` replaces.
- **`n_random_trial_vectors`** (integer): Number of random trial vectors used to initialize the micro iterations of the Davidson-family subsystem solvers.
- **`n_extra_trial_vectors`** (integer): Number of non-random trial vectors added alongside the gradient direction to initialize the micro iterations of the Davidson-family subsystem solvers, on top of the `n_random_trial_vectors` random ones. These are taken from `get_extra_trial_vectors` when that callback is provided, and are otherwise unit vectors along the lowest Hessian diagonal elements, added only where those indicate negative curvature. Vectors that vanish or are linearly dependent on the gradient direction are dropped.
- **`start_trust_radius`** (real): Initial trust radius. Defaults to 0.4 for a spherical and to 1.0 for an ellipsoidal trust region.
- **`trust_region_shape`** (string): Shape of the trust region. The Davidson-family subsystem solvers only support a spherical trust region, so the shape is only a choice for the `"tcg"` and `"gltr"` subsystem solvers. It defaults to `"spherical"` for the Davidson-family solvers and to `"ellipsoidal"` for `"tcg"` and `"gltr"`; requesting `"ellipsoidal"` for a Davidson-family solver is an error. Options include:
  - `"ellipsoidal"`: the trust region is measured in the norm induced by the positive-definite preconditioner `precond_pd` (or, by default, shaped by the absolute value of the (approximate) Hessian diagonal). This also accelerates convergence of the underlying Krylov iteration.
  - `"spherical"`: the trust region is a Euclidean sphere. For `"tcg"` this is combined with the preconditioned search direction at no extra cost. For `"gltr"`, obtaining an exact solution to the Euclidean-sphere subproblem requires falling back to no preconditioning at all, so this option sacrifices the preconditioner's acceleration of the Krylov iteration in exchange for a genuinely spherical trust region.
- **`n_macro`** (integer): Maximum number of macro iterations.
- **`n_micro`** (integer): Maximum number of micro iterations per attempted trust-region step.
- **`jacobi_davidson_start`** (integer): Number of micro iterations after which the subsystem solver switches to the Jacobi-Davidson method.
- **`global_red_factor`** (real): Reduction factor for the residual during micro iterations in the global region, also the relative size of the gradient perturbation there (see `grad_noise`).
- **`local_red_factor`** (real): Reduction factor for the residual during micro iterations in the local region, also the relative size of the gradient perturbation there (see `grad_noise`).
- **`verbose`** (integer): Controls the verbosity of output during optimization. Level 0 prints nothing, 1 prints errors, 2 also warnings, 3 also progress information and 4 also debugging information.
- **`seed`** (integer): Seed value for generating random trial vectors and the random gradient perturbation.
- **`logger`** (subroutine): Accepts a log message. Logging is otherwise routed to stdout, and error messages to stderr. Receives the host context as its last argument.
- **`context`** (unlimited polymorphic pointer in Fortran, `void *` in C, absent in Python): Opaque host data, handed back unchanged as the last argument of every callback function so that the host does not have to reach its own state through module-level variables. The library never inspects it, and the `solver` uses it only during the call. An extension factory, however, keeps the context of the settings it wires, nested stability check settings included, until the extension is released (see [Optional Extensions](#optional-extensions)), so it must stay valid that long. Two solves can therefore run at the same time, or be nested inside one another, as long as each is given its own state record and, with solver settings an extension factory wired, its own settings object. This requires a suitable build, see [Nested and Concurrent Calls](#nested-and-concurrent-calls). The factories of the extensions create a new object for every call, so solves using them can run at the same time or be nested in the same way.
- **`stability_settings`** (stability_settings_type): Settings object controlling the internal stability check that is automatically performed upon convergence when `stability` is `True` and in any case at a starting point that is already converged and at a point with an exactly vanishing gradient (see the Stability Check section below). If `stability_settings%precond`, `stability_settings%project`, `stability_settings%logger`, or `stability_settings%context` are left unset, they default to the corresponding `precond`, `project`, `logger`, and `context` supplied to `solver`. The internal stability check hands its own context to every callback function it calls, so when `stability_settings%context` is set, the Hessian linear transformation returned by `update_orbs` or `stability_hess_x` and any inherited `precond`, `project` or `logger` receive it instead of the solver's `context` and must accept it. Leaving it unset keeps the solver's `context` everywhere. `stability_settings%verbose` is raised to at least the solver's own `verbose` level, and `stability_settings%hess_symm` is switched off together with the solver's `hess_symm` unless `stability_hess_x` is set. Unless `stability_settings%conv_check` is set, the internal stability check stops as soon as an instability is certain for a symmetric Hessian (see `stop_on_instability`), and by default it uses a single random trial vector. The solver checks the nested settings, after these defaults are resolved, together with its own before the optimization starts, also without `stability`, since the internal stability check runs in any case at the points named above.

### Output
The solver keeps the state it is given up to date during the optimization, starting from the default values on every call, and leaves the final state of the optimization in it when it returns, also after an error. A callback function can therefore read the current state of the optimization: in Fortran, the state record is declared with the `target` attribute and reached through a pointer, for example one in the host context; in C, the state is copied into the host's struct before every callback function; and in Python, the state object can be read from within a callback function. The state holds the following fields:

- **`macro_iteration`** (integer): Macro iteration of the current point, 0 for the starting point.
- **`func`** (real): Objective function value at the current point.
- **`grad_rms`** (real): RMS gradient at the current point.
- **`trust_radius`** (real): Current trust radius.
- **`micro_iteration`** (integer): Number of micro iterations spent on the trust-region subproblem being solved or, until the next one starts, on the last one, including those of rejected steps, and 0 after a step along an unstable mode.
- **`jacobi_davidson_iteration`** (integer): Number of these micro iterations that were Jacobi-Davidson iterations.
- **`relative_residual`** (real): Residual of the trust-region subproblem being solved relative to the gradient it is solved for, which `grad_noise` perturbs, both measured in the norm that the stopping test of the micro iterations uses, and 0 after a step along an unstable mode.
- **`level_shift`** (real): Level shift of the last step, which the Davidson-family subsystem solvers apply to the Hessian and which is the Lagrange multiplier of the trust-region constraint for `"gltr"`, so as to keep a step within the trust region, 0 for a step inside the trust region, for a step along an unstable mode and for `"tcg"`.
- **`step_size`** (real): Size of the last step in the norm of the trust region, after the line search if one is performed. A step along an unstable mode is measured in the Euclidean norm, also for an ellipsoidal trust region.
- **`n_update_orbs`** (integer): Number of orbital updates (`update_orbs`).
- **`n_obj_func`** (integer): Number of objective function evaluations (`obj_func`), including those of the line search.
- **`n_hess_x`** (integer): Number of Hessian linear transformations (`hess_x`, and `stability_hess_x` within the internal stability check), those of an internal stability check added once it returns. The approximate Hessian linear transformations of `stability_settings%approx_hess_x`, which only warm start the stability check, are not counted.
- **`max_precision_reached`** (boolean): `true` once neither the actual nor the predicted change of the objective function in a trust-region step can be resolved at floating-point precision relative to the objective function value, once a step is rejected whose predicted change can not be resolved, as when the objective function and the orbital update disagree by more than that, or once the gradient vanishes exactly without passing the convergence check, which leaves the subsystem solvers no subproblem to solve. This ends the optimization like convergence, without an error of its own, and the returned point is the best one found; this flag lets the host program decide whether to accept it as-is, warn, or retighten and retry with a different starting point or subsystem solver. It turns `false` again if the internal stability check, which runs with `stability` set and at an exactly vanishing gradient, then finds the point unstable and the solver continues along the unstable mode, or if the convergence check (`conv_check`, or without one the criterion on the RMS gradient) passes, which takes precedence.
- **`stability_state`** (stability_state_type): State of the most recent internal stability check, see the [stability check's output](#output-1).

## Stability Check
A separate `stability_check` subroutine is available to verify whether the current solution corresponds to a minimum. It returns whether the point is stable and, optionally, the lowest Hessian eigenvalue and the corresponding eigenvector, which at an unstable point is the direction along which the objective function decreases. It takes the following arguments:

### Arguments

- **`h_diag`** (real array): Represents the Hessian diagonal at the current point.
- **`hess_x`** (subroutine): Performs Hessian-vector products at the current point:
  - Accepts a trial vector and writes the result of the Hessian transformation into an output array (real array, written in-place)
  - Returns an integer error code (0 for success, positive integers < 100 for errors)
  - Receives the host context as its last argument.
- **`stable`** (boolean): Returns whether the current point is stable, i.e. whether the lowest Hessian eigenvalue is above the stability threshold `stability_thresh` (-1e-3 by default); a less negative eigenvalue is not distinguished from numerical noise, such as that of the integration grid of a density functional.
- **`error`** (integer): An integer code indicating the success or failure of the stability check. The error code structure is explained below.
- **`settings`** (stability_settings_type, optional): Settings object which controls optional arguments as described below. The stability check only reads it. Without it, the default settings are used.
- **`state`** (stability_state_type, optional): Record that the stability check keeps up to date during its iterations and that holds its final state once it returns, as described below.
- **`kappa`** (real array, optional): If provided, receives the normalized eigenvector of the lowest Hessian eigenvalue, written in-place, whether the current point is stable or not. At an unstable point (as can be checked from `stable`), this is the descent direction along the unstable mode.

---

The following Fortran snippet demonstrates how to use the stability check interface:

```fortran
use opentrustregion, only: ip, rp, stability_settings_type, stability_state_type, &
                           hess_x_type, stability_check

real(rp), allocatable :: h_diag(:), kappa(:)
procedure(hess_x_type), pointer :: hess_x_funptr
integer(ip) :: n_param, error
logical :: stable
type(stability_settings_type) :: settings
type(stability_state_type) :: state

! set callback function pointer to existing implementation
hess_x_funptr => hess_x

! override default settings
settings%conv_tol = 1e-6_rp
settings%n_iter = 100
settings%diag_solver = "jacobi-davidson"

! hand the callback functions whatever the host needs to reach its own data
settings%context => host_data

! run stability check
call stability_check(h_diag, hess_x_funptr, stable, error, settings, state, &
                     kappa=kappa)

! read the final state
print *, "Lowest Hessian eigenvalue:", state%min_eigval
print *, "Number of Hessian linear transformations:", state%n_hess_x
```

- `hess_x_funptr` points to an existing Hessian-vector product implementation elsewhere in the program.
- `n_param` is also assumed to be defined elsewhere.
- A declared settings object holds the default settings, which can be overridden (here, `conv_tol`, `n_iter` and `diag_solver`); assigning `stability_settings_type()` resets an object to them.
- `host_data` is a variable of any type declared with the `target` attribute. Every callback function receives it as its last argument, declared as `class(*), intent(in), pointer :: context`, and recovers its own type with `select type`. The `intent(in)` applies to the pointer rather than the data, so a callback function may modify `host_data` but not point `context` elsewhere. Leaving `settings%context` unset is fine; the callback functions then receive an unassociated pointer.
- The `stable` logical output receives the result of the stability check.
- The settings, the state and the lowest eigenvector `kappa` are optional: without settings, the default settings are used, and the others are only returned if provided.
- After the call, `state` holds the final state of the stability check (here, the lowest eigenvalue `min_eigval` and `n_hess_x`), see [Output](#output-1). The settings are left unchanged.

---

The following C snippet demonstrates the equivalent usage through the C interface:

```c
#include <stdio.h>
#include <string.h>
#include "opentrustregion.h"

c_int n_param;
c_bool stable;

// set callback function pointer to existing implementation
hess_x_fp hess_x_funptr = hess_x;

// initialize settings
stability_settings_type settings = stability_settings_init();
stability_state_type state;

// override default settings
settings.conv_tol = 1e-6;
settings.n_iter = 100;
strcpy(settings.diag_solver, "jacobi-davidson");

// hand the callback functions whatever the host needs to reach its own data
settings.context = &host_data;

// pointers to Hessian diagonal and lowest eigenvector
c_real* h_diag;
c_real* kappa;

// run stability check
c_int error = stability_check(h_diag, hess_x_funptr, n_param, &stable, &settings,
                              &state, kappa);

// read the final state
printf("Lowest Hessian eigenvalue: %f\n", state.min_eigval);
printf("Number of Hessian linear transformations: %lld\n", (long long)state.n_hess_x);
```

- `hess_x_funptr` points to an existing Hessian-vector product implementation elsewhere in the program.
- `n_param` and `h_diag` are assumed to be defined elsewhere.
- `host_data` is any host object. Every callback function receives `&host_data` as its `void *context` argument and casts it back. Leaving `settings.context` as `NULL` is fine; the callback functions then receive `NULL`.
- Stability settings are initialized via a small helper function `stability_settings_init()`, which returns a struct with default values; individual settings (here, `conv_tol`, `n_iter` and `diag_solver`) can then be overridden.
- The `stable` output receives the result of the stability check which directly returns an error code in typical C fashion.
- `kappa` receives the lowest eigenvector and can be `NULL` to skip it. So can `state`, and `settings` can be `NULL` to use the default settings.
- After the call, `state` holds the final state of the stability check (here, the lowest eigenvalue `min_eigval` and `n_hess_x`), also after an error, see [Output](#output-1). The state struct needs no initialization, and the settings are left unchanged.

---

The following Python snippet demonstrates the equivalent usage through the Python interface:

```python
import numpy as np
from pyopentrustregion import StabilitySettings, StabilityState, stability_check

# initialize settings and state
settings = StabilitySettings()
state = StabilityState()

# override default settings
settings.conv_tol = 1e-6
settings.n_iter = 100
settings.diag_solver = "jacobi-davidson"

# Hessian diagonal and lowest eigenvector arrays
h_diag = np.asarray(h_diag, dtype=np.float64)
kappa = np.empty(n_param, dtype=np.float64)

# run stability check
stable = stability_check(h_diag, hess_x, n_param, settings, state, kappa=kappa)

# read the final state
print(f"Lowest Hessian eigenvalue: {state.min_eigval}")
print(f"Number of Hessian linear transformations: {state.n_hess_x}")
```

- `hess_x` is an existing Hessian-vector product implementation elsewhere in the program.
- `n_param` and `h_diag` are assumed to be defined elsewhere.
- Stability settings are initialized via the `StabilitySettings` class, which returns an object with default values; individual settings (here, `conv_tol`, `n_iter` and `diag_solver`) can then be overridden.
- There is no `context` setting in Python. Any callable works as a callback function, so a closure, a bound method or a `functools.partial` already carries whatever host data the callback function needs.
- The `stable` output receives the result of the stability check and errors can be caught in pythonic fashion in the form of a `RuntimeError`.
- The lowest eigenvector is written into `kappa` if it is provided.
- The settings and the state are optional: without settings, the default settings are used, and the state is only needed to read it.
- After the call, `state` holds the final state of the stability check (here, the lowest eigenvalue `min_eigval` and `n_hess_x`), see [Output](#output-1), whose fields are read-only attributes. The settings are left unchanged.

### Optional Settings
The stability check can be fine-tuned using the following settings:

- **`precond`** (subroutine): Applies a preconditioner to a residual vector, given a level shift `mu`, approximating `(H - mu * I)^-1`; the Davidson iterations of the stability check pass `stability_thresh` as `mu`, while its Jacobi-Davidson iterations do not call it. Writes the result in-place into a provided array and returns an integer error code (0 for success, positive integers < 100 for errors). Receives the host context as its last argument.
- **`project`** (subroutine): Applies a projection in-place to a provided vector and returns an integer error code (0 for success, positive integers < 100 for errors). Required for stability check using non-redundant parameters. When this is used, all other passed routines (`hess_x`, `precond` and `approx_hess_x`) must be self-projecting. Receives the host context as its last argument.
- **`conv_check`** (function): Returns whether the stability check has converged due to some supplied convergence criterion, which replaces the criterion on the RMS residual (`conv_tol`) and `stop_on_instability`. It has the same arguments as the convergence check of the solver and can read the `state` passed to `stability_check`, or, for the internal stability check, from `stability_state` of the state passed to `solver`, so a host that supplies it has to pass a state. Additionally, outputs an integer code indicating the success or failure of the function, positive integers less than 100 represent error conditions. Receives the host context as its last argument.
- **`approx_hess_x`** (subroutine): Applies an approximate Hessian linear transformation to a trial vector and returns an integer error code (0 for success, positive integers < 100 for errors). Used solely to obtain the non-random starting vectors counted by `n_extra_trial_vectors`, which warm start the Davidson iterations. Has to be set explicitly, and only makes sense for a Hessian that really is approximate, since the stability check would otherwise be seeded with a converged eigenpair of the very operator it diagonalizes and report whichever eigenvalue it was handed. Receives the host context as its last argument.
- **`get_extra_trial_vectors`** (subroutine): Returns extra trial vectors that seed the leading part of the initial trial space, written in-place to the provided matrix, whose column count is `n_extra_trial_vectors` and is passed to the callback itself by the C and Python interfaces, since those only receive a flat pointer. In Python, the matrix holds one trial vector per row instead. The vectors do not need to be orthonormalized or projected, and a vanishing column is dropped, which is how the callback signals that it has no vector to contribute for that slot. Additionally, outputs an integer code indicating the success or failure of the function, positive integers less than 100 represent error conditions. Receives the host context as its last argument.
- **`hess_symm`** (boolean): Determines whether the supplied Hessian is symmetric. This is sometimes not the case for approximate Hessians. A non-symmetric Hessian is supported neither by the `"jacobi-davidson"` diagonalization solver nor by `stop_on_instability`.
- **`stop_on_instability`** (boolean): When no `conv_check` is supplied, lets the default convergence policy accept a Ritz pair as soon as its eigenvalue drops below `stability_thresh`, without waiting for the residual to converge, since for a symmetric Hessian a Ritz value is an upper bound on the true lowest eigenvalue and the sign is then already certain. Not supported for a non-symmetric Hessian (`hess_symm`). Off by default for a standalone stability check. The solver switches it on for its internal check, also when it is set to false, unless `stability_settings%conv_check` is supplied or the Hessian is non-symmetric.
- **`diag_solver`** (string): Specifies which diagonalization solver to use. Options include:
  - `"davidson"`: standard Davidson method,
  - `"jacobi-davidson"`: Davidson method that switches to Jacobi-Davidson after `jacobi_davidson_start` iterations.
- **`conv_tol`** (real): Convergence criterion for the RMS residual.
- **`stability_thresh`** (real): Threshold on the lowest Hessian eigenvalue, below which a point counts as unstable (see `stable`) and at which `stop_on_instability` accepts a Ritz pair. It must not be positive, since a positive threshold would count a minimum as unstable.
- **`n_random_trial_vectors`** (integer): Number of random trial vectors used to start the Davidson iterations. A negative value, the default, selects the number automatically: 20 for a standalone stability check and 1 for the solver's internal one. The default settings in C (`stability_settings_init()`) and Python (`StabilitySettings()`) therefore show -1, and unlike the solver's `n_random_trial_vectors`, which must not be negative, any negative value means automatic.
- **`n_extra_trial_vectors`** (integer): Number of non-random trial vectors used to start the Davidson iterations, on top of the `n_random_trial_vectors` random ones. These non-random vectors are taken from `get_extra_trial_vectors` when that callback is provided, otherwise they are the lowest eigenvectors of the approximate Hessian obtained from `approx_hess_x`, and otherwise unit vectors along the lowest Hessian diagonal elements. Vectors that vanish or are linearly dependent are dropped, so the trial space holds the remaining ones plus the random ones, at most as many as there are parameters, and the stability check fails if nothing remains.
- **`n_iter`** (integer): Maximum number of Davidson iterations.
- **`jacobi_davidson_start`** (integer): Number of iterations after which the diagonalization solver switches to the Jacobi-Davidson method.
- **`verbose`** (integer): Controls the verbosity of output during the stability check. Level 0 prints nothing, 1 prints errors, 2 also warnings, 3 also progress information and 4 also debugging information.
- **`seed`** (integer): Seed value for generating random trial vectors.
- **`logger`** (subroutine): Accepts a log message. Logging is otherwise routed to stdout, and error messages to stderr. Receives the host context as its last argument.
- **`context`** (unlimited polymorphic pointer in Fortran, `void *` in C, absent in Python): Opaque host data, handed back unchanged as the last argument of every callback function so that the host does not have to reach its own state through module-level variables. The library never inspects it, and `stability_check` uses it only during the call. Nested stability check settings are the exception: once an extension factory has wired the solver settings that hold them, their context stays in use until the extension is released. Two stability checks can therefore run at the same time, or be nested inside one another, as long as each is given its own state record. This requires a suitable build, see [Nested and Concurrent Calls](#nested-and-concurrent-calls).

### Output
The stability check keeps the state it is given up to date during its iterations, starting from the default values on every call, and leaves its final state in it when it returns, also after an error. As for the solver, a callback function can therefore read the current state. The state holds the following fields:

- **`iteration`** (integer): Current iteration of the diagonalization.
- **`min_eigval`** (real): Current estimate of the lowest Hessian eigenvalue, once the check has converged the lowest eigenvalue.
- **`residual_rms`** (real): RMS of the residual of this estimate.
- **`n_hess_x`** (integer): Number of Hessian linear transformations (`hess_x`). The approximate Hessian linear transformations of `approx_hess_x` are not counted.

## Error Code Structure

The library uses structured integer return codes to indicate whether a function has encountered an error. These codes follow the format **`OOEE`**, where:

- **`OO`** = Origin of the error (which component/function reported the error)
- **`EE`** = Specific error code

### General Rules

- A return code of `0` means success.
- Return codes between `1` and `99` are currently unused.
- All current error codes start from `100` and follow the `OOEE` structure.

### Origins (`OO`)

| Code Prefix (`OO`) | Component                                                                  |
|--------------------|----------------------------------------------------------------------------|
| `01`               | `solver`                                                                   |
| `02`               | `stability_check`                                                          |
| `11`               | `obj_func`                                                                 |
| `12`               | `update_orbs`                                                              |
| `13`               | `hess_x`, also `stability_hess_x`                                          |
| `14`               | `precond`                                                                  |
| `15`               | `conv_check`, of the solver and of the stability check                     |
| `16`               | `project`                                                                  |
| `17`               | `get_extra_trial_vectors`, of the solver and of the stability check        |
| `18`               | `precond_pd`                                                               |
| `19`               | `approx_hess_x`                                                            |
| `51`               | MO extension                                                               |
| `52`               | OAO extension                                                              |
| `53`               | Density QN extension                                                       |
| `54`               | Rotation QN extension                                                      |
| `55`               | S-GEK extension                                                            |

The factory of an extension tags the errors it detects, in its arguments, in the solver settings it wires or, in C, a missing record of the state, with the origin of its extension, and its object the errors it detects itself while performing a role, such as a missing density matrix evaluating function or a failed eigendecomposition. An error in the numerics of the orbital basis object a density QN object builds on carries the origin of the orbital basis. Errors of your own callback functions keep the origin of the role, such as `1201` for an `evaluate_dm` that fails during the orbital update.

### Error Codes (`EE`)

The error field (`EE`) is `01` for a general, unspecified error. Some origins define additional, more specific codes where the host program can plausibly react differently to them (e.g. by adjusting a setting and retrying):

| Error Code | Meaning                                                                |
|------------|------------------------------------------------------------------------|
| `0101`     | General error in `solver` |
| `0102`     | Orbital optimization did not converge within the maximum number of macro iterations (`n_macro`) |
| `0201`     | General error in `stability_check` |
| `0202`     | Stability check did not converge, within the maximum number of iterations (`n_iter`) or because its trial space could not be expanded any further |

Future versions may define more specific codes for other actionable failure modes.

### Example Error Codes

| Error Code | Meaning                   |
|------------|---------------------------|
| `0101`     | General error in `solver` |
| `1201`     | Error in `update_orbs`    |

## Optional Extensions

OpenTrustRegion (OTR) supports additional optional modules that can be enabled during installation. These include **molecular orbitals (MO)**, **orthogonal atomic orbitals (OAO)**, **rotation quasi-Newton (rotation QN)**, **density quasi-Newton (density QN)**, and **subspace gradient-enhanced kriging (S-GEK)**.

A host can wrap the functions a factory returns and the routines it wires into the solver settings, for example to count their calls. A Fortran or C wrapper passes the context it receives on to the wrapped function unchanged. Where a factory takes the functions another factory returned, as the rotation QN and S-GEK factories do on settings a factory wired, Fortran accepts only the returned functions themselves, C cannot tell a wrapper from your own function, and Python accepts a wrapper in their place only if it is made with `functools.wraps`, which marks the function it wraps:

```python
import functools


@functools.wraps(update_orbs)
def counted_update_orbs(kappa, grad, h_diag):
    counter["update_orbs"] += 1
    return update_orbs(kappa, grad, h_diag)
```

### Molecular Orbitals (MO)

This extension provides orbital optimization in the molecular orbital (MO) basis for RHF and UHF, whose parameters are only the occupied-virtual rotations of every particle channel, and only those between orbitals of the same irreducible representation (irrep) if the irreps of the MOs are given. It is the foundation of density QN's MO basis, and can also be used on its own when the exact Hessian-vector product (rather than density QN's history-based approximation) is preferred.

#### Installation

Enable the extension at build time using CMake:

```sh
cmake -DENABLE_MO=ON ..
```

or, when installing the Python package:

```sh
CMAKE_FLAGS='-DENABLE_MO=ON' pip install .
```

#### Usage

The routine `mo_factory` constructs and returns MO versions of the orbital updating subroutine and the objective function and wires the solver settings with a new MO object, which provides MO versions of the preconditioning and extra trial vector functions. Since the parameters are non-redundant, the number of parameters is the sum of `n_occ * (n_mo - n_occ)` over the particle channels, or the number of occupied-virtual pairs of the same irrep if `orbsym` is given, and the MO object provides no projection, so that the projection is left to the caller. This routine takes the following arguments:

#### Required Arguments

- **`mo_coeff`** (real array): Represents the starting MO coefficients, orthonormal with respect to the AO overlap matrix and with the occupied orbitals of every particle channel first, which are rotated in place throughout the calculation, with dimension (`n_ao`, `n_mo`) for closed-shell and (`n_ao`, `n_mo`, 2) for open-shell calculations, which is (`n_ao`, `n_mo`) and (2, `n_ao`, `n_mo`) in Python. The C interface expects the matrix of every particle channel in column-major order; the Python interface accepts any memory layout and precision of a writeable floating-point array and copies the rotated MO coefficients back after every orbital update if they cannot be rotated in place.
- **`ao_overlap`** (real array): Represents the AO overlap matrix with dimension (`n_ao`, `n_ao`).
- **`n_occ`** (integer or integer array): Specifies the number of occupied orbitals, a single integer for closed-shell and one per particle channel for open-shell calculations, which may differ between the channels.
- **`n_particle`** (integer): Specifies the number of distinct particles (1 for closed-shell, 2 for open-shell).
- **`n_ao`** (integer): Specifies the number of AOs.
- **`n_mo`** (integer): Specifies the number of MOs, which must not exceed the number of AOs.
- **`evaluate_dm`** (subroutine): Accepts an AO density matrix, built from the occupied orbitals of every particle channel, and returns the same quantities as the `evaluate_dm` of the OAO extension, whose description below applies unchanged. A host that only uses the MO object through density QN, see [Density Quasi-Newton (Density QN)](#density-quasi-newton-density-qn), and does not want an internal stability check with the exact Hessian can pass none: a disassociated procedure pointer in Fortran, `NULL` in C and `None` in Python. Without an `evaluate_dm`, the factory logs a warning at the warning verbosity of its settings, and the returned functions and the Hessian linear transformation of the object, which would need it, return an error when they are called.
- **`update_orbs_mo`** (subroutine): Returned MO orbital updating subroutine as defined for the `solver` subroutine.
- **`obj_func_mo`** (function): Returned MO objective function as defined for the `solver` subroutine.
- **`extension`** (C only, `void **`): Returns the handle of the new MO object, which `release_extension(extension)` frees once neither the returned functions nor the settings are used anymore, also not by an object wired on top of it, such as that of density QN. Fortran frees the object together with the settings, which own it, and Python once the settings and the returned functions are gone.
- **`solver_settings`** (extension_solver_settings_type in Fortran, solver_settings_type in C and Python): Solver settings, which are replaced by the default settings with a warning if a C host did not initialize them. The factory wires them with the new MO object: the routines listed below replace your callback functions, while your other callback functions stay in effect and are called with your context, those of the nested stability settings with their own context if they have one. Call the `solver` with the returned functions and these settings, or pass them on to `density_qn_factory`. Their context then belongs to the object: Fortran and C hosts pass it when they call the returned functions themselves, and callback functions they set after the factory call receive it instead of their own context. A Fortran host copies settings before a factory wires them, never after: assigning wired settings duplicates their objects, but the copy's context still points at the original objects, so the copy runs those and breaks once the original settings are reset or go out of scope. Since every call creates a new object, separate calculations, each with its own settings and state record, can run at the same time or be nested, see [Nested and Concurrent Calls](#nested-and-concurrent-calls). The MO object provides the following routines:
  - **`precond`**: MO level-shifted preconditioner, based on the exact eigendecomposition of the static part of the Hessian, also for the stability check.
  - **`precond_pd`**: MO positive-definite preconditioner, based on the exact eigendecomposition of the static part of the Hessian.
  - **`get_extra_trial_vectors`**: MO extra trial vectors, the rotations between the pseudo-canonical occupied and virtual orbitals with the most negative orbital energy differences, also for the stability check.
- **`solver_state`** (solver_state_type): Record of the solver state, through which the object has access to the state of the solve, so pass the same record to the `solver`. Fortran and C hosts keep it in place as long as the object is used, while Python keeps it alive with the object. The C and Python factories reset it to the default state, as the `solver` does, and in C, `NULL` is refused with an error.
- **`error`** (integer): An integer code indicating the success or failure of the factory, which C returns and Python raises as an exception, see [Error Code Structure](#error-code-structure).
- **`settings`** (mo_settings_type, optional): Settings object which controls optional arguments as described below. The factory only reads it. Without it, the default settings are used, which `NULL` in C and `None` in Python select.
- **`orbsym`** (integer array, optional): Specifies the irreps of the MOs as the last argument, with dimension (`n_mo`) for closed-shell and (`n_mo`, 2) for open-shell calculations, which is (`n_mo`) and (2, `n_mo`) in Python, and is omitted in Fortran, `NULL` in C and `None` in Python if all MOs belong to the same irrep. Only the occupied-virtual rotations between orbitals of the same irrep are then parameters, so that symmetry-adapted starting orbitals stay symmetry-adapted. Any integers can label the irreps; only their equality matters.

---

The following Fortran snippet demonstrates how to use the MO interface:

```fortran
use opentrustregion, only: ip, rp, solver_state_type, solver, obj_func_type, &
                           update_orbs_type
use otr_common, only: extension_solver_settings_type
use otr_mo, only: mo_factory, mo_settings_type, evaluate_dm_cs_type

type(extension_solver_settings_type), target :: settings
type(solver_state_type), target :: state
type(mo_settings_type) :: mo_settings
procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_funptr
procedure(update_orbs_type), pointer :: update_orbs_mo_funptr
procedure(obj_func_type), pointer :: obj_func_mo_funptr
integer(ip) :: n_occ, n_particle, n_ao, n_mo, n_param, error
real(rp), allocatable, target :: mo_coeff(:, :)
real(rp), allocatable :: ao_overlap(:, :)

! set callback function pointers to existing implementations
evaluate_dm_funptr => evaluate_dm

! override default settings
mo_settings%verbose = 1

! get MO routines and wire the settings with a new MO object
call mo_factory(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, &
                evaluate_dm_funptr, update_orbs_mo_funptr, obj_func_mo_funptr, &
                settings, state, error, mo_settings)

! set number of parameters
n_param = n_occ * (n_mo - n_occ)

! call solver
call solver(update_orbs_mo_funptr, obj_func_mo_funptr, n_param, error, settings, state)
```

---

The following C snippet demonstrates equivalent usage through the C interface:

```C
#include <stddef.h>
#include "opentrustregion.h"
#include "opentrustregion_mo.h"

c_int n_occ[1], n_particle, n_ao, n_mo, n_param;
c_real *mo_coeff, *ao_overlap;

// set callback function pointers to existing implementations
evaluate_dm_fp evaluate_dm_funptr = (void*)evaluate_dm;

// initialize MO settings
mo_settings_type mo_settings = mo_settings_init();

// override default settings
mo_settings.verbose = 1;

// initialize settings
solver_settings_type settings = solver_settings_init();

// get MO functions and wire the settings with a new MO object
update_orbs_fp update_orbs_mo_funptr;
obj_func_fp obj_func_mo_funptr;
void *extension;
solver_state_type state;
c_int error = mo_factory(mo_coeff,
                         ao_overlap,
                         n_occ,
                         n_particle,
                         n_ao,
                         n_mo,
                         evaluate_dm_funptr,
                         &update_orbs_mo_funptr,
                         &obj_func_mo_funptr,
                         &extension,
                         &settings,
                         &state,
                         &mo_settings,
                         NULL);

// set number of parameters
n_param = n_occ[0] * (n_mo - n_occ[0]);

// call solver
error = solver(update_orbs_mo_funptr, obj_func_mo_funptr, n_param, &settings, &state);

// release the MO object
release_extension(extension);
```

---

The following Python snippet demonstrates the equivalent usage through the Python interface:

```python
import numpy as np
from pyopentrustregion import SolverSettings, SolverState, solver
from pyopentrustregion.extensions.mo import MOSettings, mo_factory

# initialize MO settings
mo_settings = MOSettings()

# override default settings
mo_settings.verbose = 1

# initialize settings and state
settings = SolverSettings()
state = SolverState()

# get MO functions and wire the settings with a new MO object
update_orbs_mo, obj_func_mo = mo_factory(
    mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, evaluate_dm, settings, state,
    mo_settings
)

# set number of parameters
n_param = sum(occ * (n_mo - occ) for occ in np.atleast_1d(n_occ))

# call solver
solver(update_orbs_mo, obj_func_mo, n_param, settings, state)
```

---

- `mo_coeff`, `ao_overlap`, `n_occ`, `n_particle`, `n_ao` and `n_mo` are assumed to be prepared elsewhere.
- `evaluate_dm` is a callback procedure provided elsewhere.
- MO settings are initialized equivalently to the `solver` and `stability_check` settings; individual settings (here, `verbose`) can then be overridden.
- `mo_factory` returns two procedures, `update_orbs` and `obj_func`, which are passed to the normal `solver` together with the wired settings and the state record passed to the factory. To use the `solver`'s default diagonal preconditioner or no extra trial vectors instead, unset the corresponding callback function, and its stability check counterpart, after the factory call.
- Unlike density QN, MO uses the exact Hessian-vector product, obtained by calling back into the `get_response` function returned by `evaluate_dm`, rather than a history-based approximation.
- Density QN in the MO basis builds on the MO object: pass the settings `mo_factory` wired to `density_qn_factory`. Keep the exact `evaluate_dm`, since the orbital update `mo_factory` returned provides the exact Hessian for the stability check, as shown for density QN, see [Density Quasi-Newton (Density QN)](#density-quasi-newton-density-qn), and leave it out only if no stability check with the exact Hessian is wanted, which is not recommended. Fortran and C hosts call the functions `mo_factory` returned with the context the settings had right after `mo_factory`.
- The MO object is freed together with the settings in Fortran and once the settings and the returned functions are gone in Python, while a C host releases it with `release_extension`. The final orbitals are the MO coefficients passed to the factory, which are rotated in place.

#### Optional Settings
The MO factory can be fine-tuned using the following settings:

- **`verbose`** (integer): Controls the verbosity of output of the MO routines.
- **`logger`** (subroutine): Accepts a log message of the MO routines. Logging is otherwise routed to stdout, and error messages to stderr. Receives the host context as its last argument.
- **`context`** (unlimited polymorphic pointer in Fortran, `void *` in C, absent in Python): Opaque host data, handed back unchanged as the last argument of `evaluate_dm`, of the `get_response` it returns and of the logging function, as the `context` of the solver settings is to the callback functions of the `solver`. The object keeps it until the extension is released, so it has to stay valid that long. A Python host carries its data in closures instead.

### Orthogonal Atomic Orbitals (OAO)

This extension provides orbital optimization in the orthogonalized atomic orbital (OAO) basis for RHF and UHF. It is the foundation of density QN's OAO basis, and can also be used on its own when the exact Hessian-vector product (rather than density QN's history-based approximation) is preferred.

#### Installation

Enable the extension at build time using CMake:

```sh
cmake -DENABLE_OAO=ON ..
```

or, when installing the Python package:

```sh
CMAKE_FLAGS='-DENABLE_OAO=ON' pip install .
```

#### Usage

The routine `oao_factory` constructs and returns OAO versions of the orbital updating subroutine and the objective function and wires the solver settings with a new OAO object, which provides OAO versions of the preconditioning, projection, and extra trial vector functions. This routine takes the following arguments:

#### Required Arguments

- **`dm_ao`** (real array): Represents the starting AO density matrix and is updated throughout the calculation with dimension (`n_ao`, `n_ao`) for closed-shell and (`n_ao`, `n_ao`, `n_particle`) for open-shell calculations, which is (2, `n_ao`, `n_ao`) in Python, where it has to be a writeable, C-contiguous `float64` array, since it is updated in place.
- **`ao_overlap`** (real array): Represents the AO overlap matrix with dimension (`n_ao`, `n_ao`).
- **`n_particle`** (integer): Specifies the number of distinct particles (1 for closed-shell, 2 for open-shell).
- **`n_ao`** (integer): Specifies the number of AOs.
- **`evaluate_dm`** (subroutine):  
  Accepts an AO density matrix and returns:
  - Energy value (real)
  - Fock matrix (real array, written in-place)
  - A `get_response` function pointer: a callback that accepts a density matrix perturbation and returns the corresponding linear response of the Fock matrix (real array, written in-place) and an integer error code, and receives the same host context as `evaluate_dm` as its last argument. This is used to build the exact Hessian-vector product. In Python, the argument is a flag saying whether this function is needed, and the callback returns it together with the energy, or `None` in its place when it is not needed.
  - An integer error code (0 for success, positive integers < 100 for errors)
  - Receives the host context of the extension's settings as its last argument.

  Every output but the energy and the error code is optional: when the caller does not need a quantity it is passed as an absent argument in Fortran, as a null pointer in C and as `None` (or, for `get_response`, `False`) in Python, and the callback must not build it.

  A host that only uses the OAO object through density QN, see [Density Quasi-Newton (Density QN)](#density-quasi-newton-density-qn), and does not want an internal stability check with the exact Hessian can pass no `evaluate_dm`: a disassociated procedure pointer in Fortran, `NULL` in C and `None` in Python. Without an `evaluate_dm`, the factory logs a warning at the warning verbosity of its settings, and the returned functions and the Hessian linear transformation of the object, which would need it, return an error when they are called.
- **`update_orbs_oao`** (subroutine): Returned OAO orbital updating subroutine as defined for the `solver` subroutine.
- **`obj_func_oao`** (function): Returned OAO objective function as defined for the `solver` subroutine.
- **`extension`** (C only, `void **`): Returns the handle of the new OAO object, which `release_extension(extension)` frees once neither the returned functions nor the settings are used anymore, also not by an object wired on top of it, such as that of density QN. Fortran frees the object together with the settings, which own it, and Python once the settings and the returned functions are gone.
- **`solver_settings`** (extension_solver_settings_type in Fortran, solver_settings_type in C and Python): Solver settings, which are replaced by the default settings with a warning if a C host did not initialize them. The factory wires them with the new OAO object: the routines listed below replace your callback functions, apart from your projection, which the OAO projection applies after its own, while your other callback functions stay in effect, and yours are called with your context, those of the nested stability settings with their own context if they have one. Call the `solver` with the returned functions and these settings, or pass them on to `density_qn_factory`. Their context then belongs to the object: Fortran and C hosts pass it when they call the returned functions themselves, and callback functions they set after the factory call receive it instead of their own context. A Fortran host copies settings before a factory wires them, never after: assigning wired settings duplicates their objects, but the copy's context still points at the original objects, so the copy runs those and breaks once the original settings are reset or go out of scope. Since every call creates a new object, separate calculations, each with its own settings and state record, can run at the same time or be nested, see [Nested and Concurrent Calls](#nested-and-concurrent-calls). The OAO object provides the following routines:
  - **`precond`**: OAO level-shifted preconditioner, based on the exact eigendecomposition of the static part of the Hessian, also for the stability check.
  - **`precond_pd`**: OAO positive-definite preconditioner, based on the exact eigendecomposition of the static part of the Hessian.
  - **`project`**: OAO projection onto the non-redundant orbital rotations, also for the stability check, followed by your projection if you set one, e.g. a symmetry projection.
  - **`get_extra_trial_vectors`**: OAO extra trial vectors, the orbital rotations between those occupied-virtual eigenvector pairs of the static part of the Hessian whose eigenvalue sums are most negative, also for the stability check.
- **`solver_state`** (solver_state_type): Record of the solver state, through which the object has access to the state of the solve, so pass the same record to the `solver`. Fortran and C hosts keep it in place as long as the object is used, while Python keeps it alive with the object. The C and Python factories reset it to the default state, as the `solver` does, and in C, `NULL` is refused with an error.
- **`error`** (integer): An integer code indicating the success or failure of the factory, which C returns and Python raises as an exception, see [Error Code Structure](#error-code-structure).
- **`settings`** (oao_settings_type, optional): Settings object which controls optional arguments as described below. The factory only reads it. Without it, the default settings are used, which `NULL` in C and `None` in Python select.

---

The following Fortran snippet demonstrates how to use the OAO interface:

```fortran
use opentrustregion, only: ip, rp, solver_state_type, solver, obj_func_type, &
                           update_orbs_type
use otr_common, only: extension_solver_settings_type
use otr_oao, only: oao_factory, oao_settings_type, evaluate_dm_cs_type

type(extension_solver_settings_type), target :: settings
type(solver_state_type), target :: state
type(oao_settings_type) :: oao_settings
procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_funptr
procedure(update_orbs_type), pointer :: update_orbs_oao_funptr
procedure(obj_func_type), pointer :: obj_func_oao_funptr
integer(ip) :: n_particle, n_ao, n_param, error
real(rp), allocatable, target :: dm_ao(:, :)
real(rp), allocatable :: ao_overlap(:, :)

! set callback function pointers to existing implementations
evaluate_dm_funptr => evaluate_dm

! override default settings
oao_settings%verbose = 1

! get OAO routines and wire the settings with a new OAO object
call oao_factory(dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm_funptr, &
                 update_orbs_oao_funptr, obj_func_oao_funptr, settings, state, error, &
                 oao_settings)

! set number of parameters
n_param = n_particle * n_ao * (n_ao - 1) / 2

! call solver
call solver(update_orbs_oao_funptr, obj_func_oao_funptr, n_param, error, settings, &
            state)
```

---

The following C snippet demonstrates equivalent usage through the C interface:

```C
#include "opentrustregion.h"
#include "opentrustregion_oao.h"

c_int n_particle, n_ao, n_param;
c_real *dm_ao, *ao_overlap;

// set callback function pointers to existing implementations
evaluate_dm_fp evaluate_dm_funptr = (void*)evaluate_dm;

// initialize OAO settings
oao_settings_type oao_settings = oao_settings_init();

// override default settings
oao_settings.verbose = 1;

// initialize settings
solver_settings_type settings = solver_settings_init();

// get OAO functions and wire the settings with a new OAO object
update_orbs_fp update_orbs_oao_funptr;
obj_func_fp obj_func_oao_funptr;
void *extension;
solver_state_type state;
c_int error = oao_factory(dm_ao,
                          ao_overlap,
                          n_particle,
                          n_ao,
                          evaluate_dm_funptr,
                          &update_orbs_oao_funptr,
                          &obj_func_oao_funptr,
                          &extension,
                          &settings,
                          &state,
                          &oao_settings);

// set number of parameters
n_param = n_particle * n_ao * (n_ao - 1) / 2;

// call solver
error = solver(update_orbs_oao_funptr, obj_func_oao_funptr, n_param, &settings, &state);

// release the OAO object
release_extension(extension);
```

---

The following Python snippet demonstrates the equivalent usage through the Python interface:

```python
import numpy as np
from pyopentrustregion import SolverSettings, SolverState, solver
from pyopentrustregion.extensions.oao import OAOSettings, oao_factory

dm_ao = np.ascontiguousarray(dm_ao, dtype=np.float64)

# initialize OAO settings
oao_settings = OAOSettings()

# override default settings
oao_settings.verbose = 1

# initialize settings and state
settings = SolverSettings()
state = SolverState()

# get OAO functions and wire the settings with a new OAO object
update_orbs_oao, obj_func_oao = oao_factory(
    dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm, settings, state, oao_settings
)

# set number of parameters
n_param = n_particle * n_ao * (n_ao - 1) // 2

# call solver
solver(update_orbs_oao, obj_func_oao, n_param, settings, state)
```

---

- `dm_ao`, `ao_overlap`, `n_particle` and `n_ao` are assumed to be prepared elsewhere.
- `evaluate_dm` is a callback procedure provided elsewhere.
- OAO settings are initialized equivalently to the `solver` and `stability_check` settings; individual settings (here, `verbose`) can then be overridden.
- `oao_factory` returns two procedures, `update_orbs` and `obj_func`, which are passed to the normal `solver` together with the wired settings and the state record passed to the factory. To use the `solver`'s default diagonal preconditioner or no extra trial vectors instead, unset the corresponding callback function, and its stability check counterpart, after the factory call.
- Unlike density QN, OAO uses the exact Hessian-vector product, obtained by calling back into the `get_response` function returned by `evaluate_dm`, rather than a history-based approximation.
- Density QN in the OAO basis builds on the OAO object: pass the settings `oao_factory` wired to `density_qn_factory`. Keep the exact `evaluate_dm`, since the orbital update `oao_factory` returned provides the exact Hessian for the stability check, as shown for density QN, see [Density Quasi-Newton (Density QN)](#density-quasi-newton-density-qn), and leave it out only if no stability check with the exact Hessian is wanted, which is not recommended. Fortran and C hosts call the functions `oao_factory` returned with the context the settings had right after `oao_factory`.
- The OAO object is freed together with the settings in Fortran and once the settings and the returned functions are gone in Python, while a C host releases it with `release_extension`. The final AO density matrix is the one passed to the factory, which is updated in place.

#### Optional Settings
The OAO factory can be fine-tuned using the following settings:

- **`verbose`** (integer): Controls the verbosity of output of the OAO routines.
- **`logger`** (subroutine): Accepts a log message of the OAO routines. Logging is otherwise routed to stdout, and error messages to stderr. Receives the host context as its last argument.
- **`context`** (unlimited polymorphic pointer in Fortran, `void *` in C, absent in Python): Opaque host data, handed back unchanged as the last argument of `evaluate_dm`, of the `get_response` it returns and of the logging function, as the `context` of the solver settings is to the callback functions of the `solver`. The object keeps it until the extension is released, so it has to stay valid that long. A Python host carries its data in closures instead.

### Rotation Quasi-Newton (Rotation QN)

This extension provides quasi-Newton approximations of the Hessian built from the history of the orbital rotations, currently with symmetric rank-1 and Broyden-Fletcher-Goldfarb-Shanno updates.

#### Installation

Enable the extension at build time using CMake:

```sh
cmake -DENABLE_ROTATION_QN=ON ..
```

or, when installing the Python package:

```sh
CMAKE_FLAGS='-DENABLE_ROTATION_QN=ON' pip install .
```

#### Usage

This allows the solver to automatically update the Hessian using a quasi-Newton scheme, so the Hessian linear transformation the host's `update_orbs` returns is never called. The routine `rotation_qn_factory` constructs and returns rotation QN versions of the orbital updating subroutine and the objective function and wires the solver settings with a new rotation QN object. This routine takes the following arguments:

#### Required Arguments

- **`n_param`** (integer): Specifies the number of parameters to be optimized.
- **`update_orbs`** (subroutine): Original orbital updating subroutine as defined for the `solver` subroutine. The `hess_x` subroutine it returns is not used. The returned orbital updating subroutine calls it with the host context of the solver settings.
- **`obj_func`** (function): Original objective function as defined for the `solver` subroutine, which the returned objective function calls with the host context of the solver settings.
- **`transport`** (subroutine): Transports a tangent vector in place along the geodesic of an orbital rotation, so that the stored steps and gradient differences refer to the current orbitals, and returns an integer error code (0 for success, positive integers < 100 for errors). Receives the host context of the rotation QN settings as its last argument.
- **`init_hess`** (subroutine): Applies the initial Hessian approximation in place to a vector and returns an integer error code (0 for success, positive integers < 100 for errors). Receives the host context of the rotation QN settings as its last argument.
- **`update_orbs_rotation_qn`** (subroutine): Returned rotation QN orbital updating subroutine as defined for the `solver` subroutine.
- **`obj_func_rotation_qn`** (function): Returned objective function as defined for the `solver` subroutine.
- **`extension`** (C only, `void **`): Returns the handle of the new rotation QN object, which `release_extension(extension)` frees once neither the returned functions nor the settings are used anymore, also not by an object wired on top of it. Fortran frees the object together with the settings, which own it, and Python once the settings and the returned functions are gone.
- **`solver_settings`** (extension_solver_settings_type in Fortran, solver_settings_type in C and Python): Solver settings, which are replaced by the default settings with a warning if a C host did not initialize them. The factory wires them with the new rotation QN object, which replaces none of your callback functions, so these stay in effect and are called with your context, those of the nested stability settings with their own context if they have one. Call the `solver` with the returned functions and these settings. Their context then belongs to the object: Fortran and C hosts pass it when they call the returned functions themselves, and callback functions they set after the factory call receive it instead of their own context. A Fortran host copies settings before a factory wires them, never after: assigning wired settings duplicates their objects, but the copy's context still points at the original objects, so the copy runs those and breaks once the original settings are reset or go out of scope. Settings a factory wired can be passed to `rotation_qn_factory` again, which wires a new object on top of the earlier one, but only together with the functions the earlier factory returned in place of `update_orbs` and `obj_func`, in C and Python also wrappers of these, see [Optional Extensions](#optional-extensions), while settings no factory wired take your own. Fortran refuses any other function, also a wrapper, and Python any other than a wrapper made with `functools.wraps`, with an error; C cannot detect them, and your own functions would then be called with the context of the earlier object instead of yours. Since every call creates a new object, separate calculations, each with its own settings and state record, can run at the same time or be nested, see [Nested and Concurrent Calls](#nested-and-concurrent-calls).
- **`solver_state`** (solver_state_type): Record of the solver state, through which the object has access to the state of the solve, so pass the same record to the `solver`. Fortran and C hosts keep it in place as long as the object is used, while Python keeps it alive with the object. The C and Python factories reset it to the default state, as the `solver` does, and in C, `NULL` is refused with an error.
- **`error`** (integer): An integer code indicating the success or failure of the factory, which C returns and Python raises as an exception, see [Error Code Structure](#error-code-structure).
- **`settings`** (rotation_qn_settings_type, optional): Settings object which controls optional arguments as described below. The factory only reads it. Without it, the default settings are used.

---

The following Fortran snippet demonstrates how to use the rotation QN interface:

```fortran
use opentrustregion, only: ip, solver_state_type, update_orbs_type, obj_func_type, &
                           solver
use otr_common, only: extension_solver_settings_type
use otr_rotation_qn, only: rotation_qn_factory, rotation_qn_settings_type, transport_type, init_hess_type

type(extension_solver_settings_type), target :: settings
type(solver_state_type), target :: state
type(rotation_qn_settings_type) :: rotation_qn_settings
procedure(update_orbs_type), pointer :: update_orbs_funptr, update_orbs_rotation_qn_funptr
procedure(obj_func_type), pointer :: obj_func_funptr, obj_func_rotation_qn_funptr
procedure(transport_type), pointer :: transport_funptr
procedure(init_hess_type), pointer :: init_hess_funptr
integer(ip) :: n_param, error

! set callback function pointers to existing implementations
update_orbs_funptr => update_orbs
obj_func_funptr => obj_func
transport_funptr => transport
init_hess_funptr => init_hess

! override default settings
rotation_qn_settings%verbose = 1
rotation_qn_settings%qn_method = "bfgs"

! get rotation QN routines and wire the settings with a new rotation QN object
call rotation_qn_factory(n_param, update_orbs_funptr, obj_func_funptr, transport_funptr, &
                init_hess_funptr, update_orbs_rotation_qn_funptr, obj_func_rotation_qn_funptr, settings, &
                state, error, rotation_qn_settings)

! call solver
call solver(update_orbs_rotation_qn_funptr, obj_func_rotation_qn_funptr, n_param, error, settings, state)
```

---

The following C snippet demonstrates the equivalent usage through the C interface:

```c
#include <string.h>
#include "opentrustregion.h"
#include "opentrustregion_rotation_qn.h"

c_int n_param;

// set callback function pointers to existing implementations
update_orbs_fp update_orbs_funptr = (void*)update_orbs;
obj_func_fp obj_func_funptr = (void*)obj_func;
transport_fp transport_funptr = (void*)transport;
init_hess_fp init_hess_funptr = (void*)init_hess;

// initialize rotation QN settings
rotation_qn_settings_type rotation_qn_settings = rotation_qn_settings_init();

// override default settings
rotation_qn_settings.verbose = 1;
strcpy(rotation_qn_settings.qn_method, "bfgs");

// initialize settings
solver_settings_type settings = solver_settings_init();

// get rotation QN functions and wire the settings with a new rotation QN object
update_orbs_fp update_orbs_rotation_qn_funptr;
obj_func_fp obj_func_rotation_qn_funptr;
void *extension;
solver_state_type state;
c_int error = rotation_qn_factory(n_param,
                         update_orbs_funptr,
                         obj_func_funptr,
                         transport_funptr,
                         init_hess_funptr,
                         &update_orbs_rotation_qn_funptr,
                         &obj_func_rotation_qn_funptr,
                         &extension,
                         &settings,
                         &state,
                         &rotation_qn_settings);

// call solver
error = solver(update_orbs_rotation_qn_funptr, obj_func_rotation_qn_funptr, n_param, &settings, &state);

// release the rotation QN object
release_extension(extension);
```

---

The following Python snippet demonstrates the equivalent usage through the Python interface:

```python
from pyopentrustregion import SolverSettings, SolverState, solver
from pyopentrustregion.extensions.rotation_qn import RotationQNSettings, rotation_qn_factory

# initialize rotation QN settings
rotation_qn_settings = RotationQNSettings()

# override default settings
rotation_qn_settings.verbose = 1
rotation_qn_settings.qn_method = "bfgs"

# initialize settings and state
settings = SolverSettings()
state = SolverState()

# get rotation QN functions and wire the settings with a new rotation QN object
update_orbs_rotation_qn, obj_func_rotation_qn = rotation_qn_factory(
    n_param, update_orbs, obj_func, transport, init_hess, settings, state, rotation_qn_settings
)

# call solver
solver(update_orbs_rotation_qn, obj_func_rotation_qn, n_param, settings, state)
```

---

- Callback functions (`update_orbs`, `obj_func`, `transport`, `init_hess`) point to existing implementations elsewhere in the program.
- `n_param` is assumed to be defined elsewhere.
- Rotation QN settings are initialized equivalently to the `solver` and `stability_check` settings; individual settings (here, `verbose` and `qn_method`) can then be overridden.
- Rotation QN versions of the orbital updating function and the objective function are obtained by calling `rotation_qn_factory`, which also wires the solver settings with a new rotation QN object.
- The `solver` is then called with the returned functions, the wired settings and the state record passed to the factory.
- The rotation QN object is freed together with the settings in Fortran and once the settings and the returned functions are gone in Python, while a C host releases it with `release_extension`.

#### Optional Settings
The rotation QN factory function can be fine-tuned using the following settings:

- **`qn_method`** (string): Specifies which quasi-Newton updating scheme to use, in any case. Options include:
  - `"sr1"`: symmetric rank-1 update,
  - `"bfgs"`: Broyden-Fletcher-Goldfarb-Shanno update.
- **`max_points`** (integer): Maximum number of points the history of steps and gradient differences keeps, at least 0, where 0 keeps none.
- **`verbose`** (integer): Controls the verbosity of output of the rotation QN routines.
- **`logger`** (subroutine): Accepts a log message of the rotation QN routines. Logging is otherwise routed to stdout, and error messages to stderr. Receives the host context as its last argument.
- **`context`** (unlimited polymorphic pointer in Fortran, `void *` in C, absent in Python): Opaque host data, handed back unchanged as the last argument of `transport` and `init_hess` and of the logging function, as the `context` of the solver settings is to the callback functions of the `solver`. The object keeps it until the extension is released, so it has to stay valid that long. A Python host carries its data in closures instead.

### Density Quasi-Newton (Density QN)

This extension provides quasi-Newton approximations of the Hessian built from the history of the density matrices, the augmented Roothaan–Hall (ARH) method and its multisecant relatives, for RHF and UHF, with the orbitals parameterized either in the molecular orbital (MO) basis, the default, or in the orthogonalized atomic orbital (OAO) basis. The density QN factory builds on the orbital basis object which `mo_factory` or `oao_factory` wired into the solver settings last, so that setting up density QN takes two calls: the factory of the basis, which needs no exact `evaluate_dm` for density QN, and then the density QN factory, called `density_qn_factory` in C and Python and `density_qn_factory_cs` and `density_qn_factory_os` for the closed- and the open-shell `evaluate_dm` in Fortran. Both bases share the same history-based approximation of the Hessian and every QN method.

#### Installation

Enable the extension at build time using CMake:

```sh
cmake -DENABLE_DENSITY_QN=ON ..
```

or, when installing the Python package:

```sh
CMAKE_FLAGS='-DENABLE_DENSITY_QN=ON' pip install .
```

#### Usage

The density QN factory constructs and returns density QN versions of the orbital updating subroutine and the objective function and wires the solver settings with a new density QN object, which moves the orbitals of the orbital basis object and provides density QN versions of the preconditioning and extra trial vector functions. The parameters are those of the orbital basis: in the MO basis, the occupied-virtual rotations of every particle channel, or only those between orbitals of the same irrep if the MO factory was given `orbsym`, and in the OAO basis, the rotations between all orthogonalized AOs, packed as the `n_ao * (n_ao - 1) / 2` elements of an antisymmetric matrix per particle channel. The factory takes the following arguments:

#### Required Arguments

- **`evaluate_dm`** (subroutine):  
  Accepts an AO density matrix of the shell of the orbitals and returns:
  - Energy value (real)
  - Fock matrix (real array, written in-place)
  - Coulomb potential (real array, written in-place): the Coulomb potential of the density matrix of every spin channel, for closed-shell calculations that of both spin channels together
  - Exact-exchange potential (real array, written in-place): the exchange matrix of the density matrix of every spin channel, weighted by the exact-exchange fractions of the functional, and zero for a functional without exact exchange
  - Non-linear (exchange-correlation) potential (real array, written in-place)
  - An integer error code (0 for success, positive integers < 100 for errors)
  - Receives the host context of the density QN settings as its last argument.

  The linear (Coulomb and exact exchange) potential of a spin channel is thereby the Coulomb potential of both spin channels minus its own exact-exchange potential, and the non-linear potential is the rest of the effective potential, i.e. the effective potential minus this linear potential.

  Every output but the energy and the error code is optional: when the caller does not need a quantity it is passed as an absent argument in Fortran, as a null pointer in C and as `None` in Python, and the callback must not build it. The arrays have the dimensions of the density matrices of the orbital basis, as for the `evaluate_dm` of the OAO extension; Fortran hosts call `density_qn_factory_cs` with the closed-shell interface `density_qn_evaluate_dm_cs_type` and `density_qn_factory_os` with the open-shell interface `density_qn_evaluate_dm_os_type`, which has to be that of the orbitals.
- **`update_orbs_density_qn`** (subroutine): Returned density QN orbital updating subroutine as defined for the `solver` subroutine.
- **`obj_func_density_qn`** (function): Returned density QN objective function as defined for the `solver` subroutine.
- **`extension`** (C only, `void **`): Returns the handle of the new density QN object, which `release_extension(extension)` frees once neither the returned functions nor the settings are used anymore, before the extension of the orbital basis object, which the density QN object uses. Fortran frees both objects together with the settings, which own them, and Python once the settings and the returned functions are gone.
- **`solver_settings`** (extension_solver_settings_type in Fortran, solver_settings_type in C and Python): Solver settings which `mo_factory` or `oao_factory` wired last, whose orbital basis object the density QN object builds on; settings another factory or no factory wired last are refused with an error, except that a C host has to pass settings some factory wired or settings without a context, since the C interface cannot tell the context a factory wrote into the settings from a host's own context, which it would misread. The factory wires them with the new density QN object: the routines listed below replace those of the orbital basis object, while the other callback functions stay in effect. Call the `solver` with the returned functions and these settings. Their context then belongs to the density QN object: Fortran and C hosts pass it when they call the returned functions themselves, and callback functions they set after the factory call receive it instead of their own context. The factory also sets `refresh_hess`, which has to stay set since every energy evaluation changes the approximate Hessian, `n_micro` to 300, since the approximate Hessian is cheap to apply, and `hess_symm` according to `qn_method`. The density QN object provides the following routines:
  - **`precond`**: Density QN level-shifted preconditioner, which applies the exact, level-shifted inverse of the full density QN approximate Hessian, also for the stability check.
  - **`precond_pd`**: Density QN positive-definite preconditioner, that of the orbital basis object, based on the exact eigendecomposition of the static part of the Hessian.
  - **`get_extra_trial_vectors`**: Density QN extra trial vectors, those of the orbital basis object, which only depend on the static part of the Hessian, also for the stability check.

  The density QN object provides no projection and leaves it to the orbital basis object, so the projection is yours, if you set one before the orbital basis factory, in the MO basis, whose parameters are non-redundant, and that of the OAO object onto the non-redundant rotations followed by yours in the OAO basis, also for the stability check. Without `orbsym` and for `n_mo` equal to `n_ao`, the approximate Hessian of the MO basis equals that of the OAO basis in MO coordinates.
- **`solver_state`** (solver_state_type): Record of the solver state, through which the object has access to the state of the solve, so pass the same record to the orbital basis factory and the `solver`. Fortran and C hosts keep it in place as long as the object is used, while Python keeps it alive with the object. The C and Python factories reset it to the default state, as the `solver` does, and in C, `NULL` is refused with an error.
- **`error`** (integer): An integer code indicating the success or failure of the factory, which C returns and Python raises as an exception, see [Error Code Structure](#error-code-structure).
- **`settings`** (density_qn_settings_type, optional): Settings object which controls optional arguments as described below. The factory only reads it, and the orbital basis object keeps its own settings. Without it, the default settings are used, which `NULL` in C and `None` in Python select.

---

The following Fortran snippet demonstrates how to use the density QN interface in the MO basis, with the exact Hessian of the MO object for the stability check:

```fortran
use opentrustregion, only: ip, rp, solver_state_type, solver, obj_func_type, &
                           update_orbs_type, hess_x_type
use otr_common, only: extension_solver_settings_type, evaluate_dm_cs_type
use otr_mo, only: mo_factory
use otr_density_qn, only: density_qn_factory_cs, density_qn_settings_type, density_qn_evaluate_dm_cs_type

type(extension_solver_settings_type), target :: settings
type(solver_state_type), target :: state
type(density_qn_settings_type) :: density_qn_settings
procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_mo_funptr
procedure(density_qn_evaluate_dm_cs_type), pointer :: evaluate_dm_funptr
procedure(update_orbs_type), pointer :: update_orbs_mo_funptr, update_orbs_density_qn_funptr
procedure(obj_func_type), pointer :: obj_func_mo_funptr, obj_func_density_qn_funptr
procedure(hess_x_type), pointer :: hess_x_mo_funptr
integer(ip) :: n_occ, n_particle, n_ao, n_mo, n_param, error
real(rp) :: func
real(rp), allocatable, target :: mo_coeff(:, :)
real(rp), allocatable :: ao_overlap(:, :), kappa(:), grad(:), h_diag(:)

! set callback function pointers to existing implementations, the exact density matrix
! evaluating function for the MO object and that of density QN
evaluate_dm_mo_funptr => evaluate_dm_mo
evaluate_dm_funptr => evaluate_dm

! override default settings
density_qn_settings%verbose = 1

! request a stability check of the converged solution
settings%stability = .true.

! wire the settings with a new MO object
call mo_factory(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, &
                evaluate_dm_mo_funptr, update_orbs_mo_funptr, obj_func_mo_funptr, &
                settings, state, error)

! set number of parameters
n_param = n_occ * (n_mo - n_occ)

! get the exact Hessian linear transformation of the MO object, which follows the
! orbitals, for the stability check from its orbital update at the starting orbitals,
! before the density QN factory, which keeps it with the context of the MO object
allocate(kappa(n_param), grad(n_param), h_diag(n_param))
kappa = 0.0_rp
call update_orbs_mo_funptr(kappa, func, grad, h_diag, hess_x_mo_funptr, error, &
                           settings%context)
settings%stability_hess_x => hess_x_mo_funptr

! wire the settings with a new density QN object built on the MO object
call density_qn_factory_cs(evaluate_dm_funptr, update_orbs_density_qn_funptr, obj_func_density_qn_funptr, &
                    settings, state, error, density_qn_settings)

! call solver
call solver(update_orbs_density_qn_funptr, obj_func_density_qn_funptr, n_param, error, settings, &
            state)
```

---

The following C snippet demonstrates equivalent usage through the C interface:

```C
#include <stddef.h>
#include "opentrustregion.h"
#include "opentrustregion_density_qn.h"

c_int n_occ[1], n_particle, n_ao, n_mo, n_param;
c_real *mo_coeff, *ao_overlap, *kappa, *grad, *h_diag, func;

// set callback function pointers to existing implementations, the exact density
// matrix evaluating function for the MO object and that of density QN
evaluate_dm_fp evaluate_dm_mo_funptr = (void*)evaluate_dm_mo;
density_qn_evaluate_dm_fp evaluate_dm_funptr = (void*)evaluate_dm;

// initialize density QN settings
density_qn_settings_type density_qn_settings = density_qn_settings_init();

// override default settings
density_qn_settings.verbose = 1;

// initialize settings and request a stability check of the converged solution
solver_settings_type settings = solver_settings_init();
settings.stability = true;

// wire the settings with a new MO object
update_orbs_fp update_orbs_mo_funptr, update_orbs_density_qn_funptr;
obj_func_fp obj_func_mo_funptr, obj_func_density_qn_funptr;
void *extension_mo, *extension_density_qn;
solver_state_type state;
c_int error = mo_factory(mo_coeff,
                         ao_overlap,
                         n_occ,
                         n_particle,
                         n_ao,
                         n_mo,
                         evaluate_dm_mo_funptr,
                         &update_orbs_mo_funptr,
                         &obj_func_mo_funptr,
                         &extension_mo,
                         &settings,
                         &state,
                         NULL,
                         NULL);

// set number of parameters
n_param = n_occ[0] * (n_mo - n_occ[0]);

// get the exact Hessian linear transformation of the MO object, which follows the
// orbitals, for the stability check from its orbital update at the starting orbitals,
// before the density QN factory, which keeps it with the context of the MO object
hess_x_fp hess_x_mo_funptr;
error = update_orbs_mo_funptr(kappa,
                              &func,
                              grad,
                              h_diag,
                              &hess_x_mo_funptr,
                              settings.context);
settings.stability_hess_x = hess_x_mo_funptr;

// wire the settings with a new density QN object built on the MO object
error = density_qn_factory(evaluate_dm_funptr,
                    &update_orbs_density_qn_funptr,
                    &obj_func_density_qn_funptr,
                    &extension_density_qn,
                    &settings,
                    &state,
                    &density_qn_settings);

// call solver
error = solver(update_orbs_density_qn_funptr, obj_func_density_qn_funptr, n_param, &settings, &state);

// release the density QN object and then the MO object
release_extension(extension_density_qn);
release_extension(extension_mo);
```

---

The following Python snippet demonstrates the equivalent usage through the Python interface:

```python
import numpy as np
from pyopentrustregion import SolverSettings, SolverState, solver
from pyopentrustregion.extensions.mo import mo_factory
from pyopentrustregion.extensions.density_qn import DensityQNSettings, density_qn_factory

# initialize density QN settings
density_qn_settings = DensityQNSettings()

# override default settings
density_qn_settings.verbose = 1

# initialize settings and state and request a stability check of the converged solution
settings = SolverSettings()
settings.stability = True
state = SolverState()

# wire the settings with a new MO object
update_orbs_mo, obj_func_mo = mo_factory(
    mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, evaluate_dm_mo, settings, state
)

# set number of parameters
n_param = sum(occ * (n_mo - occ) for occ in np.atleast_1d(n_occ))

# get the exact Hessian linear transformation of the MO object, which follows the
# orbitals, for the stability check from its orbital update at the starting orbitals
_, hess_x_mo = update_orbs_mo(np.zeros(n_param), np.empty(n_param), np.empty(n_param))
settings.stability_hess_x = hess_x_mo

# wire the settings with a new density QN object built on the MO object
update_orbs_density_qn, obj_func_density_qn = density_qn_factory(evaluate_dm, settings, state, density_qn_settings)

# call solver
solver(update_orbs_density_qn, obj_func_density_qn, n_param, settings, state)
```

---

- `mo_coeff`, `ao_overlap`, `n_occ`, `n_particle`, `n_ao` and `n_mo` are assumed to be prepared elsewhere, in C also `kappa`, a zero vector, and `grad` and `h_diag`, each with `n_param` elements.
- `evaluate_dm_mo`, the exact density matrix evaluating function of the MO extension, and `evaluate_dm` are callback procedures provided elsewhere.
- Density QN settings are initialized equivalently to the `solver` and `stability_check` settings; individual settings (here, `verbose`) can then be overridden.
- In the OAO basis, `oao_factory` takes the place of `mo_factory`, with the number of parameters of the OAO basis; the call of the density QN factory is the same. The MO and the OAO basis build on the MO and the OAO extension, which `-DENABLE_DENSITY_QN=ON` always builds.
- The stability of the converged solution has to be decided with the exact Hessian, which the orbital basis object provides through its exact `evaluate_dm`, not with the approximate Hessian of density QN. The Hessian linear transformation the orbital update of the orbital basis factory returns follows the orbitals as density QN moves them, so it is obtained once, at the starting orbitals, and set as `stability_hess_x`. Fortran and C hosts set it before calling `density_qn_factory`, which keeps it together with the context of the orbital basis object: set after it, it receives the context of the density QN object and applies the approximate Hessian of density QN instead, without an error. Since the density QN object needs neither, the orbital basis factory can be given no `evaluate_dm` and the stability check left out, which is not recommended.
- The density QN factory returns two procedures, `update_orbs` and `obj_func`, which are passed to the normal `solver` together with the wired settings and the state record. To use the `solver`'s default diagonal preconditioner or no extra trial vectors instead, unset the corresponding callback function, and its stability check counterpart, after the factory call.
- Every QN method fits the linear (Coulomb and exact-exchange) part of the response to the whole history, whose secant pairs are exact at any distance, and the non-linear (exchange-correlation) part only to the history within two trust radii of the current point, since a non-linear secant pair is only accurate for a short step: its error grows with the step length.
- Of the five `qn_method` options below, only `"arh"` breaks Hessian symmetry; `"symm_arh"`, `"ms_psb"`, `"ms_sp"`, and the default `"ms_sr1"` are symmetric. The factory sets `hess_symm` of the solver settings accordingly.
- The density QN and orbital basis objects are freed together with the settings in Fortran and once the settings and the returned functions are gone in Python, while a C host releases the density QN extension before that of the orbital basis with `release_extension`. The final orbitals are the MO coefficients or the AO density matrix passed to the orbital basis factory, which are updated in place.

#### Optional Settings
The density QN factory can be fine-tuned using the following settings:

- **`qn_method`** (string): Specifies which QN method to use. Options include:
  - `"arh"`: standard ARH method which is not symmetric but fulfills all multisecant conditions,
  - `"symm_arh"`: simple symmetrized version of ARH which does not fulfill the multisecant conditions,
  - `"ms_psb"`: multisecant Powell symmetric Broyden which is only symmetric and fulfills all multisecant conditions for HF and other methods for which the energy is quadratic in the density matrix,
  - `"ms_sp"`: subspace-projected multisecant method,
  - `"ms_sr1"`: multisecant symmetric-rank-1 method, which treats the linear (Coulomb and exact exchange) part of the response exactly as the difference of a Coulomb and an exact-exchange multisecant term, each bounded by its positive semidefinite kernel, and the non-linear (exchange-correlation) part through a separate multisecant system regularized by the multisecant analogue of the SR1 skipping criterion.
- **`verbose`** (integer): Controls the verbosity of output of the density QN routines.
- **`logger`** (subroutine): Accepts a log message of the density QN routines. Logging is otherwise routed to stdout, and error messages to stderr. Receives the host context as its last argument.
- **`context`** (unlimited polymorphic pointer in Fortran, `void *` in C, absent in Python): Opaque host data, handed back unchanged as the last argument of `evaluate_dm` and of the logging function, as the `context` of the solver settings is to the callback functions of the `solver`. The object keeps it until the extension is released, so it has to stay valid that long. A Python host carries its data in closures instead.

### Subspace Gradient-Enhanced Kriging Extension

This extension currently provides the subspace gradient-enhanced kriging method to approximate the Hessian.

#### Installation

Enable the extension at build time using CMake:

```sh
cmake -DENABLE_S_GEK=ON ..
```

or, when installing the Python package:

```sh
CMAKE_FLAGS='-DENABLE_S_GEK=ON' pip install .
```

#### Usage

This allows the solver to automatically update the Hessian by building and evaluating a surrogate model, so the Hessian linear transformation the host's `update_orbs` returns is never called. The routine `s_gek_factory` constructs and returns S-GEK versions of the orbital updating subroutine and the objective function and wires the solver settings with a new S-GEK object. This routine takes the following arguments:

#### Required Arguments

- **`n_param`** (integer): Specifies the number of parameters to be optimized.
- **`update_orbs`** (subroutine): Original orbital updating subroutine as defined for the `solver` subroutine. The `hess_x` subroutine it returns is not used. The returned orbital updating subroutine calls it with the host context of the solver settings.
- **`obj_func`** (function): Original objective function as defined for the `solver` subroutine, which the returned objective function calls with the host context of the solver settings.
- **`change_reference`** (subroutine): Expresses the parameters, local gradients and gradients of the `n_points` stored points in place relative to a new reference point and returns an integer error code (0 for success, positive integers < 100 for errors). Receives the host context of the S-GEK settings as its last argument.
- **`update_orbs_s_gek`** (subroutine): Returned S-GEK orbital updating subroutine as defined for the `solver` subroutine.
- **`obj_func_s_gek`** (function): Returned objective function as defined for the `solver` subroutine.
- **`extension`** (C only, `void **`): Returns the handle of the new S-GEK object, which `release_extension(extension)` frees once neither the returned functions nor the settings are used anymore, also not by an object wired on top of it. Fortran frees the object together with the settings, which own it, and Python once the settings and the returned functions are gone.
- **`solver_settings`** (extension_solver_settings_type in Fortran, solver_settings_type in C and Python): Solver settings, which are replaced by the default settings with a warning if a C host did not initialize them. The factory wires them with the new S-GEK object, which replaces none of your callback functions, so these stay in effect and are called with your context, those of the nested stability settings with their own context if they have one. Call the `solver` with the returned functions and these settings. Their context then belongs to the object: Fortran and C hosts pass it when they call the returned functions themselves, and callback functions they set after the factory call receive it instead of their own context. A Fortran host copies settings before a factory wires them, never after: assigning wired settings duplicates their objects, but the copy's context still points at the original objects, so the copy runs those and breaks once the original settings are reset or go out of scope. Settings a factory wired can be passed to `s_gek_factory` again, which wires a new object on top of the earlier one, but only together with the functions the earlier factory returned in place of `update_orbs` and `obj_func`, in C and Python also wrappers of these, see [Optional Extensions](#optional-extensions), while settings no factory wired take your own. Fortran refuses any other function, also a wrapper, and Python any other than a wrapper made with `functools.wraps`, with an error; C cannot detect them, and your own functions would then be called with the context of the earlier object instead of yours. Since every call creates a new object, separate calculations, each with its own settings and state record, can run at the same time or be nested, see [Nested and Concurrent Calls](#nested-and-concurrent-calls).
- **`solver_state`** (solver_state_type): Record of the solver state, through which the object has access to the state of the solve, so pass the same record to the `solver`. Fortran and C hosts keep it in place as long as the object is used, while Python keeps it alive with the object. The C and Python factories reset it to the default state, as the `solver` does, and in C, `NULL` is refused with an error.
- **`error`** (integer): An integer code indicating the success or failure of the factory, which C returns and Python raises as an exception, see [Error Code Structure](#error-code-structure).
- **`settings`** (s_gek_settings_type, optional): Settings object which controls optional arguments as described below. The factory only reads it. Without it, the default settings are used.

---

The following Fortran snippet demonstrates how to use the S-GEK interface:

```fortran
use opentrustregion, only: ip, solver_state_type, update_orbs_type, obj_func_type, &
                           solver
use otr_common, only: extension_solver_settings_type
use otr_s_gek, only: s_gek_factory, s_gek_settings_type, change_reference_type

type(extension_solver_settings_type), target :: settings
type(solver_state_type), target :: state
type(s_gek_settings_type) :: s_gek_settings
procedure(update_orbs_type), pointer :: update_orbs_funptr, update_orbs_s_gek_funptr
procedure(obj_func_type), pointer :: obj_func_funptr, obj_func_s_gek_funptr
procedure(change_reference_type), pointer :: change_reference_funptr
integer(ip) :: n_param, error

! set callback function pointers to existing implementations
update_orbs_funptr => update_orbs
obj_func_funptr => obj_func
change_reference_funptr => change_reference

! override default settings
s_gek_settings%verbose = 1
s_gek_settings%max_points = 30

! get S-GEK routines and wire the settings with a new S-GEK object
call s_gek_factory(n_param, update_orbs_funptr, obj_func_funptr, &
                   change_reference_funptr, update_orbs_s_gek_funptr, &
                   obj_func_s_gek_funptr, settings, state, error, s_gek_settings)

! call solver
call solver(update_orbs_s_gek_funptr, obj_func_s_gek_funptr, n_param, error, &
            settings, state)
```

---

The following C snippet demonstrates the equivalent usage through the C interface:

```c
#include "opentrustregion.h"
#include "opentrustregion_s_gek.h"

c_int n_param;

// set callback function pointers to existing implementations
update_orbs_fp update_orbs_funptr = (void*)update_orbs;
obj_func_fp obj_func_funptr = (void*)obj_func;
change_reference_fp change_reference_funptr = (void*)change_reference;

// initialize S-GEK settings
s_gek_settings_type s_gek_settings = s_gek_settings_init();

// override default settings
s_gek_settings.verbose = 1;
s_gek_settings.max_points = 30;

// initialize settings
solver_settings_type settings = solver_settings_init();

// get S-GEK functions and wire the settings with a new S-GEK object
update_orbs_fp update_orbs_s_gek_funptr;
obj_func_fp obj_func_s_gek_funptr;
void *extension;
solver_state_type state;
c_int error = s_gek_factory(n_param,
                            update_orbs_funptr,
                            obj_func_funptr,
                            change_reference_funptr,
                            &update_orbs_s_gek_funptr,
                            &obj_func_s_gek_funptr,
                            &extension,
                            &settings,
                            &state,
                            &s_gek_settings);

// call solver
error =
    solver(update_orbs_s_gek_funptr, obj_func_s_gek_funptr, n_param, &settings, &state);

// release the S-GEK object
release_extension(extension);
```

---

The following Python snippet demonstrates the equivalent usage through the Python interface:

```python
from pyopentrustregion import SolverSettings, SolverState, solver
from pyopentrustregion.extensions.s_gek import SGEKSettings, s_gek_factory

# initialize S-GEK settings
s_gek_settings = SGEKSettings()

# override default settings
s_gek_settings.verbose = 1
s_gek_settings.max_points = 30

# initialize settings and state
settings = SolverSettings()
state = SolverState()

# get S-GEK functions and wire the settings with a new S-GEK object
update_orbs_s_gek, obj_func_s_gek = s_gek_factory(
    n_param, update_orbs, obj_func, change_reference, settings, state, s_gek_settings
)

# call solver
solver(update_orbs_s_gek, obj_func_s_gek, n_param, settings, state)
```

---

- Callback functions (`update_orbs`, `obj_func`, `change_reference`) point to existing implementations elsewhere in the program.
- `n_param` is assumed to be defined elsewhere.
- S-GEK settings are initialized equivalently to the `solver` and `stability_check` settings; individual settings (here, `verbose` and `max_points`) can then be overridden.
- S-GEK versions of the orbital updating function and the objective function are obtained by calling `s_gek_factory`, which also wires the solver settings with a new S-GEK object.
- The `solver` is then called with the returned functions, the wired settings and the state record passed to the factory.
- The S-GEK object is freed together with the settings in Fortran and once the settings and the returned functions are gone in Python, while a C host releases it with `release_extension`.

#### Optional Settings
The S-GEK factory function can be fine-tuned using the following settings:

- **`use_subspace`** (boolean): Specifies whether the surrogate model is constructed in a subspace.
- **`max_points`** (integer): Controls the number of points for the surrogate model and subspace construction, at least 1.
- **`verbose`** (integer): Controls the verbosity of output of the S-GEK routines.
- **`logger`** (subroutine): Accepts a log message of the S-GEK routines. Logging is otherwise routed to stdout, and error messages to stderr. Receives the host context as its last argument.
- **`context`** (unlimited polymorphic pointer in Fortran, `void *` in C, absent in Python): Opaque host data, handed back unchanged as the last argument of `change_reference` and of the logging function, as the `context` of the solver settings is to the callback functions of the `solver`. The object keeps it until the extension is released, so it has to stay valid that long. A Python host carries its data in closures instead.

## AI Usage Disclosure

Recent development of OpenTrustRegion has been assisted by AI coding agents. They have been used for code generation, refactoring, writing tests and drafting documentation.

All AI-assisted contributions are reviewed, edited where necessary and validated by the developers before they are merged. The developers make all design decisions, in particular those concerning the underlying methodology, and hold AI-assisted code to the same standards as any other code in the library.
