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
| **OpenTrustRegion_BUILD_TESTING** | `BOOL` | `ON` | Build the project’s testsuite. |
| **OpenTrustRegion_INSTALL_CMAKEDIR** | `STRING` | (auto) | Project install directory. |
| **CMAKE_BUILD_TYPE** | `STRING` | `Release` | Choose the build type (`Debug`, `Release`, etc.). |
| **INTEGER_SIZE** | `STRING` | *(auto)* | Set the integer precision to `4` (32-bit) or `8` (64-bit). Required when providing custom BLAS/LAPACK libraries. Otherwise defaults to 32-bit integers and tries to locate compatible BLAS and LAPACK libraries. Falls back to 64-bit integers if 32-bit libraries cannot be found. The resulting library name reflects the chosen integer precision (`libopentrustregion_32.*` or `libopentrustregion_64.*`). C programs compiled against the header of a library with 64-bit integers need `OTR_ILP64` defined, which linking the CMake target `OpenTrustRegion::opentrustregion` does automatically. |
| **BLAS_LIBRARIES** | `PATH` | *(auto)* | Path(s) to BLAS libraries. If not provided, CMake attempts to locate a suitable BLAS automatically. |
| **LAPACK_LIBRARIES** | `PATH` | *(auto)* | Path(s) to LAPACK libraries. If not provided, CMake attempts to locate a suitable LAPACK automatically. |
| **OpenTrustRegion_HOST_PROVIDES_BLAS** | `BOOL` | `OFF` | When enabled, OpenTrustRegion will not attempt to detect or link BLAS/LAPACK and the testsuite is automatically disabled. The calling program must provide BLAS/LAPACK routines that expose the unsuffixed symbol names (for example `ddot`, `dsyev`) with an integer width matching `INTEGER_SIZE`; otherwise linking will fail loudly. |
| **OpenTrustRegion_ENABLE_XHOST** | `BOOL` | `ON` | Optimize the Release build for the current machine's instruction set (`-march=native` for GNU, `-xHost` for Intel). Automatically disabled when cross-compiling, regardless of this setting. Turn off when building for a different machine than the one compiling (e.g. packaging/conda builds). |

### Nested and Concurrent Calls

The library keeps no state of its own between calls, so a solver or stability check call can be nested inside a callback function of another call, or run concurrently with another call in a different thread, as long as each call is given its own settings object and state record. A nested call re-enters procedures of the library that are still running, which Fortran only permits for procedures compiled for recursion, and a concurrent call needs their local variables on its own thread's stack. The library then has to be built with `-frecursive` for gfortran (implied by `-fopenmp`) or `-recursive` for Intel:

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
- **`settings`** (settings_type, optional): Settings object which controls optional arguments as described below. The solver only reads it. Without it, the default settings are used.
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

! initialize settings
call settings%init(error)

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
- Solver settings are initialized using the `init()` method of the derived type, and default settings can be overridden (here, `conv_tol` and `n_macro`).
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
- Solver settings are initialized via a small helper function `solver_settings_init()`, which returns a struct with default values. Individual settings (here, `conv_tol` and `n_macro`) can then be overridden.
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
- Solver settings are initialized via the `SolverSettings` class, which returns an object with default values; individual settings (here, `conv_tol`, and `n_macro`) can then be overridden.
- There is no `context` setting in Python. Any callable works as a callback function, so a closure, a bound method or a `functools.partial` already carries whatever host data the callback function needs.
- Finally, the `solver` is called with the callback functions, the settings and a state object, and errors can be caught in pythonic fashion in the form of a `RuntimeError`. The settings and the state are optional: without settings, the default settings are used, and the state is only needed to read it.
- After the call, `state` holds the final state of the optimization (here, `n_update_orbs`), see [Output](#output), whose fields are read-only attributes. The settings are left unchanged.

### Optional Settings
The optimization process can be fine-tuned using the following settings:

- **`precond`** (subroutine): Applies a preconditioner to a residual vector, given a level shift `mu`. Writes the result in-place into a provided array and returns an integer error code (0 for success, positive integers < 100 for errors). Used by the Davidson-family subsystem solvers to approximate `(H - mu * I)^-1`; `mu` legitimately ranges over any real value, including exactly `0.0`, so `precond` is not guaranteed to be positive definite (e.g. near a saddle point). Receives the host context as its last argument.
- **`project`** (subroutine): Applies a projection in-place to a provided vector and returns an integer error code (0 for success, positive integers < 100 for errors). Required for optimization using non-redundant parameters. When this is used, all other passed routines (`update_orbs`, `hess_x`, `precond`, and `precond_pd`) must be self-projecting. Receives the host context as its last argument.
- **`precond_pd`** (subroutine): Applies a positive-definite preconditioner to a residual vector. Writes the result in-place into a provided array and returns an integer error code (0 for success, positive integers < 100 for errors). Used by the `"tcg"` and `"gltr"` subsystem solvers to define the ellipsoidal trust-region metric (see `trust_region_shape` below); unlike `precond`, this callback takes no level shift and must always return a positive-definite result. Receives the host context as its last argument.
- **`get_extra_trial_vectors`** (subroutine): Returns extra trial vectors that seed the leading part of the initial trial space, written in-place to the provided matrix, whose column count is `n_extra_trial_vectors` and is passed to the callback itself by the C and Python interfaces, since those only receive a flat pointer. The vectors do not need to be orthonormalized or projected, and a vanishing column is dropped, which is how the callback signals that it has no vector to contribute for that slot. Additionally, outputs an integer code indicating the success or failure of the function, positive integers less than 100 represent error conditions. Receives the host context as its last argument.
- **`conv_check`** (function): Returns whether the optimization has converged due to some supplied convergence criterion, which replaces the criterion on the RMS gradient (`conv_tol`). Additionally, outputs an integer code indicating the success or failure of the function, positive integers less than 100 represent error conditions. Receives the host context as its last argument.
- **`stability_hess_x`** (subroutine): Applies a different Hessian linear transformation to a trial vector for the stability check than the one `update_orbs` returns for the optimization, and returns an integer error code (0 for success, positive integers < 100 for errors). Intended for optimization with an approximate Hessian, where the stability of a converged solution has to be decided with the exact Hessian. Receives the host context as its last argument.
- **`stability`** (boolean): Determines whether a stability check is performed upon convergence. The solver checks the stability of a starting point that is already converged and of a point with an exactly vanishing gradient in any case. From an unstable point, it steps along the unstable mode until the gradient RMS of the quadratic model is ten times larger of the current one and `conv_tol`, or, if both vanish, by the trust radius.
- **`hess_symm`** (boolean): Determines whether the supplied Hessian is symmetric. This is sometimes not the case for approximate Hessians. A non-symmetric Hessian is only supported by the `"davidson_ls"` subsystem solver.
- **`line_search`** (boolean): Determines whether a line search is performed after every macro iteration.
- **`refresh_hess`** (boolean): Determines whether the subsystem solver rebuilds its Hessian information after every rejected trust-region step, since the Hessian linear transformation may have changed with the objective function evaluation at the rejected point.
- **`grad_noise`** (boolean): Determines whether the gradient is perturbed by a random vector before every subsystem solve. The perturbation keeps the solution from remaining exactly within a symmetry subspace of the starting point, which a Hessian with no negative curvature outside that subspace, or no information about it, cannot otherwise leave. Its norm is the micro-iteration reduction factor (`local_red_factor` or `global_red_factor`) of the region the previous subproblem was solved in, the first counting as global, and always `local_red_factor` for `"tcg"`, which always solves to it, times the gradient norm, so that it goes to zero together with the gradient at convergence; the micro iterations are never asked to reduce the residual below this relative size, so that they do not resolve the perturbation itself, whichever region the subproblem turns out to be solved in.
- **`subsystem_solver`** (string): Specifies which subsystem solver to use. Options include:
  - `"davidson_ls"`: generalized Davidson method applied to linear system,
  - `"jacobi-davidson_ls"`: generalized Davidson method applied to linear system with fallback to Jacobi-Davidson if convergence is difficult, or automatically after `jacobi_davidson_start` micro iterations,
  - `"davidson_ah"`: Davidson method applied diagonalization of augmented Hessian,
  - `"jacobi-davidson_ah"`: Davidson method applied diagonalization of augmented Hessian with fallback to Jacobi-Davidson if convergence is difficult, or automatically after `jacobi_davidson_start` micro iterations,
  - `"tcg"`: truncated conjugate gradient method,
  - `"gltr"`: generalized Lanczos trust region method.
- **`conv_tol`** (real): Specifies the convergence criterion for the RMS gradient.
- **`n_random_trial_vectors`** (integer): Number of random trial vectors used to initialize the micro iterations of the Davidson-family subsystem solvers.
- **`n_extra_trial_vectors`** (integer): Number of non-random trial vectors added alongside the gradient direction to initialize the micro iterations, on top of the `n_random_trial_vectors` random ones. These are taken from `get_extra_trial_vectors` when that callback is provided, and are otherwise unit vectors along the lowest Hessian diagonal elements, added only where those indicate negative curvature. Vectors that vanish or are linearly dependent on the gradient direction are dropped.
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
- **`context`** (unlimited polymorphic pointer in Fortran, `void *` in C, absent in Python): Opaque host data, handed back unchanged as the last argument of every callback function so that the host does not have to reach its own state through module-level variables. The library never inspects it and never keeps it past the call, so it only has to stay valid for the duration of the call. Two solves can therefore run at the same time, or be nested inside one another, as long as each is given its own settings object and state record. This requires a suitable build, see [Nested and Concurrent Calls](#nested-and-concurrent-calls). The optional extensions keep their orbital objects and callbacks in module-level state, so solves using their factories can neither run at the same time nor be nested.
- **`stability_settings`** (stability_settings_type): Settings object controlling the internal stability check that is automatically performed upon convergence when `stability` is `True` and in any case at a starting point that is already converged and at a point with an exactly vanishing gradient (see the Stability Check section below). If `stability_settings%precond`, `stability_settings%project`, `stability_settings%logger`, or `stability_settings%context` are left unset, they default to the corresponding `precond`, `project`, `logger`, and `context` supplied to `solver`. The internal stability check hands its own context to every callback function it calls, so when `stability_settings%context` is set, the Hessian linear transformation returned by `update_orbs` and any inherited `precond`, `project` or `logger` receive it instead of the solver's `context` and must accept it. Leaving it unset keeps the solver's `context` everywhere. `stability_settings%verbose` is raised to at least the solver's own `verbose` level. Unless `stability_settings%conv_check` is set, the internal stability check stops as soon as an instability is certain (see `stop_on_instability`), and by default it uses a single random trial vector.

### Output
The solver keeps the state it is given up to date during the optimization, starting from the default values on every call, and leaves the final state of the optimization in it when it returns, also after an error. A callback function can therefore read the current state of the optimization: in Fortran, the state record is declared with the `target` attribute and reached through a pointer, for example one in the host context; in C, the state is copied into the host's struct before every callback function; and in Python, the state object can be read from within a callback function. The state holds the following fields:

- **`macro_iteration`** (integer): Macro iteration of the current point, 0 for the starting point.
- **`func`** (real): Objective function value at the current point.
- **`grad_rms`** (real): RMS gradient at the current point.
- **`trust_radius`** (real): Current trust radius.
- **`micro_iteration`** (integer): Number of micro iterations spent on the trust-region subproblem being solved or, until the next one starts, on the last one, including those of rejected steps.
- **`jacobi_davidson_iteration`** (integer): Number of these micro iterations that were Jacobi-Davidson iterations.
- **`relative_residual`** (real): Residual of the trust-region subproblem being solved relative to the gradient it is solved for, which `grad_noise` perturbs, both measured in the norm that the stopping test of the micro iterations uses.
- **`level_shift`** (real): Level shift of the last step, which the Davidson-family subsystem solvers apply to the Hessian and which is the Lagrange multiplier of the trust-region constraint for `"gltr"`, so as to keep a step within the trust region, 0 for a step inside the trust region, for a step along an unstable mode and for `"tcg"`.
- **`step_size`** (real): Size of the last step in the norm of the trust region, after the line search if one is performed. A step along an unstable mode is measured in the Euclidean norm, also for an ellipsoidal trust region.
- **`n_update_orbs`** (integer): Number of orbital updates (`update_orbs`).
- **`n_obj_func`** (integer): Number of objective function evaluations (`obj_func`), including those of the line search.
- **`n_hess_x`** (integer): Number of Hessian linear transformations (`hess_x`, and `stability_hess_x` within the internal stability check), those of an internal stability check added once it returns. The approximate Hessian linear transformations of `stability_settings%approx_hess_x`, which only warm start the stability check, are not counted.
- **`max_precision_reached`** (boolean): `true` once neither the actual nor the predicted change of the objective function in a trust-region step can be resolved at floating-point precision relative to the objective function value, or once the gradient vanishes exactly without passing the convergence check, which leaves the subsystem solvers no subproblem to solve. This ends the optimization like convergence, without an error of its own, and the returned point is the best one found; this flag lets the host program decide whether to accept it as-is, warn, or retighten and retry with a different starting point or subsystem solver. It turns `false` again if the internal stability check, which runs with `stability` set and at an exactly vanishing gradient, then finds the point unstable and the solver continues along the unstable mode, or if a custom `conv_check` reports convergence, which takes precedence.
- **`stability_state`** (stability_state_type): State of the most recent internal stability check, see the [stability check's output](#output-1).

## Stability Check
A separate `stability_check` subroutine is available to verify whether the current solution corresponds to a minimum. It returns whether the point is stable and, optionally, the lowest Hessian eigenvalue and the corresponding eigenvector, which at an unstable point is the direction along which the objective function decreases. It takes the following arguments:

### Arguments

- **`h_diag`** (real array): Represents the Hessian diagonal at the current point.
- **`hess_x`** (subroutine): Performs Hessian-vector products at the current point:
  - Accepts a trial vector and writes the result of the Hessian transformation into an output array (real array, written in-place)
  - Returns an integer error code (0 for success, positive integers < 100 for errors)
  - Receives the host context as its last argument.
- **`stable`** (boolean): Returns whether the current point is stable.
- **`error`** (integer): An integer code indicating the success or failure of the stability check. The error code structure is explained below.
- **`settings`** (settings_type, optional): Settings object which controls optional arguments as described below. The stability check only reads it. Without it, the default settings are used.
- **`state`** (stability_state_type, optional): Record that the stability check keeps up to date during its iterations and that holds its final state once it returns, as described below.
- **`kappa`** (real array, optional): If provided, receives the normalized eigenvector of the lowest Hessian eigenvalue, written in-place, whether the current point is stable or not. At an unstable point (as can be checked from `stable`), this is the descent direction along the unstable mode.
- **`min_eigval`** (real, optional): If provided, receives the lowest Hessian eigenvalue.

---

The following Fortran snippet demonstrates how to use the stability check interface:

```fortran
use opentrustregion, only: ip, rp, stability_settings_type, stability_state_type, &
                           hess_x_type, stability_check

real(rp), allocatable :: h_diag(:), kappa(:)
real(rp) :: min_eigval
procedure(hess_x_type), pointer :: hess_x_funptr
integer(ip) :: n_param, error
logical :: stable
type(stability_settings_type) :: settings
type(stability_state_type) :: state

! set callback function pointer to existing implementation
hess_x_funptr => hess_x

! initialize settings
call settings%init(error)

! override default settings
settings%conv_tol = 1e-6_rp
settings%n_iter = 100
settings%diag_solver = "jacobi-davidson"

! hand the callback functions whatever the host needs to reach its own data
settings%context => host_data

! run stability check
call stability_check(h_diag, hess_x_funptr, stable, error, settings, state, &
                     kappa=kappa, min_eigval=min_eigval)

! read the final state
print *, "Number of Hessian linear transformations:", state%n_hess_x
```

- `hess_x_funptr` points to an existing Hessian-vector product implementation elsewhere in the program.
- `n_param` is also assumed to be defined elsewhere.
- Stability settings are initialized via the `init()` method of the derived type and can be overridden (here, `conv_tol` and `n_iter`).
- `host_data` is a variable of any type declared with the `target` attribute. Every callback function receives it as its last argument, declared as `class(*), intent(in), pointer :: context`, and recovers its own type with `select type`. The `intent(in)` applies to the pointer rather than the data, so a callback function may modify `host_data` but not point `context` elsewhere. Leaving `settings%context` unset is fine; the callback functions then receive an unassociated pointer.
- The `stable` logical output receives the result of the stability check.
- The settings, the state and the lowest eigenvector `kappa` and eigenvalue `min_eigval` are optional: without settings, the default settings are used, and the others are only returned if provided.
- After the call, `state` holds the final state of the stability check (here, `n_hess_x`), see [Output](#output-1). The settings are left unchanged.

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

// pointers to Hessian diagonal and lowest eigenvector, and lowest eigenvalue
c_real* h_diag;
c_real* kappa;
c_real min_eigval;

// run stability check
c_int error = stability_check(h_diag, hess_x_funptr, n_param, &stable, &settings,
                              &state, kappa, &min_eigval);

// read the final state
printf("Number of Hessian linear transformations: %lld\n", (long long)state.n_hess_x);
```

- `hess_x_funptr` points to an existing Hessian-vector product implementation elsewhere in the program.
- `n_param` and `h_diag` are assumed to be defined elsewhere.
- `host_data` is any host object. Every callback function receives `&host_data` as its `void *context` argument and casts it back. Leaving `settings.context` as `NULL` is fine; the callback functions then receive `NULL`.
- Stability settings are initialized via a small helper function `stability_settings_init()`, which returns a struct with default values; individual settings (here, `conv_tol` and `n_iter`) can then be overridden.
- The `stable` output receives the result of the stability check which directly returns an error code in typical C fashion.
- `kappa` and `min_eigval` receive the lowest eigenvector and eigenvalue; either can be `NULL` to skip it. So can `state`, and `settings` can be `NULL` to use the default settings.
- After the call, `state` holds the final state of the stability check (here, `n_hess_x`), also after an error, see [Output](#output-1). The state struct needs no initialization, and the settings are left unchanged.

---

The following Python snippet demonstrates the equivalent usage through the Python interface:

```python
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
stable, min_eigval = stability_check(
    h_diag, hess_x, n_param, settings, state, kappa=kappa
)

# read the final state
print(f"Number of Hessian linear transformations: {state.n_hess_x}")
```

- `hess_x` is an existing Hessian-vector product implementation elsewhere in the program.
- `n_param` and `h_diag` are assumed to be defined elsewhere.
- Stability settings are initialized via the `StabilitySettings` class, which returns an object with default values; individual settings (here, `conv_tol`) can then be overridden.
- There is no `context` setting in Python. Any callable works as a callback function, so a closure, a bound method or a `functools.partial` already carries whatever host data the callback function needs.
- The `stable` output receives the result of the stability check and errors can be caught in pythonic fashion in the form of a `RuntimeError`.
- `stability_check` always returns the stability and the lowest eigenvalue `min_eigval`.
- The lowest eigenvector is written into `kappa` if it is provided.
- The settings and the state are optional: without settings, the default settings are used, and the state is only needed to read it.
- After the call, `state` holds the final state of the stability check (here, `n_hess_x`), see [Output](#output-1), whose fields are read-only attributes. The settings are left unchanged.

### Optional Settings
The stability check can be fine-tuned using the following settings:

- **`precond`** (subroutine): Applies a preconditioner to a residual vector. Writes the result in-place into a provided array and returns an integer error code (0 for success, positive integers < 100 for errors). Receives the host context as its last argument.
- **`project`** (subroutine): Applies a projection in-place to a provided vector and returns an integer error code (0 for success, positive integers < 100 for errors). Required for stability check using non-redundant parameters. When this is used, all other passed routines (`hess_x` and `precond`) must be self-projecting. Receives the host context as its last argument.
- **`conv_check`** (function): Returns whether the optimization has converged due to some supplied convergence criterion based on the provided residual vector and current eigenvalue estimate. Additionally, outputs an integer code indicating the success or failure of the function, positive integers less than 100 represent error conditions. Receives the host context as its last argument.
- **`approx_hess_x`** (subroutine): Applies an approximate Hessian linear transformation to a trial vector and returns an integer error code (0 for success, positive integers < 100 for errors). Used solely to obtain the non-random starting vectors counted by `n_extra_trial_vectors`, which warm start the Davidson iterations. Has to be set explicitly, and only makes sense for a Hessian that really is approximate, since the stability check would otherwise be seeded with a converged eigenpair of the very operator it diagonalizes and report whichever eigenvalue it was handed. Receives the host context as its last argument.
- **`get_extra_trial_vectors`** (subroutine): Returns extra trial vectors that seed the leading part of the initial trial space, written in-place to the provided matrix, whose column count is `n_extra_trial_vectors` and is passed to the callback itself by the C and Python interfaces, since those only receive a flat pointer. The vectors do not need to be orthonormalized or projected, and a vanishing column is dropped, which is how the callback signals that it has no vector to contribute for that slot. Additionally, outputs an integer code indicating the success or failure of the function, positive integers less than 100 represent error conditions. Receives the host context as its last argument.
- **`hess_symm`** (boolean): Determines whether the supplied Hessian is symmetric. This is sometimes not the case for approximate Hessians. A non-symmetric Hessian is not supported by the `"jacobi-davidson"` diagonalization solver.
- **`stop_on_instability`** (boolean): When no `conv_check` is supplied, lets the default convergence policy accept a Ritz pair as soon as its eigenvalue drops below the instability threshold, without waiting for the residual to converge, since a Ritz value is an upper bound on the true lowest eigenvalue and the sign is then already certain. Off by default.
- **`diag_solver`** (string): Specifies which diagonalization solver to use. Options include:
  - `"davidson"`: standard Davidson method,
  - `"jacobi-davidson"`: Davidson method that switches to Jacobi-Davidson after `jacobi_davidson_start` iterations.
- **`conv_tol`** (real): Convergence criterion for the RMS residual.
- **`n_random_trial_vectors`** (integer): Number of random trial vectors used to start the Davidson iterations. Defaults to 20 for a standalone stability check and to 1 for the solver's internal one.
- **`n_extra_trial_vectors`** (integer): Number of non-random trial vectors used to start the Davidson iterations, on top of the `n_random_trial_vectors` random ones. These non-random vectors are taken from `get_extra_trial_vectors` when that callback is provided, otherwise they are the lowest eigenvectors of the approximate Hessian obtained from `approx_hess_x`, and otherwise unit vectors along the lowest Hessian diagonal elements. Vectors that vanish or are linearly dependent are dropped, so the trial space holds the remaining ones plus the random ones, at most as many as there are parameters, and the stability check fails if nothing remains.
- **`n_iter`** (integer): Maximum number of Davidson iterations.
- **`jacobi_davidson_start`** (integer): Number of iterations after which the diagonalization solver switches to the Jacobi-Davidson method.
- **`verbose`** (integer): Controls the verbosity of output during the stability check. Level 0 prints nothing, 1 prints errors, 2 also warnings, 3 also progress information and 4 also debugging information.
- **`seed`** (integer): Seed value for generating random trial vectors.
- **`logger`** (function): Accepts a log message. Logging is otherwise routed to stdout, and error messages to stderr. Receives the host context as its last argument.
- **`context`** (unlimited polymorphic pointer in Fortran, `void *` in C, absent in Python): Opaque host data, handed back unchanged as the last argument of every callback function so that the host does not have to reach its own state through module-level variables. The library never inspects it and never keeps it past the call, so it only has to stay valid for the duration of the call. Two stability checks can therefore run at the same time, or be nested inside one another, as long as each is given its own settings object and state record. This requires a suitable build, see [Nested and Concurrent Calls](#nested-and-concurrent-calls).

### Output
The stability check keeps the state it is given up to date during its iterations, starting from the default values on every call, and leaves its final state in it when it returns, also after an error. As for the solver, a callback function can therefore read the current state. The state holds the following fields:

- **`iteration`** (integer): Current iteration of the diagonalization.
- **`min_eigval`** (real): Current estimate of the lowest Hessian eigenvalue.
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
| `13`               | `hess_x`, also `stability_hess_x` and `approx_hess_x`                      |
| `14`               | `precond`                                                                  |
| `15`               | `conv_check`, of the solver and of the stability check                     |
| `16`               | `project`                                                                  |
| `17`               | `get_extra_trial_vectors`, of the solver and of the stability check        |
| `18`               | `precond_pd`                                                               |

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

OpenTrustRegion (OTR) supports additional optional modules that can be enabled during installation. These include **molecular orbitals (MO)**, **orthogonal atomic orbitals (OAO)**, **quasi-Newton**, **augmented Roothaan–Hall (ARH)**, and **subspace gradient-enhanced kriging (S-GEK)**.

### Molecular Orbitals (MO)

This extension provides orbital optimization in the molecular orbital (MO) basis for RHF and UHF, whose parameters are only the occupied-virtual rotations of every particle channel, and only those between orbitals of the same irreducible representation (irrep) if the irreps of the MOs are given. It is the foundation of ARH's MO basis, and can also be used on its own when the exact Hessian-vector product (rather than ARH's history-based approximation) is preferred.

#### Installation

Enable the extension at build time using CMake:

```sh
cmake -DENABLE_MO=ON ..
```

or, when installing the Python package:

```sh
CMAKE_FLAGS='-DENABLE_MO=ON' pip install .
```

This exposes an `mo_factory` function that prepares MO-specific callbacks for energy and orbital updates and wires MO-specific preconditioning and extra trial vectors into the solver settings.

#### Usage

The routine `mo_factory` constructs and returns MO versions of the energy and orbital updating functions, and wires MO versions of the preconditioning and extra trial vector functions into the solver settings. Since the parameters are non-redundant, the number of parameters is the sum of `n_occ * (n_mo - n_occ)` over the particle channels, or the number of occupied-virtual pairs of the same irrep if `orbsym` is given, and no projection is wired into the solver settings, so that the projection is left to the caller. This routine requires the following input arguments:

#### Required Arguments

- **`mo_coeff`** (real array): Represents the starting MO coefficients, orthonormal with respect to the AO overlap matrix and with the occupied orbitals of every particle channel first, which are rotated in place throughout the calculation, with dimension (`n_ao`, `n_mo`) for closed-shell and (`n_ao`, `n_mo`, 2) for open-shell calculations, which is (`n_ao`, `n_mo`) and (2, `n_ao`, `n_mo`) in Python. The C interface expects the matrix of every particle channel in column-major order; the Python interface accepts any memory layout and precision of a writeable floating-point array and copies the rotated MO coefficients back after every orbital update if they cannot be rotated in place.
- **`ao_overlap`** (real array): Represents the AO overlap matrix with dimension (`n_ao`, `n_ao`).
- **`n_occ`** (integer or integer array): Specifies the number of occupied orbitals, a single integer for closed-shell and one per particle channel for open-shell calculations, which may differ between the channels.
- **`n_particle`** (integer): Specifies the number of distinct particles (1 for closed-shell, 2 for open-shell).
- **`n_ao`** (integer): Specifies the number of AOs.
- **`n_mo`** (integer): Specifies the number of MOs, which must not exceed the number of AOs.
- **`evaluate_dm`** (subroutine): Accepts an AO density matrix, built from the occupied orbitals of every particle channel, and returns the same quantities as the `evaluate_dm` of the OAO extension, whose description below applies unchanged.
- **`obj_func_mo`** (subroutine): Returned MO objective function as defined for the `solver` subroutine.
- **`update_orbs_mo`** (subroutine): Returned MO orbital updating subroutine as defined for the `solver` subroutine.
- **`settings`** (solver_settings_type): Solver settings object, which is initialized first if it is not initialized yet, and into which the following MO routines are wired:
  - **`precond`**: MO level-shifted preconditioner, based on the exact eigendecomposition of the static part of the Hessian, also for the stability check.
  - **`precond_pd`**: MO positive-definite preconditioner, based on the exact eigendecomposition of the static part of the Hessian.
  - **`get_extra_trial_vectors`**: MO extra trial vectors, the rotations between the pseudo-canonical occupied and virtual orbitals with the most negative orbital energy differences, also for the stability check.
- **`error`** (integer): An integer code indicating the success or failure of the factory. The error code structure is explained below.
- **`mo_settings`** (mo_settings_type): Settings object which controls optional arguments as described below.
- **`orbsym`** (integer array, optional): Specifies the irreps of the MOs as the last argument, with dimension (`n_mo`) for closed-shell and (`n_mo`, 2) for open-shell calculations, which is (`n_mo`) and (2, `n_mo`) in Python, and is omitted in Fortran, `NULL` in C and `None` in Python if all MOs belong to the same irrep. Only the occupied-virtual rotations between orbitals of the same irrep are then parameters, so that symmetry-adapted starting orbitals stay symmetry-adapted. Any integers can label the irreps; only their equality matters.

---

The following Fortran snippet demonstrates how to use the MO interface:

```fortran
use opentrustregion, only: ip, rp, solver_settings_type, solver, obj_func_type, &
                           update_orbs_type
use otr_mo, only: mo_factory, mo_settings_type, evaluate_dm_cs_type, &
                  mo_deconstructor

type(solver_settings_type) :: settings
type(mo_settings_type) :: mo_settings
procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_funptr
procedure(obj_func_type), pointer :: obj_func_mo_funptr
procedure(update_orbs_type), pointer :: update_orbs_mo_funptr
integer(ip) :: n_occ, n_particle, n_ao, n_mo, n_param, error
real(rp), allocatable, target :: mo_coeff(:, :)
real(rp), allocatable :: ao_overlap(:, :)

! set callback function pointers to existing implementations
evaluate_dm_funptr => evaluate_dm

! initialize MO settings
call mo_settings%init(error)

! override default settings
mo_settings%verbose = 1

! initialize settings
call settings%init(error)

! get MO routines and wire the remaining ones into the settings
call mo_factory(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, &
                evaluate_dm_funptr, obj_func_mo_funptr, update_orbs_mo_funptr, &
                settings, error, mo_settings)

! set number of parameters
n_param = n_occ * (n_mo - n_occ)

! call solver
call solver(update_orbs_mo_funptr, obj_func_mo_funptr, n_param, error, settings)

! deallocate MO objects
call mo_deconstructor()
```

---

The following C snippet demonstrates equivalent usage through the C interface:

```C
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

// get callback functions and wire the remaining ones into the settings
obj_func_fp obj_func_mo_funptr;
update_orbs_fp update_orbs_mo_funptr;
c_int error = mo_factory(mo_coeff,
                         ao_overlap,
                         n_occ,
                         n_particle,
                         n_ao,
                         n_mo,
                         evaluate_dm_funptr,
                         &obj_func_mo_funptr,
                         &update_orbs_mo_funptr,
                         &settings,
                         &mo_settings,
                         NULL);

// set number of parameters
n_param = n_occ[0] * (n_mo - n_occ[0]);

// call solver
error = solver(update_orbs_mo_funptr, obj_func_mo_funptr, n_param, &settings, NULL);

// deallocate MO objects
mo_deconstructor();
```

---

The following Python snippet demonstrates the equivalent usage through the Python interface:

```python
from pyopentrustregion import SolverSettings, solver
from pyopentrustregion.extensions.mo import MOSettings, mo_factory, mo_deconstructor

# initialize MO settings
mo_settings = MOSettings()

# override default settings
mo_settings.verbose = 1

# initialize settings
settings = SolverSettings()

# get callback functions and wire the remaining ones into the settings
obj_func_mo, update_orbs_mo = mo_factory(
    mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, evaluate_dm, settings,
    mo_settings
)

# set number of parameters
n_param = sum(occ * (n_mo - occ) for occ in np.atleast_1d(n_occ))

# call solver
solver(update_orbs_mo, obj_func_mo, n_param, settings)

# deallocate MO objects
mo_deconstructor()
```

---

- `mo_coeff`, `ao_overlap`, `n_occ`, `n_particle`, `n_ao` and `n_mo` are assumed to be prepared elsewhere.
- `evaluate_dm` is a callback procedure provided elsewhere.
- MO settings are initialized equivalently to the `solver` and `stability_check` settings; individual settings (here, `verbose`) can then be overridden. `mo_settings_type` does not add any settings of its own beyond the base ones.
- `mo_factory` returns two procedures, `obj_func` and `update_orbs`, which are passed to the normal `solver`, and wires `precond`, `precond_pd`, and `get_extra_trial_vectors` into the solver settings. All other solver settings set before the call are kept, including a projection. To use the `solver`'s own default diagonal preconditioner or no extra trial vectors instead, reset the corresponding setting and its stability check counterpart after the factory call.
- Unlike ARH, MO uses the exact Hessian-vector product, obtained by calling back into the `get_response` function returned by `evaluate_dm`, rather than a history-based approximation.
- When `mo_factory` and `arh_factory_mo` are combined for the same orbitals, for example to use the exact Hessian in the stability check of an ARH calculation, call `mo_factory` first: both share the MO objects and the C callbacks of the objective function, the preconditioners and the extra trial vectors, which belong to the factory called last. In Python, pass both factories the same `mo_coeff` array object, through which they share the MO coefficients they rotate. `arh_deconstructor` then also deallocates the MO objects.
- Clean up MO resources by calling `mo_deconstructor`; the final orbitals are the MO coefficients passed to the factory, which are rotated in place.

### Orthogonal Atomic Orbitals (OAO)

This extension provides orbital optimization in the orthogonalized atomic orbital (OAO) basis for RHF and UHF. It is the foundation of ARH's OAO basis, and can also be used on its own when the exact Hessian-vector product (rather than ARH's history-based approximation) is preferred.

#### Installation

Enable the extension at build time using CMake:

```sh
cmake -DENABLE_OAO=ON ..
```

or, when installing the Python package:

```sh
CMAKE_FLAGS='-DENABLE_OAO=ON' pip install .
```

This exposes an `oao_factory` function that prepares OAO-specific callbacks for energy and orbital updates and wires OAO-specific preconditioning, projection, and extra trial vectors into the solver settings.

#### Usage

The routine `oao_factory` constructs and returns OAO versions of the energy and orbital updating functions, and wires OAO versions of the preconditioning, projection, and extra trial vector functions into the solver settings. This routine requires the following input arguments:

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
- **`obj_func_oao`** (subroutine): Returned OAO objective function as defined for the `solver` subroutine.
- **`update_orbs_oao`** (subroutine): Returned OAO orbital updating subroutine as defined for the `solver` subroutine.
- **`settings`** (solver_settings_type): Solver settings object, which is initialized first if it is not initialized yet, and into which the following OAO routines are wired:
  - **`precond`**: OAO level-shifted preconditioner, based on the exact eigendecomposition of the static part of the Hessian, also for the stability check.
  - **`precond_pd`**: OAO positive-definite preconditioner, based on the exact eigendecomposition of the static part of the Hessian.
  - **`project`**: OAO projection onto the non-redundant orbital rotations, also for the stability check.
  - **`get_extra_trial_vectors`**: OAO extra trial vectors, the orbital rotations between those occupied-virtual eigenvector pairs of the static part of the Hessian whose eigenvalue sums are most negative, also for the stability check.
- **`error`** (integer): An integer code indicating the success or failure of the factory. The error code structure is explained below.
- **`oao_settings`** (oao_settings_type): Settings object which controls optional arguments as described below.

---

The following Fortran snippet demonstrates how to use the OAO interface:

```fortran
use opentrustregion, only: ip, rp, solver_settings_type, solver, obj_func_type, &
                           update_orbs_type
use otr_oao, only: oao_factory, oao_settings_type, evaluate_dm_cs_type, &
                   oao_deconstructor

type(solver_settings_type) :: settings
type(oao_settings_type) :: oao_settings
procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_funptr
procedure(obj_func_type), pointer :: obj_func_oao_funptr
procedure(update_orbs_type), pointer :: update_orbs_oao_funptr
integer(ip) :: n_particle, n_ao, n_param, error
real(rp), allocatable, target :: dm_ao(:, :)
real(rp), allocatable :: ao_overlap(:, :)

! set callback function pointers to existing implementations
evaluate_dm_funptr => evaluate_dm

! initialize OAO settings
call oao_settings%init(error)

! override default settings
oao_settings%verbose = 1

! initialize settings
call settings%init(error)

! get OAO routines and wire the remaining ones into the settings
call oao_factory(dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm_funptr, &
                 obj_func_oao_funptr, update_orbs_oao_funptr, settings, error, &
                 oao_settings)

! set number of parameters
n_param = n_ao * (n_ao - 1) / 2

! call solver
call solver(update_orbs_oao_funptr, obj_func_oao_funptr, n_param, error, settings)

! deallocate OAO objects
call oao_deconstructor()
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

// get callback functions and wire the remaining ones into the settings
obj_func_fp obj_func_oao_funptr;
update_orbs_fp update_orbs_oao_funptr;
c_int error = oao_factory(dm_ao, 
                          ao_overlap, 
                          n_particle, 
                          n_ao,
                          evaluate_dm_funptr,
                          &obj_func_oao_funptr,
                          &update_orbs_oao_funptr,
                          &settings,
                          &oao_settings);

// set number of parameters
n_param = n_particle * n_ao * (n_ao - 1) / 2;

// call solver
error = solver(update_orbs_oao_funptr, obj_func_oao_funptr, n_param, &settings, NULL);

// deallocate OAO objects
oao_deconstructor();
```

---

The following Python snippet demonstrates the equivalent usage through the Python interface:

```python
from pyopentrustregion import SolverSettings, solver
from pyopentrustregion.extensions.oao import OAOSettings, oao_factory, oao_deconstructor

dm_ao = np.ascontiguousarray(dm_ao, dtype=np.float64)
ao_overlap = np.asarray(ao_overlap, dtype=np.float64)

# initialize OAO settings
oao_settings = OAOSettings()

# override default settings
oao_settings.verbose = 1

# initialize settings
settings = SolverSettings()

# get callback functions and wire the remaining ones into the settings
obj_func_oao, update_orbs_oao = oao_factory(
    dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm, settings, oao_settings
)

# set number of parameters
n_param = n_particle * n_ao * (n_ao - 1) // 2

# call solver
solver(update_orbs_oao, obj_func_oao, n_param, settings)

# deallocate OAO objects
oao_deconstructor()
```

---

- `dm_ao`, `ao_overlap`, `n_particle` and `n_ao` are assumed to be prepared elsewhere.
- `evaluate_dm` is a callback procedure provided elsewhere.
- OAO settings are initialized equivalently to the `solver` and `stability_check` settings; individual settings (here, `verbose`) can then be overridden. `oao_settings_type` does not add any settings of its own beyond the base ones.
- `oao_factory` returns two procedures, `obj_func` and `update_orbs`, which are passed to the normal `solver`, and wires `precond`, `precond_pd`, `project`, and `get_extra_trial_vectors` into the solver settings. All other solver settings set before the call are kept. To use the `solver`'s own default diagonal preconditioner or no extra trial vectors instead, reset the corresponding setting and its stability check counterpart after the factory call.
- Unlike ARH, OAO uses the exact Hessian-vector product, obtained by calling back into the `get_response` function returned by `evaluate_dm`, rather than a history-based approximation.
- When `oao_factory` and `arh_factory_oao` are combined for the same density matrix, for example to use the exact Hessian in the stability check of an ARH calculation, call `oao_factory` first: both share the OAO objects and the C callbacks of the objective function, the preconditioners, the projection and the extra trial vectors, which belong to the factory called last. `arh_deconstructor` then also deallocates the OAO objects.
- Clean up OAO resources by calling `oao_deconstructor`; the final AO density matrix is the one passed to the factory, which is updated in place.

### Quasi-Newton Extension

This extension currently provides symmetric rank-1 and Broyden-Fletcher-Goldfarb-Shanno Hessian updating schemes.

#### Installation

Enable the extension at build time using CMake:

```sh
cmake -DENABLE_QUASI_NEWTON=ON ..
```

or, when installing the Python package:

```sh
CMAKE_FLAGS='-DENABLE_QUASI_NEWTON=ON' pip install .
```

#### Usage

This allows the solver to automatically update the Hessian using a quasi-Newton scheme, so the `hess_x` function is never called. The routine `qn_factory` constructs and returns a quasi-Newton version of the orbital-updating subroutine. This routine requires the following input arguments:

#### Required Arguments

- **`update_orbs`** (subroutine): Original orbital updating subroutine as defined for the `solver` subroutine. The `hess_x` subroutine is not used. The returned orbital updating subroutine calls it with the host context the solver hands to it.
- **`transport`** (subroutine): Transports a tangent vector in place along the geodesic of an orbital rotation, so that the stored steps and gradient differences refer to the current orbitals, and returns an integer error code (0 for success, positive integers < 100 for errors). Receives the host context of the quasi-Newton settings as its last argument.
- **`init_hess`** (subroutine): Applies the initial Hessian approximation in place to a vector and returns an integer error code (0 for success, positive integers < 100 for errors). Receives the host context of the quasi-Newton settings as its last argument.
- **`n_param`** (integer): Specifies the number of parameters to be optimized.
- **`qn_settings`** (qn_settings_type): Settings object which controls optional arguments as described below.
- **`update_orbs`** (subroutine):  
  Returned quasi-Newton orbital updating subroutine as defined for the `solver` subroutine.

---

The following Fortran snippet demonstrates how to use the quasi-Newton interface:

```fortran
use opentrustregion, only: settings_type, update_orbs_type, solver
use opentrustregion_quasi_newton, only: qn_factory, qn_settings_type, &
                                        update_orbs_qn_deconstructor

type(settings_type) :: settings
type(qn_settings_type) :: qn_settings
procedure(update_orbs_type), pointer :: update_orbs_funptr, update_orbs_qn_funptr
procedure(obj_func), pointer :: obj_func_funptr
integer(ip) :: n_param, error

! set callback function pointers to existing implementations
update_orbs_funptr => update_orbs
obj_func_funptr => obj_func

! initialize quasi-Newton settings
call qn_settings%init(error)

! override default settings
qn_settings%verbose = 1
qn_settings%hess_update_scheme = "sr1"

! get a quasi-Newton orbital updating subroutine
call qn_factory(update_orbs_funptr, n_param, qn_settings, update_orbs_qn_funptr)

! initialize settings
call settings%init(error)

! call solver
call solver(update_orbs_qn_funptr, obj_func_funptr, n_param, error, settings)

! clean up quasi-Newton objects
call update_orbs_qn_deconstructor()
```

---

The following C snippet demonstrates the equivalent usage through the C interface:

```c
#include <string.h>
#include "opentrustregion.h"
#include "opentrustregion_quasi_newton.h"

c_int n_param;

// set callback function pointers to existing implementations
update_orbs_fp update_orbs_funptr = (void*)update_orbs;
obj_func_fp obj_func_funptr = (void*)obj_func;

// initialize quasi-Newton settings
qn_settings_type qn_settings = qn_settings_init();

// override default settings
qn_settings.verbose = 1;
strcpy(qn_settings.hess_update_scheme, "sr1");

// get a quasi-Newton orbital updating function
update_orbs_fp update_orbs_qn_funptr;
qn_factory(update_orbs, n_param, qn_settings, &update_orbs_qn_funptr);

// initialize settings
solver_settings_type settings = solver_settings_init();

// call solver
c_int error = solver(update_orbs_qn_funptr, obj_func_funptr, n_param, &settings, NULL);

// clean up quasi-Newton objects
update_orbs_qn_deconstructor();
```

---

The following Python snippet demonstrates the equivalent usage through the Python interface:

```python
from pyopentrustregion import SolverSettings, solver
from pyopentrustregion.extensions.quasi_newton import (
    QNSettings, qn_factory, update_orbs_qn_deconstructor
)

# initialize quasi-Newton settings
qn_settings = QNSettings()

# override default settings
qn_settings.verbose = 1
qn_settings.hess_update_scheme = "sr1"

# get a quasi-Newton orbital updating function
update_orbs_qn = qn_factory(update_orbs, n_param, qn_settings)

# initialize settings
settings = SolverSettings()

# call solver
solver(update_orbs_qn, obj_func, n_param, settings)

# clean up quasi-Newton objects
update_orbs_qn_deconstructor()
```

---

- Callback functions (`update_orbs`, `obj_func`) point to existing implementations elsewhere in the program.
- `n_param` is assumed to be defined elsewhere.
- Quasi-Newton settings are initialized equivalently to the `solver` and `stability_check` settings; individual settings (here, `verbose` and `hess_update_scheme`) can then be overridden.
- A quasi-Newton orbital updating function is obtained by calling `qn_factory` with the original `update_orbs` pointer, `n_param`, and `qn_settings`.
- The `solver` is then called with the quasi-Newton-wrapped `update_orbs` function.
- Clean up quasi-Newton resources by calling `update_orbs_qn_deconstructor`.

#### Optional Settings
The quasi-Newton factory function can be fine-tuned using the following settings:

- **`hess_update_scheme`** (string): Specifies which quasi-Newton updating scheme to use. Options include:
  - `"sr1"`: symmetric rank-1 update,
  - `"bfgs"`: Broyden-Fletcher-Goldfarb-Shanno update,
- **`verbose`** (integer): Controls the verbosity of output during the stability check.
- **`context`** (unlimited polymorphic pointer in Fortran, `void *` in C, absent in Python): Opaque host data, handed back unchanged as the last argument of `transport` and `init_hess` and of the logging function, as the `context` of the solver settings is to the callback functions of the `solver`. A Python host carries its data in closures instead.

### Augmented Roothaan–Hall (ARH)

This extension provides access to the augmented Roothaan-Hall method for RHF and UHF, with the orbitals parameterized either in the molecular orbital (MO) basis (`arh_factory_mo`), the default, or in the orthogonalized atomic orbital (OAO) basis (`arh_factory_oao`), both available under these names and with the same arguments in Fortran, C and Python. Fortran and Python additionally provide `arh_factory`, which accepts the arguments of either basis factory and forwards them to the one they match, the generic interface in Fortran and a helper in Python; C, which cannot overload functions, calls the basis factories directly. Both bases share the same history-based approximation of the Hessian and every ARH type.

#### Installation

Enable the extension at build time using CMake:

```sh
cmake -DENABLE_ARH=ON ..
```

or, when installing the Python package:

```sh
CMAKE_FLAGS='-DENABLE_ARH=ON' pip install .
```

This exposes the `arh_factory_mo` and `arh_factory_oao` functions, which prepare ARH-specific callbacks for energy and orbital updates and wire ARH-specific preconditioning, extra trial vectors and, in the OAO basis, projection into the solver settings.

#### Usage

The ARH factories construct and return ARH versions of the energy and orbital updating functions, and wire ARH versions of the preconditioning and extra trial vector functions into the solver settings. Both bases take the following arguments, together with the orbitals of their basis described in the two parts below:

- **`ao_overlap`** (real array): Represents the AO overlap matrix with dimension (`n_ao`, `n_ao`).
- **`n_particle`** (integer): Specifies the number of distinct particles (1 for closed-shell, 2 for open-shell).
- **`n_ao`** (integer): Specifies the number of AOs.
- **`evaluate_dm`** (subroutine):  
  Accepts an AO density matrix and returns:
  - Energy value (real)
  - Fock matrix (real array, written in-place)
  - Coulomb potential (real array, written in-place): the Coulomb potential of the density matrix of every spin channel, for closed-shell calculations that of both spin channels together
  - Exact-exchange potential (real array, written in-place): the exchange matrix of the density matrix of every spin channel, weighted by the exact-exchange fractions of the functional, and zero for a functional without exact exchange
  - Non-linear (exchange-correlation) potential (real array, written in-place)
  - An integer error code (0 for success, positive integers < 100 for errors)
  - Receives the host context of the ARH settings as its last argument.

  The linear (Coulomb and exact exchange) potential of a spin channel is thereby the Coulomb potential of both spin channels minus its own exact-exchange potential, and the non-linear potential is the rest of the effective potential, i.e. the effective potential minus this linear potential.

  Every output but the energy and the error code is optional: when the caller does not need a quantity it is passed as an absent argument in Fortran, as a null pointer in C and as `None` in Python, and the callback must not build it.
- **`obj_func_arh`** (subroutine): Returned ARH objective function as defined for the `solver` subroutine.
- **`update_orbs_arh`** (subroutine): Returned ARH orbital updating subroutine as defined for the `solver` subroutine.
- **`settings`** (solver_settings_type): Solver settings object, which is initialized first if it is not initialized yet. Its `refresh_hess` is set, since every energy evaluation adds the evaluated point to the ARH history and thereby changes the approximate Hessian, its `n_micro` is set to 300, since the linear transformations of the approximate Hessian are cheap compared to energy evaluations and a subsystem solve stopped by the micro iteration limit is rejected, its `hess_symm` is set according to `arh_type`, and the following ARH routines are wired into it:
  - **`precond`**: ARH level-shifted preconditioner, which applies the exact, level-shifted inverse of the full ARH approximate Hessian using the Sherman-Morrison-Woodbury identity, also for the stability check.
  - **`precond_pd`**: ARH positive-definite preconditioner, based on the exact eigendecomposition of the static part of the Hessian.
  - **`get_extra_trial_vectors`**: ARH extra trial vectors, which only depend on the static part of the Hessian, also for the stability check.

  How the two bases treat the projection is described in their parts below.
- **`error`** (integer): An integer code indicating the success or failure of the factory. The error code structure is explained below.
- **`arh_settings`** (arh_settings_type): Settings object which controls optional arguments as described below.

##### Orbitals in the MO basis

In the MO basis, the parameters are only the occupied-virtual rotations of every particle channel, and only those between orbitals of the same irrep if the irreps of the MOs are given as the optional last argument `orbsym`, described for the MO extension, so that the number of parameters is the sum of `n_occ * (n_mo - n_occ)` over the particle channels, or the number of occupied-virtual pairs of the same irrep. These parameters are non-redundant, so the factory wires no projection into the solver settings and leaves the projection to the caller. `precond_pd` and `get_extra_trial_vectors` are those of the MO extension. Without `orbsym` and for `n_mo` equal to `n_ao`, the approximate Hessian equals that of the OAO basis in MO coordinates. Besides the common arguments, `arh_factory_mo` takes, in the order `mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, evaluate_dm, ...`:

- **`mo_coeff`** (real array): Represents the starting MO coefficients, orthonormal with respect to the AO overlap matrix and with the occupied orbitals of every particle channel first, which are rotated in place throughout the calculation, with dimension (`n_ao`, `n_mo`) for closed-shell and (`n_ao`, `n_mo`, 2) for open-shell calculations, which is (`n_ao`, `n_mo`) and (2, `n_ao`, `n_mo`) in Python. The C interface expects the matrix of every particle channel in column-major order; the Python interface accepts any memory layout and precision of a writeable floating-point array and copies the rotated MO coefficients back after every orbital update if they cannot be rotated in place.
- **`n_occ`** (integer or integer array): Specifies the number of occupied orbitals, a single integer for closed-shell and one per particle channel for open-shell calculations, which may differ between the channels.
- **`n_mo`** (integer): Specifies the number of MOs, which must not exceed the number of AOs.

The following Fortran snippet demonstrates how to use the ARH interface in the MO basis:

```fortran
use opentrustregion, only: ip, rp, solver_settings_type, solver, obj_func_type, &
                           update_orbs_type
use otr_arh, only: arh_factory, arh_settings_type, evaluate_dm_cs_type, &
                   arh_deconstructor

type(solver_settings_type) :: settings
type(arh_settings_type) :: arh_settings
procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_funptr
procedure(obj_func_type), pointer :: obj_func_arh_funptr
procedure(update_orbs_type), pointer :: update_orbs_arh_funptr
integer(ip) :: n_occ, n_particle, n_ao, n_mo, n_param, error
real(rp), allocatable, target :: mo_coeff(:, :)
real(rp), allocatable :: ao_overlap(:, :)

! set callback function pointers to existing implementations
evaluate_dm_funptr => evaluate_dm

! initialize ARH settings
call arh_settings%init(error)

! override default settings
arh_settings%verbose = 1

! initialize settings
call settings%init(error)

! get ARH routines and wire the remaining ones into the settings
call arh_factory(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, &
                 evaluate_dm_funptr, obj_func_arh_funptr, update_orbs_arh_funptr, &
                 settings, error, arh_settings)

! set number of parameters
n_param = n_occ * (n_mo - n_occ)

! call solver
call solver(update_orbs_arh_funptr, obj_func_arh_funptr, n_param, error, settings)

! deallocate ARH objects
call arh_deconstructor()
```

---

The following C snippet demonstrates equivalent usage through the C interface:

```C
#include "opentrustregion.h"
#include "opentrustregion_arh.h"

c_int n_occ[1], n_particle, n_ao, n_mo, n_param;
c_real *mo_coeff, *ao_overlap;

// set callback function pointers to existing implementations
arh_evaluate_dm_fp evaluate_dm_funptr = (void*)evaluate_dm;

// initialize ARH settings
arh_settings_type arh_settings = arh_settings_init();

// override default settings
arh_settings.verbose = 1;

// initialize settings
solver_settings_type settings = solver_settings_init();

// get callback functions and wire the remaining ones into the settings
obj_func_fp obj_func_arh_funptr;
update_orbs_fp update_orbs_arh_funptr;
c_int error = arh_factory_mo(mo_coeff,
                             ao_overlap,
                             n_occ,
                             n_particle,
                             n_ao,
                             n_mo,
                             evaluate_dm_funptr,
                             &obj_func_arh_funptr,
                             &update_orbs_arh_funptr,
                             &settings,
                             &arh_settings,
                             NULL);

// set number of parameters
n_param = n_occ[0] * (n_mo - n_occ[0]);

// call solver
error = solver(update_orbs_arh_funptr, obj_func_arh_funptr, n_param, &settings, NULL);

// deallocate ARH objects
arh_deconstructor();
```

---

The following Python snippet demonstrates the equivalent usage through the Python interface:

```python
from pyopentrustregion import SolverSettings, solver
from pyopentrustregion.extensions.arh import ARHSettings, arh_factory, arh_deconstructor

# initialize ARH settings
arh_settings = ARHSettings()

# override default settings
arh_settings.verbose = 1

# initialize settings
settings = SolverSettings()

# get callback functions and wire the remaining ones into the settings
obj_func_arh, update_orbs_arh = arh_factory(
    mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, evaluate_dm, settings,
    arh_settings
)

# set number of parameters
n_param = sum(occ * (n_mo - occ) for occ in np.atleast_1d(n_occ))

# call solver
solver(update_orbs_arh, obj_func_arh, n_param, settings)

# deallocate ARH objects
arh_deconstructor()
```

##### Orbitals in the OAO basis

In the OAO basis, the parameters are the rotations between all orthogonalized AOs, packed as the `n_ao * (n_ao - 1) / 2` elements of an antisymmetric matrix per particle channel. These parameters are redundant, so the factory also wires the projection onto the non-redundant rotations of the OAO extension into the solver settings, also for the stability check, replacing a projection that is already set; `precond_pd` and `get_extra_trial_vectors` are those of the OAO extension. Besides the common arguments, `arh_factory_oao` takes, in the order `dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm, ...`:

- **`dm_ao`** (real array): Represents the starting AO density matrix and is updated throughout the calculation with dimension (`n_ao`, `n_ao`) for closed-shell and (`n_ao`, `n_ao`, `n_particle`) for open-shell calculations, which is (2, `n_ao`, `n_ao`) in Python, where it has to be a writeable, C-contiguous `float64` array, since it is updated in place.

The following Fortran snippet demonstrates how to use the ARH interface in the OAO basis:

```fortran
use opentrustregion, only: ip, rp, solver_settings_type, solver, obj_func_type, &
                           update_orbs_type
use otr_arh, only: arh_factory, arh_settings_type, evaluate_dm_cs_type, &
                   arh_deconstructor

type(solver_settings_type) :: settings
type(arh_settings_type) :: arh_settings
procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_funptr
procedure(obj_func_type), pointer :: obj_func_arh_funptr
procedure(update_orbs_type), pointer :: update_orbs_arh_funptr
integer(ip) :: n_particle, n_ao, n_param, error
real(rp), allocatable, target :: dm_ao(:, :)
real(rp), allocatable :: ao_overlap(:, :)

! set callback function pointers to existing implementations
evaluate_dm_funptr => evaluate_dm

! initialize ARH settings
call arh_settings%init(error)

! override default settings
arh_settings%verbose = 1

! initialize settings
call settings%init(error)

! get ARH routines and wire the remaining ones into the settings
call arh_factory(dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm_funptr, &
                 obj_func_arh_funptr, update_orbs_arh_funptr, settings, error, &
                 arh_settings)

! set number of parameters
n_param = n_ao * (n_ao - 1) / 2

! call solver
call solver(update_orbs_arh_funptr, obj_func_arh_funptr, n_param, error, settings)

! deallocate ARH objects
call arh_deconstructor()
```

---

The following C snippet demonstrates equivalent usage through the C interface:

```C
#include "opentrustregion.h"
#include "opentrustregion_arh.h"

c_int n_particle, n_ao, n_param;
c_real *dm_ao, *ao_overlap;

// set callback function pointers to existing implementations
arh_evaluate_dm_fp evaluate_dm_funptr = (void*)evaluate_dm;

// initialize ARH settings
arh_settings_type arh_settings = arh_settings_init();

// override default settings
arh_settings.verbose = 1;

// initialize settings
solver_settings_type settings = solver_settings_init();

// get callback functions and wire the remaining ones into the settings
obj_func_fp obj_func_arh_funptr;
update_orbs_fp update_orbs_arh_funptr;
c_int error = arh_factory_oao(dm_ao,
                              ao_overlap,
                              n_particle,
                              n_ao,
                              evaluate_dm_funptr,
                              &obj_func_arh_funptr,
                              &update_orbs_arh_funptr,
                              &settings,
                              &arh_settings);

// set number of parameters
n_param = n_ao * (n_ao - 1) / 2;

// call solver
error = solver(update_orbs_arh_funptr, obj_func_arh_funptr, n_param, &settings, NULL);

// deallocate ARH objects
arh_deconstructor();
```

---

The following Python snippet demonstrates the equivalent usage through the Python interface:

```python
from pyopentrustregion import SolverSettings, solver
from pyopentrustregion.extensions.arh import ARHSettings, arh_factory, arh_deconstructor

dm_ao = np.ascontiguousarray(dm_ao, dtype=np.float64)

# initialize ARH settings
arh_settings = ARHSettings()

# override default settings
arh_settings.verbose = 1

# initialize settings
settings = SolverSettings()

# get callback functions and wire the remaining ones into the settings
obj_func_arh, update_orbs_arh = arh_factory(
    dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm, settings, arh_settings
)

# set number of parameters
n_param = n_particle * n_ao * (n_ao - 1) // 2

# call solver
solver(update_orbs_arh, obj_func_arh, n_param, settings)

# deallocate ARH objects
arh_deconstructor()
```

---

- `mo_coeff` or `dm_ao`, `ao_overlap` and the dimensions are assumed to be prepared elsewhere.
- `evaluate_dm` is a callback procedure provided elsewhere.
- ARH settings are initialized equivalently to the `solver` and `stability_check` settings; individual settings (here, `verbose`) can then be overridden.
- `arh_factory` in Fortran and Python is equivalent to calling `arh_factory_mo` or `arh_factory_oao` with the same arguments.
- The MO and the OAO basis build on the MO and the OAO extension, which `-DENABLE_ARH=ON` always builds. To combine an ARH factory with the exact factory of its basis, call `mo_factory` or `oao_factory` first, as described for those extensions.
- The factories return two procedures, `obj_func` and `update_orbs`, which are passed to the normal `solver`, wire `precond`, `precond_pd`, `get_extra_trial_vectors` and, in the OAO basis, `project` into the solver settings and set its `refresh_hess`, `n_micro`, and `hess_symm`. All other solver settings set before the call are kept. To use the `solver`'s own default diagonal preconditioner or no extra trial vectors instead, reset the corresponding setting and its stability check counterpart after the factory call. `refresh_hess` has to stay set: the ARH approximate Hessian changes with every energy evaluation, so without it the subsystem solver would continue from Hessian information of two different models after a rejected step.
- Of the five `arh_type` options below, only `"arh"` breaks Hessian symmetry; `"symm_arh"`, `"ms_psb"`, `"ms_sp"`, and the default `"ms_sr1"` are symmetric. The factory sets `hess_symm` of the solver settings accordingly.
- Clean up ARH resources by calling `arh_deconstructor`; the final orbitals are the MO coefficients or AO density matrix passed to the factory, which are updated in place.

#### Optional Settings
The ARH factories can be fine-tuned using the following settings:

- **`arh_type`** (string): Specifies which ARH type to use. Options include:
  - `"arh"`: standard ARH method which is not symmetric but fulfills all multisecant conditions,
  - `"symm_arh"`: simple symmetrized version of ARH which does not fulfill the multisecant conditions,
  - `"ms_psb"`: multisecant Powell symmetric Broyden which is only symmetric and fulfills all multisecant conditions for HF and other methods for which the energy is quadratic in the density matrix,
  - `"ms_sp"`: subspace-projected multisecant method,
  - `"ms_sr1"`: multisecant symmetric-rank-1 method, which treats the linear (Coulomb and exact exchange) part of the response exactly as the difference of a Coulomb and an exact-exchange multisecant term, each bounded by its positive semidefinite kernel, and the non-linear (exchange-correlation) part through a separate multisecant system regularized by the multisecant analogue of the SR1 skipping criterion.
- **`verbose`** (integer): Controls the verbosity of output of the ARH routines.
- **`context`** (unlimited polymorphic pointer in Fortran, `void *` in C, absent in Python): Opaque host data, handed back unchanged as the last argument of `evaluate_dm` and of the logging function, as the `context` of the solver settings is to the callback functions of the `solver`. A Python host carries its data in closures instead.

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

This allows the solver to automatically update the Hessian by building and evaluating a surrogate model, so the `hess_x` function is never called. The routine `s_gek_factory` constructs and returns an S-GEK version of the orbital-updating subroutine. This routine requires the following input arguments:

#### Required Arguments

- **`update_orbs`** (subroutine): Original orbital updating subroutine as defined for the `solver` subroutine. The `hess_x` subroutine is not used. The returned orbital updating subroutine calls it with the host context the solver hands to it.
- **`change_reference`** (subroutine): Expresses the parameters, local gradients and gradients of the `n_points` stored points in place relative to a new reference point and returns an integer error code (0 for success, positive integers < 100 for errors). Receives the host context of the S-GEK settings as its last argument.
- **`n_param`** (integer): Specifies the number of parameters to be optimized.
- **`s_gek_settings`** (s_gek_settings_type): Settings object which controls optional arguments as described below.
- **`update_orbs`** (subroutine):  
  Returned S-GEK orbital updating subroutine as defined for the `solver` subroutine.

---

The following Fortran snippet demonstrates how to use the S-GEK interface:

```fortran
use opentrustregion, only: settings_type, update_orbs_type, solver
use otr_s_gek, only: s_gek_factory, s_gek_settings_type, &
                     update_orbs_s_gek_deconstructor

type(settings_type) :: settings
type(s_gek_settings_type) :: s_gek_settings
procedure(update_orbs_type), pointer :: update_orbs_funptr, update_orbs_s_gek_funptr
procedure(obj_func), pointer :: obj_func_funptr
integer(ip) :: n_param, error

! set callback function pointers to existing implementations
update_orbs_funptr => update_orbs
obj_func_funptr => obj_func

! initialize S-GEK settings
call s_gek_settings%init(error)

! override default settings
s_gek_settings%verbose = 1
s_gek_settings%max_points = 30

! get a S-GEK orbital updating subroutine
call s_gek_factory(update_orbs_funptr, n_param, s_gek_settings, update_orbs_s_gek_funptr)

! initialize settings
call settings%init(error)

! call solver
call solver(update_orbs_s_gek_funptr, obj_func_funptr, n_param, error, settings)

! clean up S-GEK objects
call update_orbs_s_gek_deconstructor()
```

---

The following C snippet demonstrates the equivalent usage through the C interface:

```c
#include <string.h>
#include "opentrustregion.h"
#include "opentrustregion_s_gek.h"

c_int n_param;

// set callback function pointers to existing implementations
update_orbs_fp update_orbs_funptr = (void*)update_orbs;
obj_func_fp obj_func_funptr = (void*)obj_func;

// initialize S-GEK settings
s_gek_settings_type s_gek_settings = s_gek_settings_init();

// override default settings
s_gek_settings.verbose = 1;
s_gek_settings.max_points = 30;

// get a S_GEK orbital updating function
update_orbs_fp update_orbs_s_gek_funptr;
s_gek_factory(update_orbs, n_param, s_gek_settings, &update_orbs_s_gek_funptr);

// initialize settings
solver_settings_type settings = solver_settings_init();

// call solver
c_int error =
    solver(update_orbs_s_gek_funptr, obj_func_funptr, n_param, &settings, NULL);

// clean up S-GEK objects
update_orbs_s_gek_deconstructor();
```

---

The following Python snippet demonstrates the equivalent usage through the Python interface:

```python
from pyopentrustregion import SolverSettings, solver
from pyopentrustregion.extensions.s_gek import (
    SGEKSettings, s_gek_factory, update_orbs_s_gek_deconstructor
)

# initialize S-GEK settings
s_gek_settings = SGEKSettings()

# override default settings
s_gek_settings.verbose = 1
s_gek_settings.max_points = 30

# get a S-GEK orbital updating function
update_orbs_s_gek = s_gek_factory(update_orbs, n_param, s_gek_settings)

# initialize settings
settings = SolverSettings()

# call solver
solver(update_orbs_s_gek, obj_func, n_param, settings)

# clean up S-GEK objects
update_orbs_s_gek_deconstructor()
```

---

- Callback functions (`update_orbs`, `obj_func`) point to existing implementations elsewhere in the program.
- `n_param` is assumed to be defined elsewhere.
- S-GEK settings are initialized equivalently to the `solver` and `stability_check` settings; individual settings (here, `verbose` and `max_points`) can then be overridden.
- A S-GEK orbital updating function is obtained by calling `s_gek_factory` with the original `update_orbs` pointer, `n_param`, and `s_gek_settings`.
- The `solver` is then called with the S-GEK-wrapped `update_orbs` function.
- Clean up S-GEK resources by calling `update_orbs_s_gek_deconstructor`.

#### Optional Settings
The S-GEK factory function can be fine-tuned using the following settings:

- **`use_subspace`** (boolean): Specifies whether the surrogate model is constructed in a subspace.
- **`max_points`** (integer): Controls the number of points for the surrogate model and subspace construction.
- **`verbose`** (integer): Controls the verbosity of output during the stability check.
- **`context`** (unlimited polymorphic pointer in Fortran, `void *` in C, absent in Python): Opaque host data, handed back unchanged as the last argument of `change_reference` and of the logging function, as the `context` of the solver settings is to the callback functions of the `solver`. A Python host carries its data in closures instead.

## AI Usage Disclosure

Recent development of OpenTrustRegion has been assisted by AI coding agents. They have been used for code generation, refactoring, writing tests and drafting documentation.

All AI-assisted contributions are reviewed, edited where necessary and validated by the developers before they are merged. The developers make all design decisions, in particular those concerning the underlying methodology, and hold AI-assisted code to the same standards as any other code in the library.
