![License](https://img.shields.io/github/license/eriksen-lab/opentrustregion)
![CI](https://github.com/eriksen-lab/opentrustregion/actions/workflows/main.yml/badge.svg)
[![codecov](https://codecov.io/github/eriksen-lab/opentrustregion/graph/badge.svg?token=NJIM6FDADD)](https://codecov.io/github/eriksen-lab/opentrustregion)
[![DOI](https://zenodo.org/badge/855894860.svg)](https://doi.org/10.5281/zenodo.17142554)

# OpenTrustRegion: A Reusable Library for Second-Order Trust Region Orbital Optimization

This library provides a robust and flexible implementation for second-order trust region orbital optimization, with extensive customization options to suit various use cases.

The following paper documents the theory and implementation of the methodology in OpenTrustRegion, and should be cited in any work using OpenTrustRegion:  

- A Reusable Library for Second-Order Orbital Optimization Using the Trust Region Method\
  Greiner, J.; Høyvik, I.-M.; Lehtola, S.; Eriksen, J. J.\
  [arXiv:2509.13931](https://arxiv.org/abs/2509.13931)

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
| **BUILD_TESTS** | `BOOL` | `ON` | Build the project’s testsuite. |
| **CMAKE_BUILD_TYPE** | `STRING` | `Release` | Choose the build type (`Debug`, `Release`, etc.). |
| **INTEGER_SIZE** | `STRING` | *(auto)* | Set the integer precision to `4` (32-bit) or `8` (64-bit). Required when providing custom BLAS/LAPACK libraries. Otherwise defaults to 32-bit integers and tries to locate compatible BLAS and LAPACK libraries. Falls back to 64-bit integers if 32-bit libraries cannot be found. The resulting library name reflects the chosen integer precision (`libopentrustregion_32.*` or `libopentrustregion_64.*`) |
| **BLAS_LIBRARIES** | `PATH` | *(auto)* | Path(s) to BLAS libraries. If not provided, CMake attempts to locate a suitable BLAS automatically. |
| **LAPACK_LIBRARIES** | `PATH` | *(auto)* | Path(s) to LAPACK libraries. If not provided, CMake attempts to locate a suitable LAPACK automatically. |

## Usage

The optimization process is initiated by calling a `solver` subroutine. This routine requires the following input arguments:

### Required Arguments

- **`update_orbs`** (subroutine):  
  Accepts and applies a variable update (e.g., orbital rotation), updates the internal state, and provides:
  - Objective function value (real)
  - Gradient (real array, written in-place)
  - Hessian diagonal (real array, written in-place)
  - A **`hess_x`** subroutine that performs Hessian-vector products:
    - Accepts a trial vector and writes the result of the Hessian transformation into an output array (real array, written in-place)
    - Returns an integer error code (0 for success, positive integers < 100 for errors)
  - Returns an integer error code (0 for success, positive integers < 100 for errors)
- **`obj_func`** (function):  
  Accepts and applies a variable update (e.g., orbital rotation) and returns:
  - Objective function value (real)
  - An integer error code (0 for success, positive integers < 100 for errors)
- **`n_param`** (integer): Specifies the number of parameters to be optimized.
- **`error`** (integer): An integer code indicating the success or failure of the solver. The error code structure is explained below.

---

The following Fortran snippet demonstrates how to use the `solver` interface:

```fortran
use opentrustregion, only: ip, rp, update_orbs_type, obj_func_type, solver_settings_type, solver

procedure(update_orbs_type), pointer :: update_orbs_funptr
procedure(obj_func_type), pointer :: obj_func_funptr
integer(ip) :: n_param, error
type(solver_settings_type) :: settings

! set callback function pointers to existing implementations
update_orbs_funptr => update_orbs
obj_func_funptr => obj_func

! initialize settings
call settings%init(error)

! override default settings
settings%conv_tol = 1e-6_rp
settings%n_macro = 100
settings%subsystem_solver = "tcg"

! run solver
call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings)
```

- Callback function pointers (`update_orbs_funptr`, `obj_func_funptr`) point to existing implementations elsewhere in the program.
- `n_param` is also assumed to be defined elsewhere.
- Solver settings are initialized using the `init()` method of the derived type, and default settings can be overridden (here, `conv_tol` and `n_macro`).
- Finally, the `solver` is called with the initialized settings and callback functions.

---

The following C snippet demonstrates the equivalent usage through the C interface:

```c
#include <stdint.h>
#include <string.h>
#include "opentrustregion.h"

c_int n_param;

// set callback function pointers to existing implementations
update_orbs_type update_orbs_funptr = (void*)update_orbs;
obj_func_type obj_func_funptr = (void*)obj_func;

// initialize settings
solver_settings_type settings = solver_settings_init();

// override default settings
settings.conv_tol = 1e-6;
settings.n_macro = 100;
strcpy(settings.subsystem_solver, "tcg");

// run solver
c_int error = solver(update_orbs_funptr, obj_func_funptr, n_param, settings);
```

- Callback function pointers (`update_orbs_funptr`, `obj_func_funptr`) point to existing implementations elsewhere in the program.
- `n_param` is also assumed to be defined elsewhere.
- Solver settings are initialized via a small helper function `solver_settings_init()`, which returns a struct with default values. Individual settings (here, `conv_tol` and `n_macro`) can then be overridden.
- Finally, the `solver` is called with the initialized settings and callback functions and directly returns an error code in typical C fashion.

---

The following Python snippet demonstrates the equivalent usage through the Python interface:

```python
from pyopentrustregion import SolverSettings, solver

# initialize settings
settings = SolverSettings()

# override default settings
settings.conv_tol = 1e-6
settings.n_macro = 100
settings.subsystem_solver = "tcg"

# run solver
solver(update_orbs, obj_func, n_param, settings)
```

- Callback functions (`update_orbs`, `obj_func`) are defined elsewhere in the program.
- `n_param` is also assumed to be defined elsewhere.
- Solver settings are initialized via the `SolverSettings` class, which returns an object with default values; individual settings (here, `conv_tol`, and `n_macro`) can then be overridden.
- Finally, the `solver` is called with the initialized settings and callback functions and errors can be caught in pythonic fashion in the form of a `RuntimeException`.

### Optional Settings
The optimization process can be fine-tuned using the following settings:

- **`precond`** (subroutine): Applies a preconditioner to a residual vector. Writes the result in-place into a provided array and returns an integer error code (0 for success, positive integers < 100 for errors).
- **`conv_check`** (function): Returns whether the optimization has converged due to some supplied convergence criterion. Additionally, outputs an integer code indicating the success or failure of the function, positive integers less than 100 represent error conditions.
- **`stability`** (boolean): Determines whether a stability check is performed upon convergence.
- **`hess_symm`** (boolean): Determines whether the supplied Hessian is symmetric. This is sometimes not the case for approximate Hessians.
- **`line_search`** (boolean): Determines whether a line search is performed after every macro iteration.
- **`subsystem_solver`** (string): Specifies which subsystem solver to use. Options include:
  - `"davidson"`: standard Davidson method,
  - `"jacobi-davidson"`: Davidson method with fallback to Jacobi-Davidson if convergence is difficult, or automatically after `jacobi_davidson_start` micro iterations,
  - `"tcg"`: truncated conjugate gradient method.
- **`conv_tol`** (real): Specifies the convergence criterion for the RMS gradient.
- **`n_random_trial_vectors`** (integer): Number of random trial vectors used to initialize the micro iterations.
- **`start_trust_radius`** (real): Initial trust radius.
- **`n_macro`** (integer): Maximum number of macro iterations.
- **`n_micro`** (integer): Maximum number of micro iterations.
- **`jacobi_davidson_start`** (integer): Number of micro iterations after which the subsystem solver switches to the Jacobi-Davidson method.
- **`global_red_factor`** (real): Reduction factor for the residual during micro iterations in the global region.
- **`local_red_factor`** (real): Reduction factor for the residual during micro iterations in the local region.
- **`verbose`** (integer): Controls the verbosity of output during optimization.
- **`seed`** (integer): Seed value for generating random trial vectors.
- **`logger`** (subroutine): Accepts a log message. Logging is otherwise routed to stdout.

## Stability Check
A separate `stability_check` subroutine is available to verify whether the current solution corresponds to a minimum. If not, it returns a boolean indicating instability and optionally, writes the eigenvector corresponding to the negative eigenvalue in-place to the provided memory.

### Required Arguments

- **`h_diag`** (real array): Represents the Hessian diagonal at the current point.
- **`hess_x`** (subroutine): Performs Hessian-vector products at the current point:
  - Accepts a trial vector and writes the result of the Hessian transformation into an output array (real array, written in-place)
  - Returns an integer error code (0 for success, positive integers < 100 for errors)
- **`stable`** (boolean): Returns whether the current point is stable.
- **`error`** (integer): An integer code indicating the success or failure of the solver. The error code structure is explained below.
- **`kappa`** (real array): If the memory is provided and the current point is not stable (as can be checked from return code of `stable`), the descent direction is written in-place in this array.

---

The following Fortran snippet demonstrates how to use the stability check interface:

```fortran
use opentrustregion, only: ip, rp, stability_settings_type, hess_x_type, stability_check

real(rp), allocatable :: h_diag(:), kappa(:)
procedure(hess_x_type), pointer :: hess_x_funptr
integer(ip) :: n_param, error
logical :: stable
type(stability_settings_type) :: settings

! set callback function pointer to existing implementation
hess_x_funptr => hess_x

! initialize settings
call settings%init(error)

! override default settings
settings%conv_tol = 1e-6_rp
settings%n_iter = 100
settings%diag_solver = "jacobi-davidson"

! run stability check
call stability_check(h_diag, hess_x_funptr, n_param, stable, error, settings, kappa=kappa)
```

- `hess_x_funptr` points to an existing Hessian-vector product implementation elsewhere in the program.
- `n_param` is also assumed to be defined elsewhere.
- Stability settings are initialized via the `init()` method of the derived type and can be overridden (here, `conv_tol` and `n_iter`).
- The `stable` logical output receives the result of the stability check.
- The descent direction `kappa` is optional and is only returned if provided.

---

The following C snippet demonstrates the equivalent usage through the C interface:

```c
#include <stdint.h>
#include <string.h>
#include "opentrustregion.h"

c_int n_param;
c_bool stable;

// set callback function pointer to existing implementation
hess_x_type hess_x_funptr = hess_x;

// initialize settings
stability_settings_type settings = stability_settings_init();

// override default settings
settings.conv_tol = 1e-6;
settings.n_iter = 100;
strcpy(settings.diag_solver, "jacobi-davidson");

// pointers to Hessian diagonal and descent direction
double* h_diag;
double* kappa;

// run stability check
c_int = stability_check(h_diag, hess_x_funptr, n_param, &stable, settings, kappa);
```

- `hess_x_funptr` points to an existing Hessian-vector product implementation elsewhere in the program.
- `n_param` and `h_diag` are assumed to be defined elsewhere.
- Stability settings are initialized via a small helper function `stability_settings_init()`, which returns a struct with default values; individual settings (here, `conv_tol` and `n_iter`) can then be overridden.
- The `stable` output receives the result of the stability check which directly returns an error code in typical C fashion.
- The descent direction `kappa` can be defined elsewhere if needed; otherwise, it can be set to `nullptr`.

---

The following Python snippet demonstrates the equivalent usage through the Python interface:

```python
from pyopentrustregion import StabilitySettings, stability_check

# initialize settings
settings = StabilitySettings()

# override default settings
settings.conv_tol = 1e-6
settings.n_iter = 100
settings.diag_solver = "jacobi-davidson"

# Hessian diagonal and descent direction arrays
h_diag = np.asarray(h_diag, dtype=np.float64)
kappa = np.empty(n_param, dtype=np.float64)

# run stability check
stable = stability_check(h_diag, hess_x, n_param, settings, kappa=kappa)
```

- `hess_x` is an existing Hessian-vector product implementation elsewhere in the program.
- `n_param` and `h_diag` are assumed to be defined elsewhere.
- Stability settings are initialized via the `StabilitySettings` class, which returns an object with default values; individual settings (here, `conv_tol`) can then be overridden.
- The `stable` output receives the result of the stability check and errors can be caught in pythonic fashion in the form of a `RuntimeException`.
- The descent direction `kappa` is optional and is only returned if provided.

### Optional Settings
The stability check can be fine-tuned using the following settings:

- **`precond`** (subroutine): Applies a preconditioner to a residual vector. Writes the result in-place into a provided array and returns an integer error code (0 for success, positive integers < 100 for errors).
- **`hess_symm`** (boolean): Determines whether the supplied Hessian is symmetric. This is sometimes not the case for approximate Hessians.
- **`diag_solver`** (string): Specifies which diagonalization solver to use. Options include:
  - `"davidson"`: standard Davidson method,
  - `"jacobi-davidson"`: Davidson method with fallback to Jacobi-Davidson if convergence is difficult, or automatically after `jacobi_davidson_start` micro iterations.
- **`conv_tol`** (real): Convergence criterion for the residual norm.
- **`n_random_trial_vectors`** (integer): Number of random trial vectors used to start the Davidson iterations.
- **`n_iter`** (integer): Maximum number of Davidson iterations.
- **`jacobi_davidson_start`** (integer): Number of micro iterations after which the subsystem solver switches to the Jacobi-Davidson method.
- **`verbose`** (integer): Controls the verbosity of output during the stability check.
- **`seed`** (integer): Seed value for generating random trial vectors.
- **`logger`** (function): Accepts a log message. Logging is otherwise routed to stdout.

## Error Code Structure

The library uses structured integer return codes to indicate whether a function has encountered an error. These codes follow the format **`OOEE`**, where:

- **`OO`** = Origin of the error (which component/function reported the error)
- **`EE`** = Specific error code

### General Rules

- A return code of `0` means success.
- Return codes between `1` and `99` are currently unused.
- All current error codes start from `100` and follow the `OOEE` structure.

### Origins (`OO`)

| Code Prefix (`OO`) | Component           |
|--------------------|---------------------|
| `01`               | `solver`            |
| `02`               | `stability_check`   |
| `11`               | `obj_func`          |
| `12`               | `update_orbs`       |
| `13`               | `hess_x`            |
| `14`               | `precond`           |
| `15`               | `conv_check`        |

### Error Codes (`EE`)

The error field (`EE`) is currently always set to `01`. Future versions may define more specific codes for different failure modes.

### Example Error Codes

| Error Code | Meaning                |
|------------|------------------------|
| `0101`     | Error in `solver`      |
| `1201`     | Error in `update_orbs` |

## Program Interfaces

The PySCF interface is available as an extension hosted at https://github.com/eriksen-lab/pyscf_opentrustregion. To install it, simply add its path to the **`PYSCF_EXT_PATH`** environment variable:
```sh
export PYSCF_EXT_PATH=path/to/pyscf_opentrustregion
```
Usage examples can be found in the **`examples`** directory of the PySCF interface repository.

The interface supports Hartree–Fock and DFT calculations via the **`mf_to_otr`** function, which wraps PySCF **`HF`** and **`KS`** objects into their OpenTrustRegion counterparts. Similarly, localization methods are available through the **`BoysOTR`**, **`PipekMezeyOTR`**, and **`EdmistonRuedenbergOTR`** classes, and state-specific CASSCF calculations are supported via the **`casscf_to_otr`** function applied to a PySCF **`CASSCF`** object. All returned objects are fully compatible with the original PySCF classes and can be used interchangeably.

Optional settings can be adjusted by modifying object attributes directly. Orbital optimization and internal stability analysis are performed using the **`kernel`** and **`stability_check`** member functions, respectively.

## Optional Extensions

OpenTrustRegion (OTR) supports additional optional modules that can be enabled during installation. Two main extensions are currently available: **quasi-Newton** and **augmented Roothaan–Hall (ARH)**.

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

This allows the solver to automatically update the Hessian using a quasi-Newton scheme, so the `hess_x` function is never called. The routine `update_orbs_qn_factory` constructs and returns a quasi-Newton version of the orbital-updating subroutine. This routine requires the following input arguments:

#### Required Arguments

- **`update_orbs`** (subroutine): Original orbital updating subroutine as defined for the `solver` subroutine. The `hess_x` subroutine is not used.
- **`n_param`** (integer): Specifies the number of parameters to be optimized.
- **`qn_settings`** (qn_settings_type): Settings object which controls optional arguments as described below.
- **`update_orbs`** (subroutine):  
  Returned quasi-Newton orbital updating subroutine as defined for the `solver` subroutine.

---

The following Fortran snippet demonstrates how to use the quasi-Newton interface:

```fortran
use opentrustregion, only: settings_type, update_orbs_type, solver
use opentrustregion_quasi_newton, only: update_orbs_qn_factory, qn_settings_type, &
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
call update_orbs_qn_factory(update_orbs_funptr, n_param, qn_settings, update_orbs_qn_funptr)

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
update_orbs_qn_factory(update_orbs, n_param, qn_settings, &update_orbs_qn_funptr);

// initialize settings
solver_settings_type settings = solver_settings_init();

// call solver
c_int error = solver(update_orbs_qn_funptr, obj_func_funptr, n_param, settings);

// clean up quasi-Newton objects
update_orbs_qn_deconstructor();
```

---

The following Python snippet demonstrates the equivalent usage through the Python interface:

```python
from pyopentrustregion import SolverSettings, solver
from pyopentrustregion.extensions.quasi_newton import (
    QNSettings, update_orbs_quasi_factory, update_orbs_quasi_deconstructor
)

# initialize quasi-Newton settings
quasi_settings = QNSettings()

# override default settings
quasi_settings.verbose = 1
quasi_settings.hess_update_scheme = "sr1"

# get a quasi-Newton orbital updating function
update_orbs_qn = update_orbs_quasi_factory(update_orbs, n_param, quasi_settings)

# initialize settings
settings = SolverSettings()

# call solver
solver(update_orbs_qn, obj_func, n_param, settings)

# clean up quasi-Newton objects
update_orbs_quasi_deconstructor()
```

---

- Callback functions (`update_orbs`, `obj_func`) point to existing implementations elsewhere in the program.
- `n_param` is assumed to be defined elsewhere.
- Quasi-Newton settings are initialized equivalently to the `solver` and `stability_check` settings; individual settings (here, `verbose` and `hess_update_scheme`) can then be overridden.
- A quasi-Newton orbital updating function is obtained by calling `update_orbs_qn_factory` with the original `update_orbs` pointer, `n_param`, and `qn_settings`.
- The `solver` is then called with the quasi-Newton-wrapped `update_orbs` function.
- Clean up quasi-Newton resources by calling `update_orbs_qn_deconstructor`.

#### Optional Settings
The quasi-Newtn factory function can be fine-tuned using the following settings:

- **`hess_update_scheme`** (string): Specifies which quasi-Newton updating scheme to use. Options include:
  - `"sr1"`: symmetric rank-1 update,
  - `"bfgs"`: Broyden-Fletcher-Goldfarb-Shanno update,
- **`verbose`** (integer): Controls the verbosity of output during the stability check.

### Augmented Roothaan–Hall (ARH)

This extension provides access to the augmented Roothaan-Hall method for RHF in the orthogonalized atomic orbital (AO) basis.

#### Installation

Enable the extension at build time using CMake:

```sh
cmake -DENABLE_ARH=ON ..
```

or, when installing the Python package:

```sh
CMAKE_FLAGS='-DENABLE_ARH=ON' pip install .
```

This exposes an `arh_factory` function that prepares ARH-specific callbacks for energy, orbital updates, and preconditioners.

#### Usage

The routine `arh_factory` constructs and returns ARH versions of energy, orbital updating, and preconditioning functions. This routine requires the following input arguments:

#### Required Arguments

- **`dm_ao`** (real array): Represents the starting AO density matrix with dimension (`n_ao`, `n_ao`).
- **`ao_overlap`** (real array): Represents the AO overlap matrix with dimension (`n_ao`, `n_ao`).
- **`n_ao`** (integer): Specifies the number of AOs.
- **`get_energy`** (function):  
  Accepts an AO density matrix and returns:
  - Energy value (real)
  - An integer error code (0 for success, positive integers < 100 for errors)
- **`get_fock`** (subroutine):  
  Accepts an AO density matrix and returns:
  - Energy value (real)
  - Fock matrix (real array, written in-place)
  - An integer error code (0 for success, positive integers < 100 for errors)
- **`obj_func_arh`** (subroutine): Returned ARH objective function as defined for the `solver` subroutine.
- **`update_orbs_arh`** (subroutine): Returned ARH orbital updating subroutine as defined for the `solver` subroutine.
- **`precond_arh`** (subroutine): Returned ARH preconditioning subroutine as defined for the `solver` subroutine.
- **`error`** (integer): An integer code indicating the success or failure of the solver. The error code structure is explained below.
- **`arh_settings`** (arh_settings_type): Settings object which controls optional arguments as described below.


---

The following Fortran snippet demonstrates how to use the ARH interface:

```fortran
use opentrustregion, only: settings_type, solver, obj_func_type, update_orbs_type, &
                           precond_type
use opentrustregion_arh, only: arh_factory, arh_settings_type, init_arh_settings, &
                               arh_deconstructor

type(settings_type) :: settings
type(arh_settings_type) :: arh_settings
procedure(obj_func_type), pointer :: obj_func_arh_funptr
procedure(update_orbs_type), pointer :: update_orbs_arh_funptr
procedure(precond_type), pointer :: precond_arh_funptr
integer(ip) :: n_ao, n_param, error
real(rp), allocatable :: dm_ao(:, :), ao_overlap(:, :)

! set callback function pointers to existing implementations
get_energy_funptr => get_energy
get_fock_funptr => get_fock

! initialize ARH settings
call init_arh_settings(arh_settings)

! override default settings
arh_settings%verbose = 1

! get ARH routines
call arh_factory(dm_ao, ao_overlap, n_ao, get_energy, get_fock, obj_func_arh_funptr, &
                 update_orbs_arh_funptr, precond_arh_funptr, error, arh_settings)

! initialize settings
call settings%init(error)

! override settings
settings%precond => precond_arh_funptr
settings%hess_symm = .false.

! set number of parameters
n_param = n_ao * (n_ao - 1) / 2

! call solver
call solver(obj_func_arh_funptr, update_orbs_arh_funptr, n_param, error, settings)

! deallocate ARH objects and get final density matrix
call arh_deconstructor(dm_ao, error)
```

---

The following C snippet demonstrates equivalent usage through the C interface:

```C
#include "opentrustregion.h"
#include "opentrustregion_arh.h"

c_int n_ao, n_param;
c_real dm_ao[n_ao][n_ao], ao_overlap[n_ao][n_ao]

// set callback function pointers to existing implementations
get_energy_fp get_energy_funptr = (void*)get_energy;
get_fock_fp get_fock_funptr = (void*)get_fock;

// initialize ARH settings
arh_settings_type arh_settings = arh_settings_init();

// override default settings
arh_settings.verbose = 1;

// get callback functions
obj_func_fp obj_func_arh_funptr;
update_orbs_fp update_orbs_arh_funptr;
precond_fp precond_arh_funptr;
c_int error = arh_factory(dm_ao, 
                          ao_overlap, 
                          n_ao,
                          get_energy, 
                          get_fock,
                          &obj_func_arh_funptr,
                          &update_orbs_arh_funptr,
                          &precond_arh_funptr
                          arh_settings);

// initialize settings
solver_settings_type settings = solver_settings_init();

// override settings
settings.precond = precond_arh_funptr
settings.hess_symm = false;

// set number of parameters
n_param = n_ao * (n_ao - 1) / 2

// call solver
error = solver(obj_func_arh_funptr, update_orbs_arh_funptr, n_param, settings);

// deallocate ARH objects and get final density matrix
error = arh_deconstructor(dm_ao);
```

---

The following Python snippet demonstrates the equivalent usage through the Python interface:

```python
from pyopentrustregion import SolverSettings, solver
from pyopentrustregion import ARHSettings, arh_factory, arh_deconstructor

# Hessian diagonal and descent direction arrays
dm_ao = np.asarray(dm_ao, dtype=np.float64)
ao_overlap = np.asarray(ao_overlap, dtype=np.float64)

# initialize ARH settings
arh_settings = ARHSettings()

# override default settings
arh_settings.verbose = 1

# get callback functions
obj_func_arh, update_orbs_arh, precond_arh = arh_factory(
    dm_ao, ao_overlap, n_ao, get_energy, get_fock, arh_settings
)

# initialize settings
settings = SolverSettings()

# override settings
settings.precond = precond_arh
settings.hess_symm = False

# set number of parameters
n_param = n_ao * (n_ao - 1) / 2

# call solver
solver(obj_func_arh, update_orbs_arh, n_param, settings)

# deallocate ARH objects and get final density matrix
arh_deconstructor(dm_ao)
```

---

- `dm_ao`, `ao_overlap` and `n_ao` are assumed to be prepared elsewhere.
- `get_energy` and `get_fock` are callback procedures provided elsewhere.
- ARH settings are initialized equivalently to the `solver` and `stability_check` settings; individual settings (here, `verbose`) can then be overridden.
- `arh_factory` returns three procedures: `obj_func`, `update_orbs`, and `precond`, which are passed to the normal `solver`.
- ARH breaks Hessian symmetry; when using ARH set `hess_symm` to `false` before calling `solver`.
- Clean up ARH resources and get final AO density matrix by calling `arh_deconstructor`.

#### Optional Settings
The ARH factory function can be fine-tuned using the following settings:

- **`verbose`** (integer): Controls the verbosity of output during the stability check.
