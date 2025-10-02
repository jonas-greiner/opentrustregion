![License](https://img.shields.io/github/license/eriksen-lab/opentrustregion)
![CI](https://github.com/eriksen-lab/opentrustregion/actions/workflows/main.yml/badge.svg)
[![codecov](https://codecov.io/github/eriksen-lab/opentrustregion/graph/badge.svg?token=NJIM6FDADD)](https://codecov.io/github/eriksen-lab/opentrustregion)
[![DOI](https://zenodo.org/badge/855894860.svg)](https://doi.org/10.5281/zenodo.17142554)

# OpenTrustRegion: A Reusable Library for Second-Order Trust Region Orbital Optimization

This library provides a robust and flexible implementation for second-order trust region orbital optimization, with extensive customization options to suit various use cases.

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

## Usage

The optimization process is initiated by calling a `solver` subroutine. This routine requires the following input arguments:

### Required Arguments

- **`update_orbs`** (subroutine):  
  Accepts and applies a variable update (e.g., orbital rotation), updates the internal state, and provides:
  - Objective function value (real)
  - Gradient (real array, written in-place)
  - Hessian diagonal (real array, written in-place)
  - A **`hess_x`** subroutine that performs Hessian-vector products:
    - Accepts a trial vector and writes the result of the Hessian transformation into an output array
    - Returns an integer error code (0 for success, positive integers < 100 for recoverable errors)
  - Returns an integer error code (0 for success, positive integers < 100 for recoverable errors)
- **`obj_func`** (function):  
  Accepts and applies a variable update (e.g., orbital rotation) and returns:
  - Objective function value (real)
  - An integer error code (0 for success, positive integers < 100 for recoverable errors)
- **`n_param`** (integer): Specifies the number of parameters to be optimized.
- **`error`** (integer): An integer code indicating the success or failure of the solver. The error code structure is explained below.

### Optional Arguments
The optimization process can be fine-tuned using the following optional arguments:

- **`precond`** (subroutine): Applies a preconditioner to a residual vector. Writes the result in-place into a provided array and returns an integer error code (0 for success, positive integers < 100 for recoverable errors).
- **`conv_check`** (function): Returns whether the optimization has converged due to some supplied convergence criterion. Additionally, outputs an integer code indicating the success or failure of the function, positive integers less than 100 represent error conditions.
- **`stability`** (boolean): Determines whether a stability check is performed upon convergence.
- **`line_search`** (boolean): Determines whether a line search is performed after every macro iteration.
- **`davidson`** (boolean): Determines whether level-shifted augmented Hessian with Davidson or truncated conjugate gradient is utilized to solve the trust-region subsystem.
- **`jacobi_davidson`** (boolean): Determines whether Jacobi-Davidson is performed whenever difficult convergence is encountered for Davidson iterations.
- **`prefer_jacobi_davidson`** (boolean): Determines whether Jacobi-Davidson should be preferred over shrinking of the trust region whenever difficult convergence is encountered for Davidson iterations.
- **`conv_tol`** (real): Specifies the convergence criterion for the RMS gradient.
- **`n_random_trial_vectors`** (integer): Number of random trial vectors used to initialize the micro iterations.
- **`start_trust_radius`** (real): Initial trust radius.
- **`n_macro`** (integer): Maximum number of macro iterations.
- **`n_micro`** (integer): Maximum number of micro iterations.
- **`global_red_factor`** (real): Reduction factor for the residual during micro iterations in the global region.
- **`local_red_factor`** (real): Reduction factor for the residual during micro iterations in the local region.
- **`verbose`** (integer): Controls the verbosity of output during optimization.
- **`seed`** (integer): Seed value for generating random trial vectors.
- **`logger`** (subroutine): Accepts a log message. Logging is otherwise routed to stdout.

## Stability Check
A separate `stability_check` subroutine is available to verify whether the current solution corresponds to a minimum. If not, it returns a boolean indicating instability and optionally an additional direction along the eigenvector corresponding to the negative eigenvalue.

### Required Arguments

- **`h_diag`** (real array): Represents the Hessian diagonal at the current point.
- **`hess_x`** (subroutine): Performs Hessian-vector products at the current point:
  - Accepts a trial vector and writes the result of the Hessian transformation into an output array
  - Returns an integer error code (0 for success, positive integers < 100 for recoverable errors)
- **`stable`** (boolean): Returns whether the current point is stable.
- **`error`** (integer): An integer code indicating the success or failure of the solver. The error code structure is explained below.

### Optional Arguments

- **`kappa`** (real array): If the memory is provided, returns descent direction if current point is not stable otherwise returns zeros.
- **`precond`** (subroutine): Applies a preconditioner to a residual vector. Writes the result in-place into a provided array and returns an integer error code (0 for success, positive integers < 100 for recoverable errors).
- **`jacobi_davidson`** (boolean): Determines whether Jacobi-Davidson is performed whenever difficult convergence is encountered for Davidson iterations.
- **`conv_tol`** (real): Convergence criterion for the residual norm.
- **`n_random_trial_vectors`** (integer): Number of random trial vectors used to start the Davidson iterations.
- **`n_iter`** (integer): Maximum number of Davidson iterations.
- **`verbose`** (integer): Controls the verbosity of output during the stability check.
- **`logger`** (function): Accepts a log message. Logging is otherwise routed to stdout.

---
Both the `solver` and `stability_check` functions can be directly accessed from Fortran, C, or Python using the same arguments but within the appropriate language.

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
