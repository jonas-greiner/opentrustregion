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

The optimization process is initiated by calling a solver function. This function requires the following input arguments:

### Required Arguments

- **`update_orbs`**: A function that accepts the variable change (i.e., the orbital rotation), updates the variable (the orbitals), and outputs:
  - The objective function value
  - The gradient
  - The Hessian diagonal
  - A `hess_x` function that performs a Hessian linear transformation for a trial function

- **`obj_func`**: A function that accepts the variable change and returns the objective function value.
- **`n_param`**: An integer specifying the number of parameters to be optimized.

### Optional Arguments
The optimization process can be fine-tuned using the following optional arguments:

- **`stability`** (boolean): Determines whether a stability check is performed upon convergence.
- **`line_search`** (boolean): Determines whether a line search is performed after every macro iteration.
- **`conv_tol`** (real): Specifies the convergence criterion for the RMS gradient.
- **`n_random_trial_vectors`** (integer): Number of random trial vectors used to initialize the micro iterations.
- **`start_trust_radius`** (real): Initial trust radius.
- **`n_macro`** (integer): Maximum number of macro iterations.
- **`n_micro`** (integer): Maximum number of micro iterations.
- **`global_red_factor`** (real): Reduction factor for the residual during micro iterations in the global region.
- **`local_red_factor`** (real): Reduction factor for the residual during micro iterations in the local region.
- **`verbose`** (integer): Controls the verbosity of output during optimization.
- **`seed`** (integer): Seed value for generating random trial vectors.

## Stability Check
A separate `stability_check` function is available to verify whether the current solution corresponds to a minimum. If not, it returns a boolean indicating instability and an additional direction along the eigenvector corresponding to the negative eigenvalue.

### Required Arguments

- **`grad`**: Array representing the gradient at the current point.
- **`h_diag`**: Array representing the Hessian diagonal at the current point.
- **`hess_x`**: A function that performs a Hessian linear transformation of a trial vector at the current point.

### Optional Arguments

- **`conv_tol`** (real): Convergence criterion for the residual norm.
- **`n_random_trial_vectors`** (integer): Number of random trial vectors used to start the Davidson iterations.
- **`n_iter`** (integer): Maximum number of Davidson iterations.
- **`verbose`** (integer): Controls the verbosity of output during the stability check.

---
Both the solver and stability check functions can be directly accessed from Fortran, C, or Python using the same arguments but within the appropriate language.

### Program Interfaces

The PySCF interface is available as an extension hosted at https://github.com/eriksen-lab/pyscf_opentrustregion. To install it, simply add its path to the **`PYSCF_EXT_PATH`** environment variable:
```sh
export PYSCF_EXT_PATH=path/to/pyscf_opentrustregion
```
Usage examples can be found in the **`examples`** directory of the PySCF interface repository.

The interface supports Hartree–Fock and DFT calculations via the **`mf_to_otr`** function, which wraps PySCF **`HF`** and **`KS`** objects into their OpenTrustRegion counterparts. Similarly, localization methods are available through the **`BoysOTR`**, **`PipekMezeyOTR`**, and **`EdmistonRuedenbergOTR`** classes, and state-specific CASSCF calculations are supported via the **`casscf_to_otr`** function applied to a PySCF **`CASSCF`** object. All returned objects are fully compatible with the original PySCF classes and can be used interchangeably.

Optional settings can be adjusted by modifying object attributes directly. Orbital optimization and internal stability analysis are performed using the **`kernel`** and **`stability_check`** member functions, respectively.
