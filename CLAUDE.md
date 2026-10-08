# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

OpenTrustRegion is a Fortran library implementing a second-order trust-region optimizer, exposed via Fortran, C, and Python (ctypes) interfaces — the same compiled shared/static library is consumed by all three. The core solver (`src/`) is agnostic to what is being optimized: it operates purely on abstract parameter/gradient/Hessian-vector-product callbacks. Domain-specific functionality is layered on top as opt-in extensions under `extensions/<name>/`. Three of these, OAO, MO and ARH, are chemistry-specific — they build density-matrix and orbital-optimization machinery on top of the core solver. OAO and MO each provide an orbital object, the orbitals parameterized in the orthogonalized AO basis and in the MO basis; ARH parameterizes the orbitals in either basis, on top of the corresponding extension, so `ENABLE_ARH` always builds OAO and MO. The MO numerics never depend on OAO's: what both bases share lives in `extensions/common` (ARH's C and Python interfaces build on OAO's or MO's for the respective basis). The other two, quasi-Newton and S-GEK are, like the core, agnostic to the underlying problem; they just add alternative Hessian-approximation strategies.

## General conventions

Apply equally to `src/` and every `extensions/<name>/`.

**Formatting.** Fortran lines max 88 columns. Layout, spacing, keyword spelling, string literals, comments and blank lines of Fortran code are defined by `tools/f90_layout.py`, whose docstring is the authoritative rule list; don't restate its rules here or apply them by hand. Run `python3 tools/f90_layout.py` after every Fortran edit: by default it checks only statements touched relative to `HEAD`, `--fix` applies the fixes, and `--all` checks whole files (with `--fix` it requires explicit paths, e.g. `git ls-files '*.f90' | xargs python3 tools/f90_layout.py --all --fix`). CI runs `--all` on every tracked `.f90` file. The default mode also checks every untracked, non-ignored `.f90` file in full, so a bare `--fix` rewrites scratch copies of Fortran files left in the tree too; move them out first. Statements containing a comment, a `;` or a continuation line starting with `&` are never reflowed (their spacing and indentation are still fixed), so shorten an overlong one by hand. Every procedure opens with a `!`-delimited comment block describing what it does; each logical step inside gets a lowercase `!` comment. Unit tests: assume success (`test_<name> = .true.`), then one `if (...) then` / `write(stderr, *) "test_<name> failed: ..."` / `test_<name> = .false.` block per assertion, no early returns except where a later assertion would crash.

**Python formatting.** All Python (`pyopentrustregion/`, its `extensions/`, `setup.py`, `tools/`) is `black`-formatted, default settings. Run `black pyopentrustregion setup.py tools` before considering Python changes done.

**C formatting.** All C headers and test sources (`include/`, `tests/*.c`) are `clang-format`-formatted per the repo-root `.clang-format` (LLVM style, 88-column limit to match the Fortran convention above). Run `clang-format -i` on touched C/H files before considering C changes done.

**Clarity over performance.** The library is not on the hot path of a quantum-chemistry calculation — the host program's integral transforms and Hessian-vector products dominate. Prefer short, obviously-correct code over fast code. Don't propose performance refactors (buffer growth, pooling, micro-optimizations) without evidence the affected code is hot for a real workload. This does not extend to the OAO, MO and ARH extensions: they run in every solver iteration and ARH's approximate Hessian replaces the host's exact Hessian-vector products, so their per-iteration work competes with what they replace. Don't add avoidable `O(n^3)` work there for simplicity, such as transforming the same matrix twice or rebuilding a matrix an object already holds.

**Norms: BLAS in production, intrinsics in tests.** Production Fortran (`src/`, `extensions/*/src/`) uses BLAS (`dnrm2`, `ddot`) for norms/dot products. Test code (`tests/`, `extensions/*/tests/`) uses intrinsics (`norm2`, `dot_product`, `matmul`, `sum`, `transpose`) instead. BLAS/LAPACK is allowed in tests only where no intrinsic exists (`dsyev`, `dgeev`, `zheev`) — don't hand-roll linear algebra to avoid the dependency.

### Unit tests must not depend on untested code

A unit test may only call the routine it is testing — never another production routine to build its inputs or expected values, or the test silently inherits that routine's correctness.

- Construct inputs directly in the test, randomly or with LAPACK, not by calling the production routine that would normally produce them.
- When a routine under test internally calls another module's routines, reimplement those as `ref_*` helpers so the expected value is computed independently — e.g. `ref_unpack_asymm` / `ref_pack_asymm` / `ref_project_asymm` / `ref_hess_x_oao` in `extensions/oao/tests/oao_unit_tests.f90`, used by both OAO's and ARH's tests.
- Where a routine merely delegates (e.g. an `update_orbs_*_callback` returning the gradient another routine produced), assert the *contract* — the returned arrays equal what the routine stored — rather than re-deriving numbers a different routine's own test already covers.

### Test code layout

Classify test code by **role** first, then place it at the most common **level** that needs it.

Roles, one file each:

| File | Role |
|---|---|
| `test_reference.f90` | Tolerances, reference values, `ref_*` reimplementations computed independently of the routine under test. |
| `<name>_unit_tests.f90` | The Fortran-level `test_*` functions, plus mocks/fixtures needed only to drive them (e.g. a stand-in `evaluate_dm`). |
| `<name>_mock.f90` | Mocks of the module's *own* production routines that are bridged to the C interface (factories, deconstructors), so a caller (typically `<name>_c_interface_unit_tests.f90`) can verify invocation without running the real logic. |
| `<name>_c_interface_unit_tests.f90` | Tests for the `bind(C)` wrapper layer; defines its own local `bind(C)` mocks rather than using `<name>_c_interface_mock.f90`. |
| `<name>_c_interface_mock.f90` | `bind(C)`-signature mocks used exclusively by the Python interface tests, dynamically loaded via `ctypes` from `libotrtestsuite`. |

Level: within `tests/`, `extensions/common/tests/`, and a specific extension's `tests/`, code lives at the most common level shared by everything that *actually* needs it — verified against real usage across every extension, not assumed from whichever extension motivated writing it.
- Needed by every extension → `extensions/common/tests/` (e.g. `common_mock.f90`'s `mock_update_orbs`/`mock_hess_x`, since `update_orbs_type`/`hess_x_type` and their C bridge are defined once in the core).
- Needed by some → the most foundational one among them. ARH depends on OAO and MO, so an OAO-specific helper shared between OAO's and ARH's tests (e.g. the `ref_*` twins of OAO's numerics) lives in `oao_unit_tests.f90`, and an MO-specific one shared between MO's and ARH's tests (e.g. `setup_mo_object` or the `ref_*` twins of the MO numerics) in `mo_unit_tests.f90`. OAO and MO do not depend on each other, so anything their tests share (shared dimensions, generators, mock-call recording, the reference value of the shared settings) lives in `common_test_reference.f90`/`common_unit_tests.f90`, and the callback mocks their C-interface tests share (`mock_obj_func`, `mock_precond`, ...) in `common_mock.f90`, together with those every extension needs.
- Needed by one → that extension's own files. The mock logger and `setup_settings`, needed by nearly every Fortran-level unit test, live in `tests/opentrustregion_unit_tests.f90`.

This recurses within a role: a helper belongs in `<name>_test_reference.f90` vs `<name>_unit_tests.f90` based on whether something *at that level* calls it — not which extension's tests happen to use it (`ref_unpack_asymm` etc. live in `oao_unit_tests.f90`, since only unit tests call them).

**`ref_*` naming.** The prefix marks a role — an independent reimplementation of what a production routine computes — not a location; it travels with the function wherever the placement rule moves it, and it is decided by what the helper computes, not by whether a call site compares against it or feeds it in as an input. A fixture that builds a *plausible* input (randomly, via LAPACK, or by reassembling a routine's own output) does not take the prefix, even in the same file. Either kind may keep an atypical name purely to dodge a collision with a production routine or local variable in scope at its call site (`ref_unpack_asymm` keeps its prefix because `unpack_asymm` is in scope) — check for a collision before insisting on the "cleaner" name.

**Module-level `use` in test files is restricted.** In `<name>_unit_tests.f90` and siblings, module-level `use` statements above `contains` are limited to `rp`, `ip`, `kw_len`, `stderr`, `stdout`, `tol`, `tol_c`, and intrinsic module bindings (`iso_c_binding` and the like) — check other modules for the exact set before adding to it. The one exception is the parent of a test-local type extension, since a type with type-bound procedures can only be declared in a module's specification part (`common_unit_tests.f90` imports `orbital_basis_type` for its `mock_orbital_basis_type`). Everything else (mock callbacks, `setup_settings`, shared dimension parameters, `ref_*` helpers) is imported locally inside the specific `test_*` function that needs it, even if several functions in the file need the same symbol. A module-level `use` beyond this set pollutes every procedure in the file and hides which test actually relies on what.

**Prefer shared dimension-parameter fixtures over local literals**, unless it complicates the test a lot. The test reference modules define canonical dimensions (`otr_common_test_reference` those both orbital bases use, OAO's and MO's reference modules the derived and MO-specific ones) that most unit tests should import rather than redeclaring.
- Every dimension is defined once, without a basis suffix (`n_ao`, not `n_ao_mo`), and one quantity has one name (the number of occupied orbitals of a channel, which is also the rank of its density matrix, is `n_occ` everywhere, never `n_electrons`).
- A test whose local variable already carries the name renames the import on its `use` line only, with the suffix `_ref` (`n_occ_ref => n_occ`), never the prefix `ref_`, which marks independent reimplementations.
- The shared `n_param` is only right for a test that works at the shared `n_particle`; a test pinning its own `n_particle` (closed-shell-only, or walking both shells in one body) needs its own local value.
- Exception: a test whose expected values are hand-derived for a specific small/structured case (a closed-form 2×2 rotation, 3×3 matrices with hand-computed results, a matrix crafted to trigger a rank-detection path) may keep local literals and declare its own dimensions, since forcing the shared size would mean re-deriving the numbers. When correctness is checked via algebraic properties rather than exact values, use the shared fixture.

**Prefer generated data over hand-typed literals** (`generate_random_density_matrix` / `generate_random_symm_matrix` / `call random_number`), unless the specific values are load-bearing. A hand-typed `reshape([...])` is justified only when the test checks a specific hand-derivable closed-form target depending on those exact numbers (a weighted-symmetrization blend, a dependency-screening path, a regularization threshold crossed at a specific eigenvalue). When every assertion is relational (a copy, a fixed scaling, an algebraic identity, or a value the test recomputes independently from the same input) the specific values never mattered, and random generation exercises the same code path while being harder to accidentally pass with a bug. Before hand-typing a matrix meant to be "just some valid state" (often signaled by the test calling it "arbitrary"), check whether random generation works instead. Where a matrix must stay block-diagonal/diagonal/structurally constrained to match what production produces, keep that structure but randomize within it.

### Shared architectural patterns

**Settings types.** `solver_settings_type` / `stability_settings_type` (core) and each extension's own settings type (`s_gek_settings_type`, `qn_settings_type`, and `oao_settings_type`, `mo_settings_type` and `arh_settings_type`, which extend `otr_common`'s `orbital_settings_type`, the settings the orbital objects hold, since ARH passes its settings to the setup of either orbital object) all extend the abstract `settings_type`, with an `init` type-bound procedure and an `initialized` flag. The owning routine checks `settings%initialized` on entry and calls `init` if needed; C/Python wrappers pre-populate defaults via `*_init()`/`__init__` so the flag is `.true.` by the time Fortran sees it. Default values are defined once in a `default_*_settings` parameter and replicated in the C `*_init()` function and the Python `Structure` defaults — keep them in sync.

**The Python → C → Fortran call chain** is the same shape everywhere, illustrated by the core solver: `solver(...)` in `python_interface.py` builds a `SolverSettings` ctypes struct, wraps Python callbacks with `CFUNCTYPE`, and calls into `libopentrustregion`'s C entry point (`src/c_interface.f90`), which stores the C function pointers in module-level `procedure(..)_before_wrapping` slots and invokes `standard_solver` (the real `solver`) with Fortran-shaped wrappers that call back through those pointers. Every extension's `<name>_c_interface.f90` follows the identical shape. **Consequence:** none of these C interfaces are re-entrant across threads — module-level pointer state is shared. A signature change must keep all layers (Fortran abstract interface → C wrapper + abstract interface → C typedef → Python `CFUNCTYPE` + `Structure`) in lockstep.

## Build & test

```sh
# Fortran/C build
mkdir build && cd build
cmake ..              # add -DBUILD_SHARED_LIBS=ON for shared
cmake --build .

# Python install (invokes CMake under the hood via setup.py)
pip install .
pip install -e .       # editable
```

```sh
# full suite (Python driver calling into libotrtestsuite: Fortran unit + system tests)
python3 -m pyopentrustregion.testsuite                  # from an installed/editable build
python3 pyopentrustregion/testsuite.py                  # from the source tree against ./build

# single test class or method (stdlib unittest)
python3 -m unittest pyopentrustregion.testsuite.SystemTests
python3 -m unittest pyopentrustregion.testsuite.OpenTrustRegionTests.test_solver
```

The Python suite runs Fortran- and C-side tests (via symbols dynamically loaded from `libotrtestsuite`) and pure-Python wrapper tests. System tests need `pyopentrustregion/test_data/*.bin`.

**Where to add coverage for a new setting.** As assertions inside the existing unit test for the routine the setting affects (`tests/opentrustregion_unit_tests.f90`, or the analogous `<name>_unit_tests.f90`), not as a new system test — system tests exercise real chemistry data end-to-end, and one per settings combination would explode combinatorially. When a setting can push a routine down a materially different code path (interior vs. boundary-crossing branch in a trust-region solver), exercise both — a check that only hits the "easy" branch (e.g. always starting near a minimum) can pass while a branch reached only near a saddle point stays untested.

### Registering a new Fortran unit test

A `logical(c_bool) function test_<routine>() bind(C)` is only reachable once wired into all four:

1. The source file in the CMake test source list (`OPENTRUSTREGION_TESTS`, inside the matching `if(ENABLE_<EXT>)` block for an extension test, unconditionally for a core test).
2. The `fortran_tests["<ext>_tests"]` name list in `pyopentrustregion/extensions/<ext>/tests.py` (or the core equivalent in `pyopentrustregion/tests.py`), alphabetically ordered, without the `test_` prefix.
3. An `@add_tests class <Ext>Tests` in that file, whose `tests` attribute is that list.
4. The class import in `pyopentrustregion/testsuite.py`.

A name in the Python list with no matching Fortran symbol fails at import with `AttributeError: dlsym(...): symbol not found`.

**Exactly one registered test per production routine** — never `test_<routine>_<case_a>`, `test_<routine>_<case_b>`. Cover multiple configurations (each `arh_type`, closed vs. open shell, each basis) as cases *inside* that one test: loop where setup can be shared, or call a per-case `check_*` helper that is deliberately not `bind(C)` and not registered. Keep the case name in the failure message (`"test_<routine> failed for <case>: ..."`) since the harness only reports the registered name. Worked example: `test_inv_hess_x_arh` — four cases (MO/OAO basis × closed/open shell), one registered entry, driven by its internal procedure `check_inv_hess_x_arh_case`, parameterized over the basis and the shell. Loop only over the configurations the routine itself distinguishes: `inv_hess_x_arh` depends on the ARH type only through the model it inverts, which its test sets up directly, so its test does not loop over the ARH types.

Because a routine's cases live behind one entry, a silently-skipped case is invisible in suite output — have the driving test run every case unconditionally and `and` the results together, rather than returning early on the first failure.

**Helper procedures in test files.**
- A `check_*` (or other non-`bind(C)`) helper earns its existence by eliminating real duplication, never merely to keep one case's code out of the registered test's body — five near-identical one-case functions just relocate the duplication; write one function parameterized over the case, or a loop with a shared body.
- A helper that would be called once without arguments, only reading and writing its host's variables, is part of the host's body and stays inline. A single call justifies a separate procedure only for a mock callback, which has to be a procedure, or for a function of its arguments used in an expression, such as a `ref_*` reimplementation.
- A helper only one procedure uses is an internal procedure of that procedure, after its own `contains`, unless it is one of a family of `ref_*` twins whose other members live at module level (`ref_project_symm` next to `ref_project_asymm`). An internal procedure cannot contain another, so helpers calling each other are nested side by side in their common caller. Mocks in `<name>_mock.f90`/`<name>_c_interface_mock.f90` stay module procedures however many callers they have, since those files exist to export them.
- An internal procedure uses its host's variables and imports directly instead of redeclaring them or taking them as arguments that every call fills with the host's same-named variable. What must stay its own (loop indices, which would otherwise overwrite a host loop calling it, and dummies receiving a slice or a different array) gets a name the host does not use, so that it never masks a host name. A nested mock callback suffixes its outputs with `_out` (`energy_out`, `fock_out`, `error_out`); it may be pointer-assigned to a procedure pointer of `oao_object`/`arh_object` as long as the object is deallocated before the host returns.
- Helpers used by several procedures sit together near the top of the file, before the first `test_*` function, grouped by role so that each group only uses the ones above it: mocks and their support, then generators of random or structured data, then fixtures that set up objects or collect their state, then `ref_*` reimplementations, then `check_*` helpers. Within a group, helpers concerning the same object stay next to each other, in the order of the production code they serve (variants MO before OAO, closed before open shell, deterministic constructors ahead of random generators); `check_*` helpers follow the order of the tests they drive.
- Names follow the role: `generate_random_<what>` for a generator of random data, `<what>_matrix` for a deterministic constructor (`identity_matrix`, `diagonal_matrix`), `setup_<what>` for a fixture that sets up a module-global object.

**Verify a new test actually fails when the routine is broken.** Mutate the routine (flip a sign, swap an index, drop a term), rebuild, confirm the test fails, then restore. Mutating several routines at once and checking the failure set matches one-to-one is efficient for a whole suite. Some mutations are legitimately benign for a given input (scaling safeguards, guards on paths the test doesn't reach) — pick a different mutation rather than concluding the test is weak.

**Before finishing:** build with `-DCMAKE_BUILD_TYPE=Debug` (`-Wall -Wextra -fcheck=all`) and run the suite — it catches out-of-bounds accesses release silently tolerates. Because extension test sources are added per-`ENABLE_<EXT>` block, a new shared test module can break configurations you weren't working on — after touching `extensions/common/`, build with each extension enabled alone, plus none at all.

### CMake options that matter

- `INTEGER_SIZE` (`4` or `8`): library integer width. Unset → CMake autodetects, trying 32-bit BLAS/LAPACK first. Output library named `libopentrustregion_32.*` / `_64.*`. `USE_ILP64` (auto-set when `INTEGER_SIZE=8`) switches Fortran `ip` and C `c_ip` to 64-bit and remaps BLAS/LAPACK symbols to `_64` variants when `check_fortran_function_exists` finds them.
- `BLAS_LIBRARIES` / `LAPACK_LIBRARIES`: must be set together with `INTEGER_SIZE` if overriding autodetection — one without the other is a fatal error.
- `OpenTrustRegion_BUILD_TESTING` (default `ON` when top-level): builds `libotrtestsuite`, which Python loads to drive the Fortran tests.
- `CONDA_BUILD=1` env var: `setup.py` skips the embedded CMake invocation (the conda recipe builds the C library separately).

### Preprocessing and integer kinds

All Fortran sources (core and extensions) compile with `Fortran_PREPROCESS ON`. Integer kind selection and the BLAS/LAPACK 64-bit symbol remap (`ddot=ddot_64`, etc.) happen via `#ifdef USE_ILP64` and `add_compile_definitions` in CMake — never hardcode integer kinds.

BLAS/LAPACK are called through implicit interfaces (`external :: dgemm`), so gfortran infers each dummy argument's kind from the first call site and rejects a later call that disagrees. Both rules below are invisible in the default 32-bit build, where `ip` is `int32` and `1_ip == 1`:

- **Every integer argument to a BLAS/LAPACK routine must be `ip`-kinded** — `1_ip` not `1` for increments/leading dimensions, `size(x, kind=ip)` when passing on a `size()` result. Same for integers passed to project routines with an explicit `integer(ip)` dummy. Locals used as LAPACK arguments must be declared `integer(ip)`.
- **A newly used BLAS/LAPACK symbol must be added to the remap lists in `CMakeLists.txt`** (`add_compile_definitions("name=name_64")`, guarded by `BLAS_64`/`LAPACK_64`) — including test sources (`zheev` reaches the build only via `tests/opentrustregion_system_tests.f90`). A missing entry links the 32-bit-integer symbol from an ILP64 build, corrupting arguments at runtime rather than failing to build.

**Verifying an ILP64 build without an ILP64 BLAS.** Most dev machines only have 32-bit-integer BLAS, so `check_fortran_function_exists("sgemm_64")` fails and the remap is never exercised. Force it by pre-seeding the cache variables and checking the symbols the objects actually reference:

```sh
cc -shared -o /tmp/blas64stub.dylib /tmp/stubs.c   # one empty function per name_64_ symbol
cmake -S . -B /tmp/build_ilp64 -DINTEGER_SIZE=8 \
      -DBLAS_LIBRARIES=/tmp/blas64stub.dylib -DLAPACK_LIBRARIES=/tmp/blas64stub.dylib \
      -DBLAS_64=1 -DLAPACK_64=1 -DENABLE_OAO=ON -DENABLE_ARH=ON -DBUILD_SHARED_LIBS=ON
cmake --build /tmp/build_ilp64
find /tmp/build_ilp64 -name '*.o' | xargs nm -u | grep -E '_(d|z)[a-z]+_'   # none may lack _64
```

A clean link proves nothing is left unmapped; it does **not** verify numerical behaviour, which needs a real ILP64 BLAS.

### Known gotchas

- **A `size()`/literal-kind mistake in a BLAS call passes the default build and only breaks under `INTEGER_SIZE=8`**, as `Error: Type mismatch between actual argument at (1) and actual argument at (2) (INTEGER(4)/INTEGER(8))` pointing at two unrelated call sites of the same routine — the two that disagree, not the one that's wrong. Fix by making every integer argument `ip`-kinded, not by changing the named site.
- **`gfortran -fsyntax-only -I<moddir>` gives false confidence.** It checks only the pointed-to file against whatever `.mod` files already exist — it doesn't re-verify those `.mod`s. Editing a `type` whose fields are used across modules can leave a stale consumer `.mod` "passing" syntax-only checks while a real build breaks with `Fatal Error: Mismatch in components of derived type '...': expecting 'X', but got 'Y'`. Always confirm interface changes with `cmake --build`, not `-fsyntax-only`.
- **A `build/` directory created by `pip install` can't be rebuilt directly later.** pip's ephemeral `cmake` path gets baked into `build/CMakeCache.txt` (`CMAKE_COMMAND`) and generated Makefile stamp rules. Once pip's temp env is gone, `cmake --build build` fails with `<temp-path>/cmake: No such file or directory`. Diagnose with `grep CMAKE_COMMAND build/CMakeCache.txt`. Fix: re-run `pip install -e .`, or maintain a separate manually-configured build directory with the system `cmake`.
- **The Python driver only looks for `../build`.** `python_interface.py` searches site-packages, then `pyopentrustregion/`, then `<repo>/../build` — a manually-configured `build_manual/` is invisible to it, and it silently loads whatever's in `build/` instead, possibly with different extensions enabled. To drive a custom build directory, run from a scratch directory with symlinks named `pyopentrustregion` and `build`:

  ```sh
  mkdir -p /tmp/run && cd /tmp/run
  ln -sfn <repo>/pyopentrustregion pyopentrustregion
  ln -sfn <repo>/build_manual build
  python3 -m pyopentrustregion.testsuite
  ```
- **The core test name list in `pyopentrustregion/tests.py` has no import-time guard.** `testsuite` wraps each extension's test import in `try/except AttributeError`, so a disabled extension is silently skipped — but the core list does `getattr(lib, f"test_{name}")` unguarded at module import time. A drifted list (a Fortran test renamed/removed without updating it) raises `AttributeError: dlsym(...): symbol not found` and aborts the *entire* run, including unrelated extension tests. Keep the name list and the Fortran `test_*` functions in exact sync.
- **`ref_settings_type_c` in `test_reference.f90` has an implicit, unenforced field-order convention.** `PyInterfaceTests.setUpClass` in `pyopentrustregion/tests.py` rebuilds a matching struct by walking `SolverSettings.c_struct._fields_ + StabilitySettings.c_struct._fields_` (that concatenation order, bools-then-reals-then-ints, then character arrays in relative order) and passes it into Fortran's `get_reference_values` via a raw `POINTER(RefSettingsC)` cast. So `ref_settings_type_c` must declare all `solver_settings_type_c`-owned character/string fields before all `stability_settings_type_c`-owned ones, regardless of what feels natural. Getting it wrong causes no crash — the two structs read across each other's byte ranges, so a field silently reports another field's value, or raises `UnicodeDecodeError` if the misaligned bytes aren't valid UTF-8. A new settings field failing with an inexplicable value mismatch or decode error → suspect field order here first.

## Core library

### Source layout

- `src/opentrustregion.f90` — the numerical core: `solver`, `stability_check`, the Davidson/Jacobi-Davidson/truncated-CG subsystem solvers, settings derived types, callback abstract interfaces, error-code constants. Single module, several thousand lines.
- `src/c_interface.f90` — `bind(C)` wrapper module. Stores Fortran procedure pointers to C callbacks at module scope (`update_orbs_before_wrapping`, etc.) and adapts C-style `(*)` arrays + return-code functions into Fortran-style `(:)` arrays + `intent(out) :: error` subroutines.
- `include/opentrustregion.h` — C header mirroring `solver_settings_type` / `stability_settings_type` as C structs, plus `*_init()` helpers and `solver`/`stability_check` prototypes.
- `pyopentrustregion/python_interface.py` — ctypes wrapper. `SolverSettings`/`StabilitySettings` as `ctypes.Structure` mirrors of the C structs, Python callbacks wrapped with `CFUNCTYPE`, integer error codes converted to `RuntimeError`.
- `tests/` — `opentrustregion_unit_tests.f90` / `c_interface_unit_tests.f90` (Fortran/C interface unit tests, both in `libotrtestsuite`), `opentrustregion_system_tests.f90` (system tests against `pyopentrustregion/test_data/`), `c_system_tests.c` (system tests through the C interface), `*_mock.f90` (mock callbacks for both unit suites), `test_reference.f90` (tolerances + reference values).

### Fortran/C/Python interfaces must stay consistent

The same callback signatures, settings fields, and defaults are described in seven places that must agree — any change to a callback signature, a settings field (add/remove/reorder), or a default value must land in all seven in the same PR:

1. Fortran abstract interfaces and `solver_settings_type`/`stability_settings_type` in `src/opentrustregion.f90`.
2. C abstract interfaces and `bind(C)` `solver_settings_type_c`/`stability_settings_type_c` in `src/c_interface.f90`.
3. C struct layouts/typedefs in `include/opentrustregion.h`.
4. `SolverSettingsC`/`StabilitySettingsC` ctypes `_fields_` and `CFUNCTYPE` declarations in `pyopentrustregion/python_interface.py`.
5. `default_solver_settings`/`default_stability_settings` (Fortran) ↔ `solver_settings_init`/`stability_settings_init` (C) ↔ the Python `Settings` wrapper defaults.
6. Argument lists and snippets in `README.md`.
7. `ref_settings_type`/`ref_settings_type_c` in `tests/test_reference.f90`, used by the settings round-trip tests in `pyopentrustregion/tests.py` — see the field-*order* gotcha above.

Error-origin codes (`error_obj_func` etc.) must also stay synchronized with the README error-code table.

### Error codes

Encoded as `OOEE` integers (origin × 100 + specific code, see README). The Fortran core uses `add_error_origin` to tag a non-zero callback error with its origin (e.g. `error_update_orbs = 1200`). Extensions reuse these core callback-origin codes rather than defining their own. C entry points return the integer directly; Python wrappers raise a `RuntimeError` that includes the raw code as-is — they don't decode it into origin/specific parts, so callers must read the README table themselves. New origins: add as a parameter in `opentrustregion.f90` and keep the README table in sync.

## Extensions

Opt-in, default OFF (`ENABLE_OAO`, `ENABLE_MO`, `ENABLE_ARH`, `ENABLE_QUASI_NEWTON`, `ENABLE_S_GEK`). `ENABLE_ARH` forces `ENABLE_OAO ON` and `ENABLE_MO ON` in `CMakeLists.txt`. Enable via `CMAKE_FLAGS='-DENABLE_ARH=ON -DBUILD_SHARED_LIBS=ON' pip install -e .` — without `BUILD_SHARED_LIBS=ON` there's no `.dylib`/`.so` for ctypes to load, and the extension's Python wrapper raises `RuntimeError: Please reinstall the package with: CMAKE_FLAGS=...`.

### Source layout

- `extensions/<name>/tests/` — mirrors core's `tests/` per extension: `<name>_unit_tests.f90`, `<name>_c_interface_unit_tests.f90`, `<name>_mock.f90` (mocks matching that extension's abstract interfaces), `<name>_test_reference.f90` (reference settings, funptr checkers, `ref_*` reimplementations of the extension's own routines).
- `extensions/common/src/common_c_interface.f90` (module `otr_common_c_interface`) — shared *production* C-interface code: the `update_orbs`/`hess_x` C-wrapper implementations reused as-is by every extension's `<name>_c_interface.f90`, since those types and their C bridge are defined once in the core, and the `*_impl` routines of the exact-evaluation callbacks (`evaluate_dm`, `get_response`) and of the objective function and solver-setting callbacks of the orbital bases, which take the extension's module-level callback slot as an argument. Their C typedefs are in `extensions/common/include/opentrustregion_common.h`, which `opentrustregion_oao.h` and `opentrustregion_mo.h` include. On the real call path, not test support.
- `extensions/common/src/common.f90` (module `otr_common`) — shared *production* numerics of the OAO and the MO basis, above all the abstract `orbital_basis_type` both orbital objects extend and the `orbital_settings_type` they hold. Listed only in the `ENABLE_OAO`, `ENABLE_MO` and `ENABLE_ARH` blocks, together with its tests in `extensions/common/tests/common_unit_tests.f90`, registered in `pyopentrustregion/extensions/common/tests.py` (skipped via `AttributeError` when none is built).
- `extensions/common/tests/` — `common_mock.f90`, the callback mocks of the C-interface tests shared by every extension or by the orbital bases, listed in every `if(ENABLE_<EXT>)` CMake block (CMake dedupes the repeated source entry); `common_test_reference.f90` and `common_unit_tests.f90`, the tests of `otr_common` plus the test code OAO's, MO's and ARH's tests share (see Test code layout above), listed only in the `ENABLE_OAO`, `ENABLE_MO` and `ENABLE_ARH` blocks; `common_c_interface_unit_tests.f90`, the `bind(C)` mocks of the exact-evaluation callbacks the C-interface tests of the orbital bases share, listed in the same blocks, since their C types exist only with an orbital basis.
- `pyopentrustregion/extensions/common/python_interface.py` — the Python counterparts of the shared exact-evaluation interfaces and the factory-returned callback classes both orbital bases and ARH build on (`EvaluateDMInterface`, `ObjFuncPyInterface`, `UpdateOrbsPyInterface`, `attach_wired_callbacks`).

### Extensions repeat the core's six-location interface pattern

Each extension (`arh`, `mo`, `oao`, `quasi_newton`, `s_gek`) has its own interface chain — a callback signature change touches all six:

1. Abstract interface + settings type in `extensions/<name>/src/<name>.f90` (the exact-evaluation interfaces every orbital basis shares in `extensions/common/src/common.f90`).
2. `bind(C)` interface + wrapper in `extensions/<name>/src/<name>_c_interface.f90` (the shared ones in `extensions/common/src/common_c_interface.f90`).
3. C typedef/struct in `extensions/<name>/include/opentrustregion_<name>.h` (the exact-evaluation typedefs every orbital basis shares in `extensions/common/include/opentrustregion_common.h`).
4. `CFUNCTYPE`/`Structure` in `pyopentrustregion/extensions/<name>/python_interface.py` (the shared ones in `pyopentrustregion/extensions/common/python_interface.py`).
5. Mock callbacks in `extensions/<name>/tests/<name>_mock.f90` / `<name>_c_interface_mock.f90`.
6. Reference values in `extensions/<name>/tests/<name>_test_reference.f90` and the Python mock in `pyopentrustregion/extensions/<name>/tests.py`.

If an extension doesn't actually need a given callback, remove it from all six rather than leaving it defined-but-unused — a present-but-ignored parameter looks load-bearing to callers, which will build and pass a closure for nothing.

When a factory-style C interface must accept either of two distinct callback signatures for the same slot (e.g. two shells whose callbacks take different arguments), expose it as a named union of the two typedefs with members named for the two cases, not a generic `void*` or a single opaque typedef — keeps the header type-safe while giving the factory one parameter slot. Where the two signatures coincide, as the closed- and open-shell `evaluate_dm` of the ARH extension do, a single typedef serves both.

### OAO, MO and ARH

The chemistry-specific extensions: OAO builds orthogonal-atomic-orbital machinery on density matrices; MO holds the orbitals parameterized by the occupied-virtual rotations of MO coefficients, with exact-Hessian callbacks like OAO's; ARH (augmented Roothaan-Hall) is a history-based approximate-Hessian method in either the OAO or the MO basis. Quasi-Newton and S-GEK have no equivalent chemistry-specific state — the conventions below don't apply to them.

**Naming.**
- `_cs`/`_os` for closed-shell/open-shell variants, consistently in Fortran, C, and Python (not `_closed_shell`/`_open_shell` or rank-based `_2d`/`_3d`).
- When two routines, types, variables or tests do the same thing in the two orbital bases, both carry the basis after the base name and before the spin and role suffixes: `<name>_<basis>_<spin>_<role>` (`arh_factory_mo_cs`, `arh_factory_mo_c_wrapper`). Implementations of the `orbital_basis_type` methods are `<method>_oao`/`<method>_mo`, those of the ARH methods `<method>_arh_oao`/`<method>_arh_mo`. The `oao_*` prefix of `otr_oao`'s own routines (`oao_factory`, `oao_sanity_check`) is a namespace, not a basis qualifier.
- Procedures wired into the solver settings or returned by the factories are not type-bound and reach their object through the module-level global, which the suffix `_callback` marks: `<setting>_<extension>_callback` (`update_orbs_oao_callback`, `obj_func_arh_callback`); their interface-checking procedure pointers are `<callback>_ptr`, their tests `test_<callback>`. A callback dispatching to a method shares the method's name up to the suffix.

**Ordering.** Routines, tests, mocks, declarations and C/Python entry points: variants of one routine stay next to each other, MO before OAO (MO is the default basis), closed shell before open shell, and a variant with a role suffix after those without (`arh_factory_mo_cs`, `arh_factory_mo_os`, `arh_factory_mo_common`, `arh_factory_oao_cs`, ...); tests follow the order of the routines they test.

**Factories read the same in all languages.** `arh_factory_mo`, `arh_factory_oao`, `mo_factory` and `oao_factory` exist under these names in Fortran, C and Python, with the same arguments in the same order, the sizes included even where Python could infer them from the array shapes — as generics over the two shells in Fortran and one function each in C and Python. On top of them, Fortran's generic `arh_factory` covers the four specific ARH factories (`arh_factory_mo_cs`/`_os`, `arh_factory_oao_cs`/`_os`) and Python's `arh_factory` helper forwards to whichever basis factory's signature its arguments bind to; C cannot overload and has no `arh_factory`.

**Python checks only what the C interface cannot.** The C interface passes bare pointers, so the Python factories check the array shapes against the sizes, the length of the occupation list (ctypes pads a short one with zeros), the memory layout (a density matrix the factory updates in place must be a writeable, contiguous `float64` array; MO coefficients must be a writeable real array) and that `evaluate_dm` can be called with the positional arguments the factory's callback interface passes, which `check_callback_arguments` (`pyopentrustregion/extensions/common/python_interface.py`) checks through `inspect.signature(...).bind`, so that defaulted and variadic arguments are accepted. MO coefficients are passed through a column-major `float64` buffer: the caller's array itself if it is stored like that, which every factory rotates through its own view of it, and otherwise a copy copied back after every orbital update, which `check_mo_coeff` (`pyopentrustregion/extensions/mo/python_interface.py`) keeps per caller array object, referenced only weakly so that it never keeps the array alive, and hands to every MO-basis factory called with it, refreshed from the array: `mo_object` points to the buffer of the latest factory call, so a second buffer would let the copy-back of an earlier factory's orbital updating function overwrite the caller's orbitals with stale ones. Both factories therefore have to be passed the same array object, not two views of it. Value ranges (`n_particle` being 1 or 2, a positive `n_ao`, occupations fitting into the MOs) are checked only by the Fortran sanity checks and tested only in their tests, never again in Python.

**Testing routines with module-global state.** `oao_object`/`mo_object`/`arh_object` are module-level allocatables whose components are mostly pointers. Populate them directly from test-local `target` arrays rather than calling the factory (which would violate the no-untested-dependencies rule and drag in the whole setup chain). Setting `s_inv_sqrt` to the identity makes the AO and OAO bases coincide, removing the basis transformation from expected values. Call the `orbital_basis_type` methods on the object directly, and build ARH objects from the structure constructors (`arh_mo_type(mo_object)`, `arh_oao_type(oao_object)`) only after the object's quantities are allocated, since the constructors associate only allocated ones. Always `deallocate` the global object at the end of the test — otherwise its pointers dangle into freed test locals and later tests see leftover state. gfortran can't see the later deallocation and emits `Warning: Pointer at (1) in pointer assignment might outlive the pointer target [-Wtarget-lifetime]` in a Debug build for each such assignment; expected, not a defect, as long as the matching `deallocate` is present on every exit path.

**Orbital bases and ARH.** The quantities and operations every orbital basis has are components and type-bound procedures of `otr_common`'s abstract `orbital_basis_type`: the settings, the staleness flags and the exact-evaluation callbacks it holds, deferred procedures every basis implements, and shared procedures built on those (`refresh_response`, `precond`, `precond_pd`, `fill_extra_trial_vectors`), which the callbacks of every basis and ARH call instead of repeating them. It is extended by `oao_type` (`otr_oao`) and `mo_type` (`otr_mo`). ARH is the abstract `arh_type`, holding the settings, the history, the approximate Hessian model and `orbitals`, a pointer at `oao_object` or `mo_object`; `arh_mo_type`/`arh_oao_type` extend it with the ARH-only operations. `arh_object` is `class(arh_type), allocatable`, so shared code in `otr_arh` dispatches automatically and never branches on the basis. The ARH types cannot extend `oao_type`/`mo_type` instead, since the history state would have no common ancestor to dispatch on and in the OAO basis ARH shares `oao_object` with OAO's exact callbacks.
- Methods read the dimensions from `self` and write their own components; a component a method writes is never also passed as an argument, which would alias it with `self`.
- `dm_ao` is a pointer: OAO points it at the caller's density matrix, which it updates in place, while `mo_type` allocates it through the pointer and frees it in its finalizer, so an `mo_type`'s `dm_ao` must never be pointed at another array.
- In the MO basis, the factories never touch the projection of the solver settings, so that a caller's projection (e.g. a symmetry projection) is kept, consistently in Fortran, C and Python; solver settings still carrying OAO's projection are a usage error, which is not detected. The OAO basis still replaces a caller's projection with `project_oao_callback`.
- Symmetry in the MO basis is the optional last argument `orbsym` of `mo_factory` and `arh_factory_mo`, the irreps of the MOs of every channel, in Fortran, C (`NULL` if absent) and Python (`None`). Every `mo_channel_type` stores its `irreps` and the `param_mask` of its same-irrep occupied-virtual pairs (all true without symmetry), and every packing of parameters goes through `pack_ov`/`unpack_ov`/`mo_param_rows`, never through `reshape` or `n_occ * n_virt` counts, so that symmetry needs no other code path; `count_mo_params` counts the parameters also where no MO object holds them, in the sanity check and the C wrappers. `diagonalize_per_irrep` diagonalizes the Fock blocks per irrep and stores every eigenvector at the orbitals of its irrep, so that rotations into the eigenbasis keep the mask pattern and degenerate eigenvalues of different irreps cannot mix them. The MO tests cover symmetry through the shared occupation cases (`case_irreps`); their `ref_*` helpers derive the same-irrep pairs from the channels' irreps, never from `param_mask`.
- The MO history is stored in the AO basis (the MO basis changes every update), densities as `S D S` so that they transform to the MO basis like the potentials (`C^T M C`). Its history columns must stay the *full* `n_mo x n_mo` difference matrices, not their occupied-virtual blocks: differences a finite step apart have second-order oo/vv blocks the potentials respond to, and an ov-only history converged several times slower. Only the packed columns (the parameter space) are ov blocks.

**Re-entry: the orbital objects are reused, the ARH object never is.** `oao_factory_common` and `mo_factory_common` keep an existing object's allocations if it was already set up for the same dimensions (and, for MO, occupations and irreps), otherwise set it up anew; either way they take over the starting orbitals and reset the evaluation state. They detect a completed setup from a quantity allocated only once setup succeeded, and read the stored dimensions only after that, since the dimensions have no default. The reuse is not an optimization: `oao_object` is shared between `oao_factory` and `arh_factory_oao`, which PySCF calls one after the other, so the second setup must not discard what the first stored (above all the `evaluate_dm` callbacks the exact callbacks need, which only `oao_factory_cs`/`oao_factory_os` and `mo_factory_cs`/`mo_factory_os` set, each clearing the other shell's, since the exact callbacks dispatch on which one is associated); `mo_object`, shared in the same way between `mo_factory` and `arh_factory_mo`, follows the same rules. Since an orbital basis and ARH in that basis also share the C callback slots of the objective function, the preconditioners, the extra trial vectors and, in the OAO basis, the projection, the factory called last owns these callbacks for both, so a caller combining them calls `oao_factory`/`mo_factory` before `arh_factory_oao`/`arh_factory_mo`. The ARH factories instead always deallocate and reallocate `arh_object`, because leftover history from an unrelated trajectory must never survive, even if the dimensions match. All Fortran factories run their sanity checks before changing any state. Don't simplify one re-entry rule to match the other.

**Evaluation staleness is an explicit flag, never inferred from values.** The orbital updating callbacks evaluate without an orbital rotation only when `evaluation_stale` (in the orbital objects also `response_stale`) is set: it starts `.true.`, is set before the density is rotated and cleared only once the rotated density has been evaluated, so a failed evaluation leaves it set, and ARH adds the current point to its history only while it is clear. The orbital objects' `response_stale` likewise starts `.true.` and is set again by the factory and whenever `rotate_orbitals` moves the orbitals without rebuilding the response (as ARH does), which makes the exact Hessian refresh the response first. Don't reintroduce proxies such as a vanishing energy or unallocated arrays, which can be legitimate states, and don't give other components default initializers to serve as such markers.
