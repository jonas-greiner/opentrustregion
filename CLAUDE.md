# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

OpenTrustRegion is a Fortran library implementing a second-order trust-region optimizer, exposed via Fortran, C, and Python (ctypes) interfaces — the same compiled shared/static library is consumed by all three. The core solver (`src/`) is agnostic to what is being optimized: it operates purely on abstract parameter/gradient/Hessian-vector-product callbacks. Domain-specific functionality is layered on top as opt-in extensions under `extensions/<name>/`. Two of these, OAO and ARH, are chemistry-specific — they build density-matrix and orbital-optimization machinery on top of the core solver. The other two, quasi-Newton and S-GEK are, like the core, agnostic to the underlying problem; they just add alternative Hessian-approximation strategies.

## General conventions

Apply equally to `src/` and every `extensions/<name>/`.

**Formatting.** Fortran lines max 88 columns, continued with a trailing `&` aligned to the opening parenthesis of the call. Every procedure opens with a `!`-delimited comment block describing what it does; each logical step inside gets a lowercase `!` comment. Unit tests: assume success (`test_<name> = .true.`), then one `if (...) then` / `write (stderr, *) "test_<name> failed: ..."` / `test_<name> = .false.` block per assertion, no early returns except where a later assertion would crash.

**Python formatting.** All Python (`pyopentrustregion/`, its `extensions/`, `setup.py`) is `black`-formatted, default settings. Run `black pyopentrustregion setup.py` before considering Python changes done.

**C formatting.** All C headers and test sources (`include/`, `tests/*.c`) are `clang-format`-formatted per the repo-root `.clang-format` (LLVM style, 88-column limit to match the Fortran convention above). Run `clang-format -i` on touched C/H files before considering C changes done.

**Clarity over performance.** The library is not on the hot path of a quantum-chemistry calculation — the host program's integral transforms and Hessian-vector products dominate. Prefer short, obviously-correct code over fast code. Don't propose performance refactors (buffer growth, pooling, micro-optimizations) without evidence the affected code is hot for a real workload.

**Norms: BLAS in production, intrinsics in tests.** Production Fortran (`src/`, `extensions/*/src/`) uses BLAS (`dnrm2`, `ddot`) for norms/dot products. Test code (`tests/`, `extensions/*/tests/`) uses intrinsics (`norm2`, `dot_product`, `matmul`, `sum`, `transpose`) instead. BLAS/LAPACK is allowed in tests only where no intrinsic exists (`dsyev`, `dgeev`, `zheev`) — don't hand-roll linear algebra to avoid the dependency.

### Unit tests must not depend on untested code

A unit test may only call the routine it is testing — never another production routine to build its inputs or expected values, or the test silently inherits that routine's correctness.

- Construct inputs directly in the test, or with LAPACK, not by calling the production routine that would normally produce them (populate `arh_object`'s cached history projections with `ref_cache_dirs`, not `cache_history_projections`; inject a random symmetric metric pseudoinverse and `A_sym` core directly instead of deriving them with `get_arh_metric_inv` / `build_a_sym_cs`).
- When a routine under test internally calls another module's routines, reimplement those as `ref_*` helpers so the expected value is computed independently — e.g. `ref_unpack_asymm` / `ref_pack_asymm` / `ref_project_asymm` / `ref_project_symm` / `ref_hess_x` in `extensions/oao/tests/oao_unit_tests.f90`, used by both OAO's and ARH's Hessian tests.
- Where a routine merely delegates (e.g. `update_orbs_*` returning the gradient `calculate_grad_h_diag` produced), assert the *contract* — the returned arrays equal what the routine stored — rather than re-deriving numbers a different routine's own test already covers.

### Test code layout

Classify test code by **role** first, then place it at the most common **level** that needs it.

Roles, one file each:

| File | Role |
|---|---|
| `test_reference.f90` | Tolerances, reference values, `ref_*` reimplementations computed independently of the routine under test. |
| `<name>_unit_tests.f90` | The Fortran-level `test_*` functions, plus mocks/fixtures needed only to drive them (e.g. a stand-in `get_energy`/`update_dm`). |
| `<name>_mock.f90` | Mocks of the module's *own* production routines that are bridged to the C interface (factories, deconstructors), so a caller (typically `<name>_c_interface_unit_tests.f90`) can verify invocation without running the real logic. |
| `<name>_c_interface_unit_tests.f90` | Tests for the `bind(C)` wrapper layer; defines its own local `bind(C)` mocks rather than using `<name>_c_interface_mock.f90`. |
| `<name>_c_interface_mock.f90` | `bind(C)`-signature mocks used exclusively by the Python interface tests, dynamically loaded via `ctypes` from `libotrtestsuite`. |

Level: within `tests/`, `extensions/common/tests/`, and a specific extension's `tests/`, code lives at the most common level shared by everything that *actually* needs it — verified against real usage across every extension, not assumed from whichever extension motivated writing it. Needed by every extension → `extensions/common/tests/` (e.g. `common_mock.f90`'s `mock_update_orbs`/`mock_hess_x`, since `update_orbs_type`/`hess_x_type` and their C bridge are defined once in the core). Needed by some → the most foundational one among them (ARH depends on OAO, so a helper shared only between OAO's and ARH's tests lives in `oao_unit_tests.f90`, not `extensions/common/`). Needed by one → that extension's own files. The mock logger and `setup_settings`, needed by nearly every Fortran-level unit test, live in `tests/opentrustregion_unit_tests.f90`.

This recurses within a role: a helper belongs in `<name>_test_reference.f90` vs `<name>_unit_tests.f90` based on whether something *at that level* still calls it — not which extension's tests happen to use it. `ref_unpack_asymm` etc. live in `oao_unit_tests.f90` rather than `oao_test_reference.f90` because nothing in `test_reference.f90` calls them anymore; only `oao_unit_tests.f90` and `arh_unit_tests.f90` (importing from `otr_oao_unit_tests`) do.

**`ref_*` naming.** The prefix marks a role — an independent reimplementation used to check a routine's expected value — not a location; it travels with the function wherever the placement rule moves it. What decides the prefix is what the helper *computes*, not how a call site uses it: `ref_cache_dirs` independently reproduces what `cache_history_projections` would return, so it keeps the prefix even where a test feeds its result straight into `arh_object` as an input fixture rather than comparing against a production call. A pure fixture that builds a *plausible* input (via LAPACK, or by reassembling a routine's own output) does not take the prefix, even in the same file — e.g. ARH's `generate_fock_partition`, `generate_nonredundant_vector`, `generate_random_dm_diff`, `generate_random_symm_hessian`. Either kind may keep an atypical name purely to dodge a collision with an identically-named production routine or local variable in scope at its call site (`ref_unpack_asymm` etc. keep the prefix in `oao_unit_tests.f90` because bare `unpack_asymm` etc., and in `arh_unit_tests.f90` local `hess_x` variables, are already in scope) — check for a collision before insisting on the "cleaner" name.

**Module-level `use` in test files is restricted.** In `<name>_unit_tests.f90` and siblings, module-level `use` statements above `contains` are limited to `rp`, `ip`, `kw_len`, `stderr`, `stdout`, `tol`, `tol_c`, and intrinsic module bindings (`iso_c_binding` and the like) — check other modules for the exact set before adding to it. Everything else (mock callbacks, `setup_settings`, shared dimension parameters, `ref_*` helpers) is imported locally inside the specific `test_*` function that needs it, even if several functions in the file need the same symbol. A module-level `use` beyond this set pollutes every procedure in the file and hides which test actually relies on what.

**Prefer shared dimension-parameter fixtures over local literals**, unless it complicates the test a lot. Extension `<name>_test_reference.f90` modules define canonical dimensions (e.g. OAO's `n_ao = 3`, `n_particle = 2`, derived `n_param` in `otr_oao_test_reference`) that most unit tests in that extension should import rather than redeclaring. `n_param` specifically is only right for a test that works at the shared `n_particle` — a test pinning its own `n_particle` (closed-shell-only, or one walking a mutable local `n_particle`/`n_param` through both shells in one body, e.g. `test_calculate_grad_h_diag`, `test_hess_x_oao`, `test_hess_x_arh`) needs its own local value. Exception: a test whose expected values are hand-derived for a specific small/structured case (a closed-form 2×2 rotation, a matrix crafted to trigger a rank-detection path) may keep local literals, since forcing the shared size would mean re-deriving the numbers from scratch. When correctness is checked via algebraic properties rather than exact values, there's no such barrier — use the shared fixture.

**Prefer generated data over hand-typed literals** (`generate_random_density_matrix` / `generate_random_symm_matrix` / `call random_number`), unless the specific values are load-bearing. A hand-typed `reshape([...])` is justified only when the test checks a specific hand-derivable closed-form target depending on those exact numbers (a weighted-symmetrization blend, a dependency-screening path, a regularization threshold crossed at a specific eigenvalue). When every assertion is relational (a copy, a fixed scaling, an algebraic identity, or a value the test recomputes independently from the same input) the specific values never mattered, and random generation exercises the same code path while being harder to accidentally pass with a bug. Before hand-typing a matrix meant to be "just some valid state" (often signaled by the test calling it "arbitrary"), check whether random generation works instead. Where a matrix must stay block-diagonal/diagonal/structurally constrained to match what production produces, keep that structure but randomize within it.

### Shared architectural patterns

**Settings types.** `solver_settings_type` / `stability_settings_type` (core) and each extension's own settings type (`oao_settings_type`, `s_gek_settings_type`, `qn_settings_type`, `arh_settings_type` — which extends `oao_settings_type` since ARH depends on OAO) all extend the abstract `settings_type`, with an `init` type-bound procedure and an `initialized` flag. The owning routine checks `settings%initialized` on entry and calls `init` if needed; C/Python wrappers pre-populate defaults via `*_init()`/`__init__` so the flag is `.true.` by the time Fortran sees it. Default values are defined once in a `default_*_settings` parameter and replicated in the C `*_init()` function and the Python `Structure` defaults — keep them in sync.

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

**Exactly one registered test per production routine** — never `test_<routine>_<case_a>`, `test_<routine>_<case_b>`. Cover multiple configurations (each `arh_type`, closed vs. open shell, each preconditioner) as cases *inside* that one test: loop where setup can be shared, or call a per-case `check_*` helper that is deliberately not `bind(C)` and not registered. Keep the case name in the failure message (`"test_<routine> failed for <case>: ..."`) since the harness only reports the registered name. Worked example: `test_inv_hess_x_arh` — ten cases, one registered entry, driven by `check_inv_hess_x_arh_cs` / `check_inv_hess_x_arh_os`, each parameterized over `arh_type` and looped over all five types.

A `check_*` (or other non-`bind(C)`) helper earns its existence by eliminating real duplication, never merely to keep one case's code out of the registered test's body — five near-identical one-case functions just relocate the duplication; write one function parameterized over the case, or a loop with a shared body. Reserve separate helper routines for pieces reused *across* cases (a shared fixture, a shared assertion). Place any such helper — fixture or `check_*`/`ref_*` — together near the top of the file, right after the mocks and before the first `test_*` function, not interleaved next to the one registered test that calls it.

Because a routine's cases now live behind one entry, a silently-skipped case is invisible in suite output — have the driving test run every case unconditionally and `and` the results together, rather than returning early on the first failure.

**Verify a new test actually fails when the routine is broken.** Mutate the routine (flip a sign, swap an index, drop a term), rebuild, confirm the test fails, then restore. Mutating several routines at once and checking the failure set matches one-to-one is efficient for a whole suite. Some mutations are legitimately benign for a given input (scaling safeguards, guards on paths the test doesn't reach) — pick a different mutation rather than concluding the test is weak.

**Before finishing:** build with `-DCMAKE_BUILD_TYPE=Debug` (`-Wall -Wextra -fcheck=all`) and run the suite — it catches out-of-bounds accesses release silently tolerates. Because extension test sources are added per-`ENABLE_<EXT>` block, a new shared test module can break configurations you weren't working on — after touching `extensions/common/tests/`, build with each extension enabled alone, plus none at all.

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

Opt-in, default OFF (`ENABLE_OAO`, `ENABLE_ARH`, `ENABLE_QUASI_NEWTON`, `ENABLE_S_GEK`). `ENABLE_ARH` forces `ENABLE_OAO ON` in `CMakeLists.txt`. Enable via `CMAKE_FLAGS='-DENABLE_ARH=ON -DBUILD_SHARED_LIBS=ON' pip install -e .` — without `BUILD_SHARED_LIBS=ON` there's no `.dylib`/`.so` for ctypes to load, and the extension's Python wrapper raises `RuntimeError: Please reinstall the package with: CMAKE_FLAGS=...`.

### Source layout

- `extensions/<name>/tests/` — mirrors core's `tests/` per extension: `<name>_unit_tests.f90`, `<name>_c_interface_unit_tests.f90`, `<name>_mock.f90` (mocks matching that extension's abstract interfaces), `<name>_test_reference.f90` (reference settings, funptr checkers, `ref_*` reimplementations of the extension's own routines).
- `extensions/common/src/common_c_interface.f90` (module `otr_common_c_interface`) — shared *production* C-interface code: the `update_orbs`/`hess_x` C-wrapper implementations reused as-is by every extension's `<name>_c_interface.f90`, since those types and their C bridge are defined once in the core. On the real call path, not test support.
- `extensions/common/tests/` — test support genuinely shared by every extension: `common_mock.f90`'s `mock_update_orbs`/`mock_hess_x`. Listed in every `if(ENABLE_<EXT>)` CMake block (CMake dedupes the repeated source entry). Numerical fixtures needed only by OAO and ARH live in `oao_unit_tests.f90`/`arh_unit_tests.f90` instead (see Test code layout above).

### Extensions repeat the core's six-location interface pattern

Each extension (`arh`, `oao`, `quasi_newton`, `s_gek`) has its own interface chain — a callback signature change touches all six:

1. Abstract interface + settings type in `extensions/<name>/src/<name>.f90`.
2. `bind(C)` interface + wrapper in `extensions/<name>/src/<name>_c_interface.f90`.
3. C typedef/struct in `extensions/<name>/include/opentrustregion_<name>.h`.
4. `CFUNCTYPE`/`Structure` in `pyopentrustregion/extensions/<name>/python_interface.py`.
5. Mock callbacks in `extensions/<name>/tests/<name>_mock.f90` / `<name>_c_interface_mock.f90`.
6. Reference values in `extensions/<name>/tests/<name>_test_reference.f90` and the Python mock in `pyopentrustregion/extensions/<name>/tests.py`.

If an extension doesn't actually need a given callback, remove it from all six rather than leaving it defined-but-unused — a present-but-ignored parameter looks load-bearing to callers, which will build and pass a closure for nothing. Example: ARH's `get_response` was dropped entirely, since ARH approximates the Hessian response from its own history rather than calling a supplied function.

When a factory-style C interface must accept either of two distinct callback signatures for the same slot (e.g. closed-shell vs. open-shell `update_dm`), expose it as a named union of the two typedefs with members named for the two cases (`update_dm_fp` with `.cs`/`.os` members in `extensions/arh/include/opentrustregion_arh.h`), not a generic `void*` or a single opaque typedef — keeps the header type-safe while giving the factory one parameter slot.

### OAO and ARH

The chemistry-specific extensions: OAO builds orthogonal-atomic-orbital machinery on density matrices; ARH (augmented Roothaan-Hall) is a history-based approximate-Hessian method built on top of OAO. Quasi-Newton and S-GEK have no equivalent chemistry-specific state — the conventions below don't apply to them.

**`_cs`/`_os` suffix for closed-shell/open-shell variants**, consistently, in OAO/ARH Fortran, C, and Python. The only spelling in use — don't reintroduce `_closed_shell`/`_open_shell`, rank-based names like `_2d`/`_3d`, or an extension-local pair like ARH's old `_nonlinear`/`_spin`.

**Testing routines with module-global state.** `oao_object`/`arh_object` are module-level allocatables whose components are mostly pointers. Populate them directly from test-local `target` arrays rather than calling the factory (which would violate the no-untested-dependencies rule and drag in the whole OAO setup chain). Setting `s_inv_sqrt` to the identity makes the AO and OAO bases coincide, removing the basis transformation from expected values. Always `deallocate` the global object at the end of the test — otherwise its pointers dangle into freed test locals and later tests see leftover state. gfortran can't see the later deallocation and emits `Warning: Pointer at (1) in pointer assignment might outlive the pointer target [-Wtarget-lifetime]` in a Debug build for each such assignment; expected, not a defect, as long as the matching `deallocate` is present on every exit path.

**`arh_factory_common` vs `oao_factory_common` re-entry handling differs on purpose.** `arh_factory_common` unconditionally deallocates and reallocates `arh_object` (no `allocated()` check) because ARH's factory only ever starts a fresh calculation — leftover history from a previous, unrelated trajectory must never survive, even if dimensions happen to match. `oao_factory_common` keeps a more nuanced reuse check because it can legitimately be invoked twice within the same calculation (directly, and via ARH's own setup) and redoing that setup would be wasted. Don't simplify one to match the other without checking which reuse case applies.
