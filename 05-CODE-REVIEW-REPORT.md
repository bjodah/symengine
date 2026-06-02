# Code Review Report: F5B, MoGVW, and M4GB Groebner Changes

## Scope

Reviewed the unstaged branch changes against `04-IMPROVED-IMPL-PLAN-MOGVW-AND-M4GB.md`, with focus on dispatch behavior, GF(p) arithmetic, signature algorithm invariants, M4GB polymatrix invariants, tests, and benchmarks.

The review covered:

- `symengine/polys/groebner.{h,cpp}`
- `symengine/polys/groebner_internal.h`
- `symengine/polys/groebner_f5b.h`
- `symengine/polys/groebner_mogvw.h`
- `symengine/polys/groebner_m4gb.h`
- `symengine/polys/groebner_signature.h`
- new Groebner tests and `benchmarks/groebner_bench.cpp`

## Fixes Applied

### GF(p) interreduction used a modulus-zero domain

`interreduce()` compared coefficients using a default-constructed `Domain`. For `GFpCoeffDomain`, that meant `p_ == 0`, so modular subtraction could silently miscompare coefficients or underflow.

Fixed by comparing through the polynomial's runtime domain (`a.dom`) instead of a default domain.

### Modular precondition failures escaped as exceptions

Modular mode with symbolic coefficients or rational denominators divisible by the modulus could throw `SymEngineException` instead of returning `GroebnerStatus::UnsupportedCoefficientDomain`.

Fixed by adding an internal `GroebnerUnsupportedDomain` sentinel and returning the public unsupported status for:

- symbolic coefficients with `options.modulus > 0`
- non-prime/composite/too-large modulus
- rational coefficient denominators that vanish modulo `p`
- non-rational coefficients in modular mode

### M4GB did not propagate options into the polymatrix

`groebner_basis_m4gb_impl()` copied options into `ctx.options` but not into `ctx.pm.options`, so cancellation behavior inside `get_u_g()` used default options rather than the caller's token/interval.

Fixed by assigning `ctx.pm.options = options`.

### Cancellation interval zero was unsafe

The new matrix/M4GB code used `% options.cancellation_check_interval` and a bit-mask expression derived from it. If a caller set `cancellation_check_interval = 0`, that could divide by zero or underflow.

Fixed by adding `should_check_cancellation()` and using it in the matrix reducer and M4GB `get_u_g()`.

### F5B ignored public resource limits

Explicit F5B calls incremented S-pair and F5-reduction counters but did not honor `max_s_pairs` or `max_reduction_steps`.

Fixed by throwing `GroebnerLimitExceeded` when those limits are reached, producing `GroebnerStatus::ResourceLimitExceeded`.

### M4GB leaked debug output

`insert_basis_polynomial()` and the M4GB tests printed basis internals during normal test runs.

Removed those prints.

### Stats were incomplete for new algorithms

F5B/MoGVW/M4GB could return successful bases with unset `input_polys`, `output_polys`, or matrix stats.

Fixed by normalizing success stats in `groebner_basis()` and by populating matrix/polymatrix stats in the reducer and M4GB path.

### Benchmarks did not expose new counters

`benchmarks/groebner_bench.cpp` still printed only Buchberger-era counters.

Extended benchmark output with signature counters, matrix counters, and polymatrix entry counts. Added benchmark rows for F5B, MoGVW, and a small GF(31) M4GB Katsura-3 run.

## Test Coverage Added

Added coverage for:

- F5B resource-limit behavior
- composite modulus rejection
- symbolic modular coefficient rejection
- rational coefficient denominator rejection modulo `p`
- M4GB explicit precondition failures
- M4GB `cancellation_check_interval = 0`
- M4GB parity with modular Buchberger on the small GF(5) fixture
- M4GB parity with both F5B and modular Buchberger on Katsura-3 over GF(31)

## Remaining Review Findings / Risks

These are not fixed in this pass and should be treated as follow-up work:

- `GroebnerAlgorithm::Auto` still resolves to Buchberger only. The heuristic table from the plan is not implemented.
- `GroebnerOptions::max_milliseconds`, `max_degree`, and `verify_with_buchberger` remain effectively unused.
- M4GB chooses `max_codec_deg = max_input_lm_degree + 30` internally. If a computation needs higher-degree ordinals, it can still fail rather than returning a controlled public status.
- The M4GB implementation is still a small/runtime port with a dense reducer and eager polymatrix entries. It is not yet the full allocation-conscious implementation described in the plan.
- Public `is_groebner()` and `is_reduced_basis()` validate over QQ/Expression, not GF(p). Modular tests therefore use algorithm parity checks rather than a true public modular Groebner predicate.
- MoGVW is scalar-only and not independently cross-checked against a reference implementation beyond canonical basis tests.

## Verification Performed

Commands run from `/work`:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target test_groebner groebner_bench -j2
./build/symengine/tests/polynomial/test_groebner
./build/benchmarks/groebner_bench 1
cmake --build build -j2
ctest --output-on-failure
```

Results:

- `test_groebner`: 172 assertions in 30 test cases passed.
- `groebner_bench 1`: completed successfully and printed populated F5B/MoGVW/M4GB counters.
- Full build completed.
- Full CTest suite: 60/60 tests passed.

The full build still emits existing warnings from `mp_wrapper`/basic arithmetic test code; no warnings remained from the touched Groebner files after the final patch.
