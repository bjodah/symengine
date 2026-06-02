# Implementation Plan: Remaining Groebner Work

## Purpose

This plan covers the follow-up work identified in
`05-CODE-REVIEW-REPORT.md` after the initial F5B, MoGVW, and M4GB pass. It
also includes additional improvements visible from the current implementation
that should be addressed before treating the Groebner branch as production
ready.

The remaining work falls into five themes:

- make public option semantics real and consistent across algorithms;
- implement a meaningful `GroebnerAlgorithm::Auto` resolver;
- turn M4GB from a small runtime port into a controlled, bounded algorithm;
- add true modular validation helpers instead of relying only on parity checks;
- strengthen independent correctness, performance, and maintainability coverage.

## Current State Summary

The current branch has working dispatch for explicit `Buchberger`, `F5B`,
`MoGVW`, and `M4GB`, with finite-field arithmetic and basic tests in place.
The latest review pass also fixed several correctness and API-status issues:
modular domain error handling, F5B resource limits, M4GB option propagation,
cancellation interval zero handling, debug-output leakage, and missing stats.

Remaining risks from the review:

- `GroebnerAlgorithm::Auto` always resolves to Buchberger.
- `GroebnerOptions::max_milliseconds`, `max_degree`, and
  `verify_with_buchberger` are declared but effectively unused.
- M4GB internally chooses `max_codec_deg = max_input_lm_degree + 30`; degree
  overflow can still escape as an internal failure instead of a public status.
- M4GB is still a dense, eager, small-runtime port rather than the
  allocation-conscious implementation described in
  `04-IMPROVED-IMPL-PLAN-MOGVW-AND-M4GB.md`.
- Public `is_groebner()` and `is_reduced_basis()` validate over QQ/Expression
  only, not GF(p).
- MoGVW is scalar-only and is not independently cross-checked beyond canonical
  basis tests.

Additional issues visible in the implementation:

- `groebner_signature.h` is only a forwarding include while shared signature
  types live in `groebner_internal.h`; the boundary between common internals
  and algorithm-specific internals should be clarified.
- `groebner_m4gb.h` still includes development-oriented headers such as
  `<iostream>` and `<cassert>` even though normal debug prints were removed.
- The public header documents option fields, but does not clearly state which
  algorithms honor each option and which preconditions cause
  `UnsupportedCoefficientDomain`.
- The modular primality check in `groebner_basis()` is a naive trial-division
  lambda, which is acceptable for tests but can become unexpectedly expensive
  near the `2^31` limit.
- The shared matrix reducer is dense and classical Gaussian elimination only;
  it does not yet support signature-constrained one-sided elimination.

## Guiding Rules

- Explicit algorithm requests must never silently fall back to another
  algorithm. If preconditions fail, return a controlled public
  `GroebnerStatus`.
- `Auto` may fall back internally, but the selected algorithm should be
  deterministic, documented, and visible through tests or stats.
- Public options must either be honored or documented as unsupported. Do not
  leave API fields that appear functional but have no effect.
- Algorithm verification should prefer independent predicates or reference
  cross-checks over comparing only output formatting.
- Keep changes private to `symengine/polys/groebner*`, tests, benchmarks, and
  documentation unless a wider API change is explicitly needed.

## Phase 1: Option Semantics and Public Contracts

### 1.1 Document Algorithm Preconditions

Files:

- `symengine/polys/groebner.h`
- tests under `symengine/tests/polynomial/`

Implementation tasks:

- Add concise comments near `GroebnerAlgorithm` and `GroebnerOptions` stating:
  - `Buchberger` supports QQ/rational coefficients and symbolic coefficient
    expressions; with `modulus > 0`, it supports GF(p).
  - `F5B` supports the same coefficient domains as Buchberger in the current
    implementation.
  - `MoGVW` currently supports rational and GF(p) scalar coefficients, with
    `DegRevLex` and `GrLex`; symbolic coefficients and `Lex` are unsupported.
  - `M4GB` currently supports only GF(p), `DegRevLex`, and prime
    `2 <= p < 2^31`.
  - `Auto` chooses a concrete algorithm by heuristic and may fall back.
- Document that `cancellation_check_interval == 0` disables interval-based
  checks, while direct loop-entry checks can still observe cancellation.
- Document which resource controls apply to which engines:
  - `max_s_pairs`: Buchberger and F5B now; M4GB critical-pair batches in Phase
    3; MoGVW lifted-monomial work in Phase 2.
  - `max_reduction_steps`: Buchberger, F5B, normal forms; MoGVW top reductions
    in Phase 2.
  - `max_matrix_rows` / `max_matrix_columns`: matrix reducer and M4GB.
  - `max_degree`: all algorithms after Phase 2.
  - `max_milliseconds`: all algorithms after Phase 2.

Acceptance criteria:

- Header comments match actual implementation behavior.
- Tests exist for each explicit algorithm's unsupported precondition path.

### 1.2 Replace Ad Hoc Primality Check

Files:

- `symengine/polys/groebner.cpp`
- optionally `symengine/polys/groebner_internal.h`

Implementation tasks:

- Move modulus validation into a small helper, for example:
  `bool is_supported_prime_modulus(uint64_t p)`.
- Replace trial division by a deterministic Miller-Rabin implementation for
  32-bit inputs, or use a bounded optimized trial division with documented
  worst-case behavior. Since the valid range is `< 2^31`, deterministic bases
  for 32-bit integers are sufficient.
- Return `UnsupportedCoefficientDomain` for `p < 2`, composite `p`, and
  `p >= 2^31`.

Acceptance criteria:

- Existing composite and too-large modulus tests still pass.
- Add tests for boundary values `0`, `1`, `2`, `2147483647`, and
  `2147483648`.
- No algorithm body has to duplicate modulus validation.

## Phase 2: Honor Remaining Public Options

### 2.1 Implement `max_milliseconds`

Files:

- `symengine/polys/groebner.h`
- `symengine/polys/groebner_internal.h`
- `symengine/polys/groebner.cpp`
- `symengine/polys/groebner_f5b.h`
- `symengine/polys/groebner_mogvw.h`
- `symengine/polys/groebner_m4gb.h`
- FGLM code paths in `groebner.cpp`

Implementation tasks:

- Add a private runtime context that stores:
  - the user options;
  - a `std::chrono::steady_clock::time_point` start time;
  - a helper that checks cancellation token and elapsed milliseconds together.
- Preserve the public `GroebnerOptions` type; do not add mutable timing state
  to it.
- Replace direct `check_cancellation(options)` calls in long-running loops with
  a helper that also enforces `max_milliseconds`.
- Keep exception translation as `GroebnerStatus::ResourceLimitExceeded` for
  elapsed-time exhaustion. This matches the existing resource-limit status.
- Ensure `cancellation_token` still takes precedence when both cancellation and
  time limit are hit.

Acceptance criteria:

- A tiny `max_milliseconds = 1` or similarly low test can force
  `ResourceLimitExceeded` on a workload that normally succeeds.
- Cancellation tests still return `Cancelled`.
- Buchberger, F5B, MoGVW, M4GB, normal form, interreduce, and FGLM paths all
  share the same timeout behavior.

### 2.2 Implement `max_degree`

Files:

- `symengine/polys/groebner_internal.h`
- `symengine/polys/groebner.cpp`
- `symengine/polys/groebner_f5b.h`
- `symengine/polys/groebner_mogvw.h`
- `symengine/polys/groebner_m4gb.h`

Implementation tasks:

- Define one helper:
  `void enforce_max_degree(unsigned degree, const GroebnerOptions &options)`.
- Call it before:
  - creating or processing Buchberger/F5B critical pairs by LCM degree;
  - accepting a nonzero reduced polynomial into a basis;
  - MoGVW variable lifting and post-top-reduction insertion;
  - M4GB codec construction, `increase_upper_bound`, `selection`, and
    `get_u_g` monomial multiplication.
- Treat degree exhaustion as `GroebnerLimitExceeded`, producing
  `GroebnerStatus::ResourceLimitExceeded`.
- For M4GB, use `options.max_degree` as the codec degree when nonzero and
  reject any computation that exceeds it instead of silently extending by an
  internal cushion.

Acceptance criteria:

- `max_degree = 0` preserves current behavior.
- Tests show each algorithm returns `ResourceLimitExceeded` when the bound is
  too small.
- M4GB no longer throws a raw runtime error for `degree > D`.

### 2.3 Implement `verify_with_buchberger`

Files:

- `symengine/polys/groebner.cpp`
- `symengine/polys/groebner.h`
- tests under `symengine/tests/polynomial/`

Implementation tasks:

- When `options.verify_with_buchberger` is true and the selected algorithm is
  not Buchberger:
  - run the selected algorithm normally;
  - run Buchberger with the same variables, order, reduction, sorting, modulus,
    and resource/cancellation options, but with `verify_with_buchberger = false`;
  - compare both reduced bases after canonical sorting.
- Use true modular comparison once Phase 5 adds modular predicates. Until then,
  for GF(p), compare canonical output vectors produced by the same public
  conversion path.
- Avoid recursive verification when Buchberger itself is selected.
- Decide the public failure status:
  - preferred: add `GroebnerStatus::VerificationFailed` if API expansion is
    acceptable;
  - fallback: throw/translate to `UnsupportedCoefficientDomain` only if the API
    cannot be extended. This is less precise and should be documented.

Acceptance criteria:

- Verification is off by default.
- Tests cover success verification for F5B, MoGVW, and M4GB.
- A test-only injected mismatch path, if feasible, demonstrates the selected
  failure status.

### 2.4 Extend Limits to MoGVW Work Units

Files:

- `symengine/polys/groebner_mogvw.h`

Implementation tasks:

- Increment a consistent work counter for MoGVW lifted monomials and top
  reductions.
- Honor `max_s_pairs` as a cap on lifted monomial candidates, or document a
  separate interpretation if the name is too Buchberger-specific.
- Honor `max_reduction_steps` inside `top_reduce()`.
- Update stats so rejected candidates and reduction steps are visible.

Acceptance criteria:

- Tests can force `ResourceLimitExceeded` in MoGVW with both a low candidate
  limit and a low reduction-step limit.
- Existing MoGVW correctness tests remain green.

## Phase 3: Implement Auto Algorithm Resolution

Files:

- `symengine/polys/groebner.cpp`
- `symengine/polys/groebner.h`
- benchmarks and tests

Implementation tasks:

- Replace the current hardcoded Buchberger return in `resolve_algorithm()`.
- Compute cheap input features:
  - number of variables;
  - number of input polynomials;
  - total number of terms;
  - max input total degree;
  - coefficient class: symbolic, rational, or GF(p);
  - monomial order;
  - whether a finite-field modulus is present.
- Implement conservative heuristics:
  - symbolic coefficients: Buchberger.
  - explicit modulus with `DegRevLex`, no symbolic coefficients, and moderate
    variable/degree counts: M4GB if bounds are compatible.
  - rational/GF(p) scalar systems with `DegRevLex` or `GrLex`: F5B for small
    systems where signature overhead is acceptable; MoGVW only when its scalar
    implementation is known to be competitive on the benchmark corpus.
  - `Lex`: Buchberger or F5B; do not choose M4GB or MoGVW until support exists.
  - If an Auto-selected algorithm fails preconditions before computation,
    fall back to Buchberger rather than surfacing unsupported status.
- Keep explicit requests strict: this fallback applies only to `Auto`.
- Add an internal selected-algorithm variable visible in tests. Preferred:
  add a `selected_algorithm` field to `GroebnerStats`; fallback: expose the
  decision through deterministic test cases and algorithm-specific counters.

Suggested initial heuristic table:

| Input shape | Auto choice |
| --- | --- |
| symbolic coefficients | Buchberger |
| GF(p), `DegRevLex`, variables <= 8, max input degree <= effective max degree | M4GB |
| GF(p), `GrLex`, non-symbolic | MoGVW |
| rational, `DegRevLex`/`GrLex`, polynomials >= 3 | F5B |
| rational, `Lex` | F5B for small systems, otherwise Buchberger |
| unsupported or uncertain | Buchberger |

Acceptance criteria:

- Auto no longer always selects Buchberger.
- Tests assert representative decisions for symbolic, rational, GF(p), Lex,
  GrLex, and DegRevLex inputs.
- Explicit algorithm behavior remains unchanged.
- Benchmarks include `Auto` rows so heuristic regressions are visible.

## Phase 4: M4GB Correctness, Bounds, and Performance

### 4.1 Replace `max_input_lm_degree + 30`

Files:

- `symengine/polys/groebner_m4gb.h`
- `symengine/polys/groebner_internal.h`

Implementation tasks:

- Remove the internal fixed cushion:
  `unsigned max_codec_deg = max_deg + 30`.
- Define an effective codec degree:
  - if `options.max_degree > 0`, use that exact value;
  - otherwise choose an adaptive starting degree based on input degree and
    grow the codec when required, or compute a documented conservative bound
    for the selected batch.
- Convert codec overflow, binomial overflow, and `to_index()` degree overflow
  into `GroebnerLimitExceeded` or `GroebnerUnsupportedDomain` as appropriate.
- Prefer `ResourceLimitExceeded` when the user's degree/resource bound is hit.
  Prefer `UnsupportedCoefficientDomain` when the monomial universe cannot be
  represented in `uint64_t`.

Acceptance criteria:

- No raw `std::runtime_error("degree > D!")` can escape public
  `groebner_basis()`.
- Tests cover degree overflow and explicit `max_degree` exhaustion.
- M4GB parity tests still pass.

### 4.2 Make Polymatrix Growth Bounded

Files:

- `symengine/polys/groebner_m4gb.h`
- `symengine/polys/groebner_internal.h`

Implementation tasks:

- Enforce `max_matrix_rows` and `max_matrix_columns` before dense matrix
  materialization and during polymatrix growth, not only inside the reducer.
- Track and check:
  - number of dense-index columns;
  - number of matrix entries;
  - selected batch rows;
  - generated multiples from `increase_upper_bound()` and `get_u_g()`.
- Use `matrix_batch_size` consistently to cap selected critical pairs with the
  same LCM degree.
- Add stats for rejected/deferred pairs if needed.

Acceptance criteria:

- Low row/column limits reliably return `ResourceLimitExceeded` before large
  allocations.
- Stats show maximum rows and columns even on limit failure when possible.
- `matrix_batch_size` has a regression test that changes the number of
  matrices built or batches processed.

### 4.3 Reduce Allocation Pressure

Files:

- `symengine/polys/groebner_m4gb.h`
- `symengine/polys/groebner_internal.h`

Implementation tasks:

- Avoid repeated sparse/dense conversions where a row can be updated in place.
- Reuse temporary exponent vectors in `RuntimeDegRevLexCodec::multiply`,
  `divide`, `divides`, and `lcm`, or add overloads that operate on caller-owned
  buffers.
- Reserve dense row capacity based on known `dense_index` size.
- Avoid repeated full scans of `pm.basis` where a leading-monomial index can
  answer divisor queries.
- Remove unnecessary includes such as `<iostream>` from M4GB and internal
  headers.

Acceptance criteria:

- Benchmarks show no regression for existing small M4GB fixtures.
- Allocation-sensitive benchmarks, if available through heap profiling, show
  fewer allocations in `get_u_g()`, `reduce_sparse()`, and codec operations.
- Code remains deterministic and test results are unchanged.

### 4.4 Add M4GB Invariant Tests

Files:

- `symengine/tests/polynomial/test_groebner_m4gb.cpp`
- optionally a dedicated internal test file

Implementation tasks:

- Test codec round-trips for every monomial up to small `(num_vars, degree)`
  pairs.
- Test ordinal ordering agrees with `compare_monomial_less()` for DegRevLex.
- Test `multiply`, `divide`, `divides`, and `lcm` against exponent-vector
  operations.
- Add targeted polymatrix tests for:
  - dense index/inverse index consistency;
  - matrix entries referencing valid basis leading monomials;
  - shrink preserving sparse polynomial values;
  - no duplicate or stale basis leading monomials after update.

Acceptance criteria:

- M4GB internals have property-style tests independent of final basis parity.
- Tests run as part of `test_groebner`.

## Phase 5: Modular Public Predicates

Files:

- `symengine/polys/groebner.h`
- `symengine/polys/groebner.cpp`
- `symengine/tests/polynomial/test_groebner.cpp`
- `symengine/tests/polynomial/test_groebner_m4gb.cpp`

Implementation tasks:

- Add overloads or option-aware variants for public predicates:
  - `normal_form(..., const GroebnerOptions &options)`;
  - `is_groebner(..., const GroebnerOptions &options)`;
  - `is_reduced_basis(..., const GroebnerOptions &options)`.
- Preserve existing overloads as QQ/Expression-compatible defaults.
- When `options.modulus > 0`, parse coefficients through `GFpCoeffDomain`
  using the same conversion and unsupported-domain rules as `groebner_basis()`.
- Ensure modular monic checks compare against GF(p) one, not rational one.
- Use the modular predicate in M4GB tests instead of only comparing against
  Buchberger/F5B output.

Acceptance criteria:

- Public modular `is_groebner()` and `is_reduced_basis()` tests pass for GF(5)
  and GF(31) fixtures.
- Tests include a basis that is valid over QQ but not over a selected GF(p), or
  vice versa, to prove the modulus is applied.
- Existing public predicate overloads remain source-compatible.

## Phase 6: MoGVW Independent Validation and Scope Decision

Files:

- `symengine/polys/groebner_mogvw.h`
- `symengine/tests/polynomial/test_groebner_mogvw.cpp`
- benchmarks

Implementation tasks:

- Port small reference examples from the MoGVW paper where feasible.
- Add randomized small scalar systems over QQ and GF(p), comparing:
  - MoGVW reduced basis against Buchberger;
  - MoGVW reduced basis against F5B where F5B supports the same input;
  - public modular predicates after Phase 5.
- Add tests for each rejection criterion:
  - LCM;
  - syzygy;
  - rewritten;
  - collision resolution.
- Decide whether the current scalar-only MoGVW should remain an exposed public
  algorithm or be documented as experimental until the full monomial-oriented
  implementation is complete.

Acceptance criteria:

- MoGVW has correctness tests that fail if each major criterion is disabled.
- Benchmark output includes MoGVW on at least one case where it is expected to
  win or remain competitive.
- Header documentation states the current MoGVW scope.

## Phase 7: Shared Internal Cleanup

Files:

- `symengine/polys/groebner_internal.h`
- `symengine/polys/groebner_signature.h`
- `symengine/polys/groebner_f5b.h`
- `symengine/polys/groebner_mogvw.h`
- `symengine/polys/groebner_m4gb.h`

Implementation tasks:

- Move `Signature`, `SignatureLess`, and `signature_mul()` into
  `groebner_signature.h`, or remove the forwarding header if the team prefers
  one internal header. Avoid the current ambiguous split.
- Move modular-domain helpers and modulus validation into a cohesive internal
  section.
- Remove unused development includes from Groebner headers.
- Audit for raw `std::runtime_error` throws inside public algorithm paths and
  translate them into internal sentinel exceptions.
- Keep comments that describe invariants; remove comments that restate obvious
  field names.

Acceptance criteria:

- No public algorithm path leaks a raw C++ exception for expected input/domain
  failures.
- Header dependencies remain minimal.
- Full build produces no new warnings from touched Groebner files.

## Phase 8: Benchmarks and Performance Tracking

Files:

- `benchmarks/groebner_bench.cpp`
- optional benchmark documentation

Implementation tasks:

- Add `Auto` benchmark rows for existing fixtures.
- Print selected algorithm if a `GroebnerStats::selected_algorithm` field is
  added.
- Add benchmark cases that exercise:
  - rational F5B on small and medium systems;
  - MoGVW scalar systems in `GrLex` and `DegRevLex`;
  - M4GB over multiple primes, including GF(31) and a larger valid prime;
  - low `matrix_batch_size` versus default batch size.
- Add optional benchmark verification mode using `verify_with_buchberger`.
- Keep benchmark output compact and stable for comparison.

Acceptance criteria:

- `./build/benchmarks/groebner_bench 1` completes and prints all new counters.
- Benchmark rows make Auto decisions and option effects visible.

## Phase 9: Verification Matrix

Run these after each phase that changes behavior:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target test_groebner groebner_bench -j2
./build/symengine/tests/polynomial/test_groebner
./build/benchmarks/groebner_bench 1
```

Run these before declaring the branch complete:

```bash
cmake --build build -j2
ctest --output-on-failure
```

Additional targeted verification:

- Run debug builds for invariant assertions once M4GB invariant tests are added.
- Run sanitizer builds if the local toolchain supports them, especially after
  M4GB allocation and codec changes.
- Compare benchmark output before and after Auto/M4GB changes to detect silent
  performance regressions.

## Suggested Implementation Order

1. Phase 1: document contracts and centralize modulus validation.
2. Phase 2: implement `max_milliseconds`, `max_degree`,
   `verify_with_buchberger`, and MoGVW limits.
3. Phase 4.1: fix M4GB codec-degree overflow and public status translation.
4. Phase 5: add modular public predicates.
5. Phase 3: implement `Auto` once algorithm behavior and validation are more
   stable.
6. Phase 4.2-4.4: make M4GB bounded, faster, and independently tested.
7. Phase 6: strengthen MoGVW validation and decide whether to mark it
   experimental.
8. Phase 7: internal cleanup after behavior is stable.
9. Phase 8 and final Phase 9 verification.

This order avoids making `Auto` prefer algorithms whose limits and validation
semantics are still incomplete.
