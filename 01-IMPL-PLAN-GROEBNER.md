# Standalone Groebner Basis Implementation Plan

## Summary

Add native SymEngine support for Groebner basis computation and the order-conversion workflow needed for polynomial system solving:

1. Convert input equations into a sparse multivariate polynomial representation over a field.
2. Compute a Groebner basis, typically in degree reverse lexicographic order.
3. For zero-dimensional ideals, convert the basis to lexicographic order with FGLM.
4. Use the lex basis as a triangular system for univariate solving and back-substitution.

Singular under `./singular-4.4.1+ds/` is reference material only. The implementation must not link to Singular, shell out to Singular, wrap Singular APIs, or introduce a runtime dependency on Singular.

## Repository Investigation

### Existing SymEngine Polynomial Support

Relevant files:

- `symengine/polys/msymenginepoly.h`
- `symengine/polys/msymenginepoly.cpp`
- `symengine/polys/basic_conversions.h`
- `symengine/rings.cpp`
- `symengine/solve.cpp`
- `symengine/tests/polynomial/*`

Findings:

- SymEngine has native multivariate sparse polynomial wrappers:
  - `MIntPoly`: multivariate integer coefficients.
  - `MExprPoly`: multivariate expression coefficients.
- There is no dedicated multivariate rational polynomial type comparable to `URatPoly`.
- The existing multivariate containers are unordered maps from exponent vectors to coefficients. That is good for basic arithmetic, but not ideal for Groebner algorithms, which repeatedly need leading monomial/leading coefficient queries under a chosen monomial order.
- `basic_conversions.h` can convert `Basic` expressions into `MIntPoly` and `MExprPoly` with explicit generators. This is useful for API-level conversion but should not dictate the internal Groebner representation.
- `rings.cpp` has older sparse integer polynomial helpers, but they are too limited for Groebner work.
- `solve.cpp` currently handles:
  - univariate polynomial solving up to quartic by formulas;
  - optional FLINT factorization for some higher-degree univariate cases;
  - linear systems through matrix methods.
- There is no nonlinear polynomial-system solver and no Groebner API.

### Build and Test Layout

Relevant files:

- `symengine/CMakeLists.txt`
- `symengine/tests/polynomial/CMakeLists.txt`
- `symengine/tests/basic/test_solve.cpp`

Findings:

- New core code should be added to `symengine/CMakeLists.txt` `SRC` and `HEADERS`.
- Unit tests for polynomial internals belong under `symengine/tests/polynomial/`.
- Solver-level tests belong under `symengine/tests/basic/test_solve.cpp` or a new solver test if the surface grows enough.
- The project targets C++11 by default. Avoid C++17 APIs unless gated behind existing LLVM-related standard changes.

### `symcse` Prototype References

Relevant files:

- `symcse/symcse/polynomials/_groebner.py`
- `symcse/symcse/polynomials/_groebner_sage.py`
- `symcse/symcse/polynomials/_singular_groebner.py`
- `symcse/symcse/polynomials/mpoly.py`
- `symcse/symcse/polynomials/mpolysys.py`
- `symcse/share/polynomials/_singular_groebner_postprocess.cpp`
- `symcse/tests/polynomials/test__groebner.py`
- `symcse/tests/polynomials/test__groebner_sage.py`
- `symcse/tests/polynomials/test__singular_groebner.py`

Findings:

- `symcse` treats unknowns and parameters separately. This matches the desired `C1`, `C2`, ... workflow: variables are explicit, and all other symbols are coefficients.
- Its preferred expensive path is `degrevlex` computation followed by FGLM conversion to `lex`.
- It uses an augmented variable `zaug` with an equation `zaug - z[i]` to force a chosen solved variable into the triangular output. This is important for practical back-substitution because it lets the caller choose which original variable becomes the univariate root variable.
- Its post-processing assumes a lex basis shape consisting of one univariate polynomial plus remaining equations that are linear in the other variables over the chosen root. SymEngine should support this as a first practical triangular-solving mode before attempting fully general algebraic back-substitution.
- It applies explicit resource controls for external computations (`ulimit -v`, `ulimit -t`). In native SymEngine these should become cooperative cancellation plus optional work/time/memory budgets.
- It has tests that are directly useful as future acceptance scenarios:
  - `{x**2 - y, x**3 - x}` should produce a lex basis containing `x**2 - y`, `x*y - x`, and `y**2 - y`.
  - `{x**2 + K1*x*y, x*y + 2*y**3 - K2}` tests coefficient parameters.
  - chemical-equilibrium-like systems test the augmented-variable flow and numeric residual validation.
- `symcse` has prototypes for C++ F5 and moGVW command-line tools, but only as external experimental calls. This plan should keep F5/moGVW as later native algorithm options, not as dependencies.

### Singular as Algorithm Reference Only

Relevant Singular reference files:

- `singular-4.4.1+ds/kernel/GBEngine/kstd1.h`
- `singular-4.4.1+ds/kernel/GBEngine/kstd2.cc`
- `singular-4.4.1+ds/kernel/fglm/fglm.h`
- `singular-4.4.1+ds/libpolys/polys/monomials/ring.h`
- `singular-4.4.1+ds/libpolys/misc/options.h`

Useful observations:

- Singular separates ring/order setup, polynomial arithmetic, standard basis computation, and FGLM conversion.
- Its standard basis engine has many strategies and criteria. SymEngine should start with a smaller, maintainable Buchberger implementation, then add F4-style linear algebra once correctness is established.
- Singular exposes an interrupt option around long Groebner computations. SymEngine should provide a native cooperative cancellation mechanism rather than using global process signals internally.

## Curated Test Cases From Existing Suites

Sources inspected:

- SymPy installed tests:
  - `/opt-3/cpython-v3.13-apt-deb/lib/python3.13/site-packages/sympy/polys/tests/test_orderings.py`
  - `/opt-3/cpython-v3.13-apt-deb/lib/python3.13/site-packages/sympy/polys/tests/test_groebnertools.py`
  - `/opt-3/cpython-v3.13-apt-deb/lib/python3.13/site-packages/sympy/polys/tests/test_polytools.py`
  - `/opt-3/cpython-v3.13-apt-deb/lib/python3.13/site-packages/sympy/solvers/tests/test_polysys.py`
- `symcse/tests/polynomials/test__groebner.py`
- `symcse/tests/polynomials/test__groebner_sage.py`
- `symcse/tests/polynomials/test__singular_groebner.py`
- Singular reference examples:
  - `singular-4.4.1+ds/kernel/fglm/test.cc`
  - `singular-4.4.1+ds/kernel/combinatorics/test.cc`

These should be copied into SymEngine tests selectively. Exact-output tests are good for early correctness; property tests are better where ordering, reduction normalization, or symbolic-parameter simplification could reasonably differ.

### Monomial Ordering Tests

Add these in Phase 0 before implementing any Groebner algorithm:

- `lex((1, 2, 3)) == (1, 2, 3)`.
- `grlex((1, 2, 3))` compares as total degree `6` followed by `(1, 2, 3)`.
- `degrevlex((1, 2, 3))` compares as total degree `6` followed by reverse negative exponents `(-3, -2, -1)`.
- Equal-degree comparison checks:
  - `(0, 1, 1) > (0, 0, 2)` under `grlex` and `degrevlex`.
  - `(0, 3, 1) < (2, 2, 1)` under `grlex` and `degrevlex`.
- With generator order `[z, y, x]`, sort
  `{x, x^2*z^2, x*y, x^2, 1, y^2, x^3, y, z, x*y^2*z, x^2*y^2}` and verify:
  - `lex`: `{1, x, x^2, x^3, y, x*y, y^2, x^2*y^2, z, x*y^2*z, x^2*z^2}`;
  - `grlex`: `{1, x, y, z, x^2, x*y, y^2, x^3, x^2*y^2, x*y^2*z, x^2*z^2}`;
  - `degrevlex`: `{1, x, y, z, x^2, x*y, y^2, x^3, x^2*y^2, x^2*z^2, x*y^2*z}`.

### Phase 1 Exact Groebner Basis Tests

Small rational exact-output cases from SymPy:

- Over `QQ[x, y]`, `lex`:
  - input `{x^2 + 2*x*y^2, x*y + 2*y^3 - 1}`;
  - expected `{x, y^3 - 1/2}`.
- Over `QQ[y, x]`, `lex`:
  - input `{2*x^2*y + y^2, 2*x^3 + x*y - 1}`;
  - expected `{y, x^3 - 1/2}`.
- Over `QQ[x, y]`, `grlex`:
  - input `{x^3 - 2*x*y, x^2*y + x - 2*y^2}`;
  - expected `{x^2, x*y, -1/2*x + y^2}`.
- Over `QQ[x, y, z]`, `lex`:
  - input `{-x^2 + y, -x^3 + z}`;
  - expected `{x^2 - y, x*y - z, x*z - y^2, y^3 - z^2}`.
- Same input over `grlex`:
  - expected `{y^3 - z^2, x^2 - y, x*y - z, x*z - y^2}`.
- Over `QQ[x, y, z]`, `lex`:
  - input `{x - y^2, -y^3 + z}`;
  - expected `{x - y^2, y^3 - z}`.
- Same input over `grlex`:
  - expected `{x^2 - y*z, x*y - z, -x + y^2}`.
- Over `QQ[x, y, z]`, `lex`:
  - input `{x - z^2, y - z^3}`;
  - expected `{x - z^2, y - z^3}`.
- Same input over `grlex`:
  - expected `{x^2 - y*z, x*z - y, -x + z^2}`.
- Over `QQ[x, y, z]`, `lex`:
  - input `{-y^2 + z, x - y^3}`;
  - expected `{x - y*z, y^2 - z}`.
- Same input over `grlex`:
  - expected `{-x^2 + z^3, x*y - z^2, y^2 - z, -x + y*z}`.
- Univariate reduction behavior:
  - `{x^3 - 1, x^2 - 1}` should reduce to `{x - 1}`;
  - `{x^2 - 1, x^3 + 1}` should reduce to `{x + 1}`.
- Unit and empty ideals:
  - empty input returns an empty basis;
  - `{1}` returns `{1}`.

Medium exact-output cases for a second correctness wave:

- Over `QQ[x, y]`, `lex`:
  - input `{4*x^2*y^2 + 4*x*y + 1, x^2 + y^2 - 1}`;
  - expected `{x - 4*y^7 + 8*y^5 - 7*y^3 + 3*y, y^8 - 2*y^6 + 3/2*y^4 - 1/2*y^2 + 1/16}`.
- Katsura-3 over `ZZ` or `QQ`, variables `{x0, x1, x2}`, `lex`:
  - input `{x0 + 2*x1 + 2*x2 - 1, x0^2 + 2*x1^2 + 2*x2^2 - x0, 2*x0*x1 + 2*x1*x2 - x1}`;
  - expected `{ -7 + 7*x0 + 8*x2 + 158*x2^2 - 420*x2^3, 7*x1 + 3*x2 - 79*x2^2 + 210*x2^3, x2 + x2^2 - 40*x2^3 + 84*x2^4 }`, normalized to the implementation's canonical leading coefficient convention.

Future finite-field/domain test from SymPy:

- Over `GF(7)[x, y, z]`, input `{3*x^2 + y*z - 5*x - 1, 2*x + 3*x*y + y^2, x - 3*y + x*z - 2*z^2}`.
- Expected basis has four generators:
  - `1 + x + y + 3*z + 2*z^2 + 2*z^3 + 6*z^4 + z^5`;
  - `1 + 3*y + y^2 + 6*z^2 + 3*z^3 + 3*z^4 + 3*z^5 + 4*z^6`;
  - `1 + 4*y + 4*z + y*z + 4*z^3 + z^4 + z^6`;
  - `6 + 6*z + z^2 + 4*z^3 + 3*z^4 + 6*z^5 + 3*z^6 + z^7`.
- Use this initially as a skipped fixture documenting future finite-field behavior.

### Reduction and Basis-Property Tests

These should not require exact output order:

- `normal_form(2*x^4 + y^2 - x^2 + y^3, {x^3 - x, y^3 - y})` should produce quotient `{2*x, 1}` and remainder `x^2 + y^2 + y`.
- For `G = groebner({x^2 + y^2 - 1, x*y - 2})`, reducing `2*x^3 + y^3 + 3*y` by `G` should have zero remainder.
- `contains` behavior:
  - the same `G` contains `2*x^3 + y^3 + 3*y`;
  - it does not contain that polynomial plus `1`.
- `is_groebner` behavior:
  - `{x^2, x*y, -1/2*x + y^2}` is a Groebner basis under the tested ring/order;
  - `{x^3, x*y, -1/2*x + y^2}` is not.
- `is_reduced` behavior:
  - `{2*x^2 + y^2, x^2 + y^2}` is not reduced;
  - the computed reduced basis for that ideal is reduced.
- Adding zero generators should not change the computed basis.

### FGLM and Zero-Dimensional Tests

Zero-dimensional detection:

- `{x, y}` over `{x, y}` is zero-dimensional.
- `{x^3 + y^2}` over `{x, y}` is not zero-dimensional.
- `{x, y, z}` over `{x, y, z}` is zero-dimensional.
- `{x, y, z}` over `{x, y, z, t}` is not zero-dimensional.
- `{x*y - z, y*z - x, x*y - y}` over `{x, y, z}` is zero-dimensional.
- `{x^2 - 2*x*z + 5, x*y^2 + y*z^3, 3*y^2 - 8*z^2}` over `{x, y, z}` is zero-dimensional.

FGLM exact-output tests:

- For `F = {x^2 - x - 3*y + 1, -2*x + y^2 + y - 1}`:
  - compute a `lex` basis, then convert to `grlex`;
  - expected `grlex` basis is `{x^2 - x - 3*y + 1, y^2 - 2*x + y - 1}`.
- For `F = {a + b + c + d, a*b + a*d + b*c + b*d, a*b*c + a*b*d + a*c*d + b*c*d, a*b*c*d - 1}`:
  - compute in `grlex`, convert to `lex`;
  - expected lex basis:
    - `4*a + 3*d^9 - 4*d^5 - 3*d`;
    - `4*b + 4*c - 3*d^9 + 4*d^5 + 7*d`;
    - `4*c^2 + 3*d^10 - 4*d^6 - 3*d^2`;
    - `4*c*d^4 + 4*c - d^9 + 4*d^5 + 5*d`;
    - `d^12 - d^8 - d^4 + 1`.
- Include a later large-coefficient FGLM stress test from SymPy's `test_fglm` for the two-variable `{t, x}` system with degree-8 univariate polynomial. This should be marked slow or benchmark-only until linear algebra performance is acceptable.

### Polynomial System Solving Tests

Use these after the Groebner and FGLM layers exist:

- Single variable: `{x - 1}` over `{x}` gives `{(1)}`.
- Inconsistent system: `{y - x, y - x - 1}` gives no solutions.
- Shared zero: `{y - x^2, y + x^2}` gives `{(0, 0)}`.
- Linear 3-variable system:
  - `{2*x - 3, 3/2*y - 2*x, z - 5*y}` gives `{(3/2, 2, 10)}`.
- Finite nonlinear system:
  - `{x*y - 2*y, 2*y^2 - x^2}` gives `{(0, 0), (2, -sqrt(2)), (2, sqrt(2))}`.
- Complex roots:
  - `{y - x^2, y + x^2 + 1}` gives `{(-I*sqrt(1/2), -1/2), (I*sqrt(1/2), -1/2)}`.
- Symmetric 3-variable system:
  - `{x^2 + y + z - 1, x + y^2 + z - 1, x + y + z^2 - 1}`;
  - exact algebraic solutions are permutations of `(1, 0, 0)` plus `(sqrt(2)-1, sqrt(2)-1, sqrt(2)-1)` and `(-sqrt(2)-1, -sqrt(2)-1, -sqrt(2)-1)`;
  - triangular solving over `QQ` may initially return only the three rational points, with algebraic-extension support adding the two irrational points later.
- Positive-dimensional rejection:
  - `{x^3 - y^3}` over `{x, y}` should fail with a clear not-zero-dimensional or unsupported status;
  - `{x - 1}` over variables `{x, y}` should also be rejected for finite-solution solving.
- Unsolvable high-degree univariate:
  - `{x^5 - x + 1}` should either return a structured `RootOf`/`ConditionSet` style result or a clear unsupported status; it should not be reported as having no mathematical roots.

### Symbolic-Parameter and `symcse` Acceptance Tests

Parameter coefficients:

- Variables `{x, y}`, parameters `{K1, K2}`:
  - input `{x^2 + K1*x*y, x*y + 2*y^3 - K2}`;
  - exact string output can differ by coefficient normalization, but the basis should have two generators and reduce both inputs to zero;
  - the expected Sage/Singular-shaped lex forms include:
    - `K1*K2*y - 2*K1*y^4 + K2*x - 2*K2*y^2 + 4*y^5`;
    - `K1*K2*y^2 - 2*K1*y^5 + K2^2 - 4*K2*y^3 + 4*y^6`.
- Variables `{z, y, x}`, parameters `{K1, ux}`:
  - input `{y*z - K1*(ux - x), y - x, z - x}`;
  - lex basis has three generators in the current `symcse` workflow;
  - solving with `K1 = 1e-14/55.4`, `ux = 55.4 + 1e-7` should produce approximately `(z, y, x) = (1e-7, 1e-7, 1e-7)`.
- Variables `{z, y, x}`, parameters `{K1, b0, b1, ux}`:
  - input `{y*z - K1*(ux - x), (ux - x) + y - b0, (ux - x) + z - b1}`;
  - with `K1 = 1e-14/55.4`, `b0 = b1 = 55.4 + 1e-7`, `ux = min(b0, b1)`, numeric solving should again produce approximately `(1e-7, 1e-7, 1e-7)`.

Augmented-variable workflow:

- For `{x^2 - y, x^3 - x}` over `{x, y}`:
  - direct `lex`, `degrevlex` plus FGLM, and reverse-lex-style input plus FGLM should all yield a basis containing `x^2 - y`, `x*y - x`, and `y^2 - y`.
- Chemical-equilibrium system from `symcse`:
  - variables `{NH4p, NH3, Hp, H2O, OHm}`;
  - parameters `{Ka, Kw, iNH3, iOHm, iHp}`;
  - equations:
    - `NH4p*Ka - NH3*Hp`;
    - `Kw*H2O - OHm*Hp`;
    - `NH4p + NH3 - iNH3`;
    - `H2O + OHm - iOHm`;
    - `H2O + Hp + NH4p - iHp`;
  - for each selected variable, add `zaug - selected_variable`, compute `degrevlex`, convert to `lex`, solve, and verify residuals numerically with `Ka = 10^-9.24`, `Kw = 10^-14`, `iNH3 = 1.0`, `iOHm = 55.4`, `iHp = 55.4`.
- Larger augmented stress case from `symcse/tests/polynomials/test__singular_groebner.py`:
  - variables `{z0, z1, z2, z3, z4, z5, z6, z7, zaug}`;
  - parameters `{K0n1, K1p1, K1m1, K1p2, b0, b1, b2, b3}`;
  - equations:
    - `-z2*z3 + z1*z0*K0n1`;
    - `z4*K1p1 - z5*z3`;
    - `-z6 + z1*z5*K1m1`;
    - `z0*K1p2 - z3*z7`;
    - `-b0 + z0 + z2 + z7`;
    - `-b1 + z4 + z5 + z6`;
    - `-b2 + z0 + z3 + z4`;
    - `-b3 + z1 + z2 + z6`;
    - `-z1 + zaug`;
  - initially require only clean completion/cancellation and residual validation after numeric substitution, not a fixed exact basis.

### Singular-Derived Property and Benchmark Cases

Use these only as standalone mathematical examples; do not compare against a Singular runtime:

- Toric/binomial ideal from Singular kernel tests:
  - variables `{w, x, y, z}`, order `degrevlex`;
  - input `{w^2 - x*z, w*x - y*z, x^2 - w*y, x*y - z^2, y^2 - w*z}`;
  - use as a property test: every input reduces to zero by the computed basis, every S-polynomial of the basis reduces to zero, and the reduced basis is canonical under SymEngine normalization.
- FGLM engine examples in `kernel/fglm/test.cc` should be mined for additional zero-dimensional benchmark systems after the native FGLM implementation exists. The first implementation should prefer the smaller SymPy FGLM exact-output tests above.

## Goals

- Provide a standalone exact Groebner basis implementation in SymEngine.
- Support monomial orders:
  - lexicographic (`lex`);
  - graded lexicographic (`grlex`);
  - degree reverse lexicographic (`degrevlex`).
- Support common workflow:
  - compute basis in `degrevlex`;
  - convert zero-dimensional ideals to `lex` using FGLM;
  - expose enough lex-basis structure for triangular solving.
- Support coefficient domains needed by practical symbolic systems:
  - exact rationals first;
  - symbolic parameter coefficients such as `C1*x*y + C2*x*z**2 - 42*y**3*z` when variables are explicitly `{x, y, z}`;
  - finite fields later for modular methods and performance.
- Provide cancellation that exits computations cleanly through normal C++ control flow.

## Non-Goals

- No Singular bridge.
- No subprocess calls to `Singular`.
- No dependency on Singular headers, libraries, global state, or data structures.
- No first-pass attempt to reproduce Singular's full optimized GB engine.
- No support for noncommutative Groebner bases in the initial implementation.
- No guaranteed complete parametric case splitting in the first implementation.

## Proposed Public API

Add a new public header:

```text
symengine/polys/groebner.h
```

Initial API:

```cpp
namespace SymEngine {

enum class MonomialOrder {
    Lex,
    GrLex,
    DegRevLex
};

enum class GroebnerAlgorithm {
    Buchberger,
    F5B,
    MoGVW,
    Auto
};

enum class GroebnerStatus {
    Success,
    Cancelled,
    ResourceLimitExceeded,
    NotZeroDimensional,
    UnsupportedCoefficientDomain
};

class GroebnerCancellationToken {
public:
    void cancel();
    bool is_cancelled() const;
};

struct GroebnerOptions {
    MonomialOrder order = MonomialOrder::DegRevLex;
    GroebnerAlgorithm algorithm = GroebnerAlgorithm::Auto;
    bool reduced = true;
    bool interreduce_input = true;
    bool sort_output = true;
    unsigned cancellation_check_interval = 1024;
    const GroebnerCancellationToken * cancellation_token = nullptr;
    unsigned max_s_pairs = 0;
    unsigned max_reduction_steps = 0;
    unsigned max_milliseconds = 0;
};

struct GroebnerStats {
    unsigned input_polys = 0;
    unsigned output_polys = 0;
    unsigned s_pairs_processed = 0;
    unsigned reductions_to_zero = 0;
    unsigned max_basis_size = 0;
};

struct GroebnerResult {
    vec_basic basis;
    vec_sym variables;
    MonomialOrder order;
    GroebnerStatus status;
    GroebnerStats stats;
};

GroebnerResult groebner_basis(const vec_basic &polys,
                              const vec_sym &variables,
                              const GroebnerOptions &options = GroebnerOptions());

GroebnerResult fglm_convert(const GroebnerResult &source,
                            MonomialOrder target_order,
                            const GroebnerOptions &options = GroebnerOptions());

GroebnerResult augmented_groebner_basis(const vec_basic &polys,
                                        const vec_sym &variables,
                                        const RCP<const Symbol> &selected_variable,
                                        const RCP<const Symbol> &auxiliary_variable,
                                        const GroebnerOptions &options = GroebnerOptions());

}
```

Notes:

- `variables` is explicit and ordered. This avoids ambiguity between unknowns and coefficient parameters.
- In `lex`, `variables[0] > variables[1] > ... > variables[n - 1]`. This matches the expected triangular solving convention where the last variable is typically eliminated least aggressively.
- `basis` returns ordinary `Basic` expressions for public API stability. The internal representation remains private.
- `GroebnerStatus::Cancelled` is for non-throwing future variants. The first implementation may throw `SymEngineException("Groebner basis computation cancelled")` if that is more consistent with existing code. The plan should settle this before implementation.
- `max_s_pairs`, `max_reduction_steps`, and `max_milliseconds` are zero for unlimited. They are native equivalents of the external time/memory controls used in `symcse`; they should exit through the same cleanup path as cancellation.
- `augmented_groebner_basis` is inspired by `symcse`'s `zaug - z[i]` pattern. It can initially be implemented as a documented helper around `groebner_basis(polys + {auxiliary_variable - selected_variable}, variables + {auxiliary_variable})`, followed by FGLM when requested.

Later API:

```cpp
GroebnerResult groebner_basis(const std::vector<RCP<const MIntPoly>> &polys,
                              const vec_sym &variables,
                              const GroebnerOptions &options);

RCP<const Set> solve_poly_system(const vec_basic &equations,
                                 const vec_sym &variables,
                                 const GroebnerOptions &options = GroebnerOptions());

RCP<const Set> solve_poly_system_via_univariate_root(const vec_basic &equations,
                                                     const vec_sym &variables,
                                                     const RCP<const Symbol> &selected_variable,
                                                     const RCP<const Symbol> &auxiliary_variable,
                                                     const GroebnerOptions &options = GroebnerOptions());
```

## Internal Representation

Add private implementation files:

```text
symengine/polys/groebner.cpp
symengine/polys/groebner_internal.h
```

Core internal types:

```cpp
using exponent_t = unsigned int;
using Monomial = std::vector<exponent_t>;

template <typename Coeff>
struct Term {
    Monomial monomial;
    Coeff coeff;
};

template <typename Coeff, typename Domain>
class GPoly {
public:
    std::vector<Term<Coeff>> terms; // sorted descending by active monomial order
};
```

Rationale:

- Keep terms sorted by the active monomial order. This makes leading term access O(1), which is central to Buchberger and reduction.
- Do not use `MIntDict` or `MExprDict` directly for the main algorithm because unordered maps require repeated scans to find leading terms.
- Use hash maps temporarily only during polynomial addition/multiplication assembly, then canonicalize into a sorted vector.
- Keep the internal type private so the representation can change later to a heap, packed monomial vector, or F4 matrix format without breaking public API.

Required monomial operations:

- `compare_monomial(a, b, order)`.
- `degree(m)`.
- `divides(a, b)`.
- `quotient(b, a)` for `a | b`.
- `lcm(a, b)`.
- `multiply(m, n)` with overflow checks.

Required polynomial operations:

- canonicalization: remove zero terms, merge equal monomials, sort;
- leading monomial and leading coefficient;
- monomial multiplication;
- scalar multiplication;
- polynomial addition/subtraction;
- S-polynomial construction;
- normal form reduction by a basis;
- interreduction and reduced-basis normalization.

## Coefficient Domains

Groebner algorithms require coefficient division, so the primary abstraction should be a field interface.

Recommended internal concept:

```cpp
template <typename Coeff>
struct CoeffDomain {
    Coeff zero() const;
    Coeff one() const;
    bool is_zero(const Coeff &) const;
    bool is_one(const Coeff &) const;
    Coeff neg(const Coeff &) const;
    Coeff add(const Coeff &, const Coeff &) const;
    Coeff sub(const Coeff &, const Coeff &) const;
    Coeff mul(const Coeff &, const Coeff &) const;
    Coeff div(const Coeff &, const Coeff &) const;
    Coeff normalize(const Coeff &) const;
};
```

### Phase 1: Rational Coefficients

Implement `RationalCoeffDomain` using `rational_class`.

Input handling:

- Accept integer and rational coefficients.
- Convert all integer coefficients to rationals.
- Clear denominators only as an optional optimization, not as the semantic representation.

This gives reliable exact arithmetic and deterministic zero tests.

### Phase 2: Symbolic Coefficients

Support expressions such as:

```text
C1*x*y + C2*x*z**2 - 42*y**3*z
```

when `variables = {x, y, z}` and `C1`, `C2` are coefficient parameters.

Initial practical implementation:

- Add `ExpressionCoeffDomain` using `Expression` or `RCP<const Basic>`.
- Treat all symbols not in `variables` as coefficient expressions.
- Use `add`, `sub`, `mul`, `div`, `neg`, `expand`, and `cancel` where appropriate.
- Use `eq(*expr, *zero)` for zero testing.
- If division by an expression that may be zero is required, proceed formally and optionally record denominators in a diagnostics field later.

Limitations:

- Complete parametric Groebner bases require case splitting on coefficient vanishing conditions. That should not be in the first implementation.
- The first symbolic-coefficient implementation is a generic fraction field over expressions, with conservative simplification.
- A later exact implementation should add multivariate rational functions over parameter symbols, probably backed by a new `MRatPoly`/rational-function layer.

### Phase 3: Finite Fields and Modular Methods

Add `GFpCoeffDomain` using modular integer arithmetic.

Uses:

- faster Groebner computation for numeric systems;
- modular checking;
- future rational reconstruction;
- F4 linear algebra over fields.

## Conversion From SymEngine Expressions

Add converter code in `groebner.cpp`, not in `basic_conversions.h` initially.

Input:

- `vec_basic polys`
- ordered `vec_sym variables`
- coefficient-domain selection inferred from coefficients and options

Process:

1. Expand each input polynomial.
2. Traverse `Add`, `Mul`, `Pow`, `Symbol`, `Integer`, `Rational`.
3. Split each term into:
   - monomial exponents for symbols in `variables`;
   - coefficient expression for everything else.
4. Reject negative or non-integer exponents in variables.
5. Reject functions of variables for polynomial input.
6. Canonicalize zero polynomials away.

Important behavior:

- If a symbol is not in `variables`, it is a coefficient parameter.
- If a user wants `C1` solved as a variable, they must include it in `variables`.

## Buchberger Implementation

Start with Buchberger because it is smaller and easier to verify.

Algorithm outline:

1. Convert and canonicalize inputs.
2. Make all polynomials monic.
3. Interreduce inputs if requested.
4. Initialize unordered pairs `(i, j)`.
5. Select pairs by increasing lcm total degree, then by monomial order. This is a simple sugar-like strategy.
6. Apply Buchberger criteria:
   - product criterion: skip if leading monomials are relatively prime;
   - chain criterion: skip if an existing leading monomial divides the pair lcm and relevant pairs have already been processed.
7. Compute S-polynomial.
8. Reduce S-polynomial by current basis.
9. If remainder is nonzero:
   - make monic;
   - reduce existing basis by the new polynomial if maintaining interreduced basis incrementally;
   - add new pairs.
10. After pair queue exhaustion, compute the reduced basis:
   - reduce each basis polynomial by all others;
   - make monic;
   - remove redundant polynomials whose leading monomial is divisible by another leading monomial;
   - sort output.

Correctness checks to implement in tests:

- Every input polynomial reduces to zero by the output basis.
- Every S-polynomial of output basis elements reduces to zero.
- Output basis is monic.
- In reduced mode, no non-leading term is divisible by another basis leading monomial.

## F4 Roadmap

Do not implement F4 first. Add it after Buchberger is correct.

F4 design:

1. Select a batch of critical pairs by degree.
2. Build symbolic preprocessing set.
3. Construct sparse coefficient matrix.
4. Row-reduce over the active coefficient domain.
5. Convert nonzero rows back to polynomials.
6. Insert new polynomials into the basis.

This will need a reusable exact sparse matrix layer over coefficient domains. The existing matrix classes are `Basic`-oriented and not suitable for high-volume coefficient-domain linear algebra without extra overhead.

## F5B and moGVW Roadmap

`symcse/symcse/polynomials/_groebner_cxx.py` references prototype command-line implementations of F5 and moGVW. Treat these as algorithm inspiration only.

Native SymEngine roadmap:

- Add `GroebnerAlgorithm::F5B` only after the Buchberger path has a stable internal polynomial representation and reduction API.
- Add `GroebnerAlgorithm::MoGVW` after there is an internal signature-basis framework shared with F5-style algorithms.
- Keep these implementations behind `GroebnerAlgorithm::Auto` heuristics only after correctness and deterministic output ordering are well tested.
- Do not call external executables from SymEngine.

## FGLM Implementation

FGLM is required for the target workflow: compute in `degrevlex`, convert to `lex`, then triangular solve.

Initial scope:

- Source basis must be reduced.
- Source ideal must be zero-dimensional.
- Source and destination rings must have the same variables and coefficient domain.
- Destination order initially supports `lex`; later allow any supported order.

Zero-dimensional check:

- For each variable `x_i`, some leading monomial in the source basis must be a pure power of `x_i`.
- If this fails, return or throw `NotZeroDimensional`.

Standard monomial basis:

- Enumerate monomials not divisible by any leading monomial of the source basis.
- Use the pure-power bounds from the zero-dimensional check to keep enumeration finite.
- Sort according to source order for reduction operations and target order for FGLM candidate processing where needed.

Core FGLM loop:

1. Maintain a list of accepted standard monomials in target order.
2. For each candidate monomial:
   - compute its normal form modulo the source basis;
   - express that normal form as a coordinate vector in the source standard monomial basis;
   - test linear dependence against prior vectors.
3. If dependent, produce a new target-order basis polynomial.
4. If independent, add candidate to the target monomial basis and enqueue its products by variables.
5. Stop once the number of accepted basis monomials reaches the vector-space dimension and all necessary border relations are found.

Linear algebra:

- Phase 1 can use dense rational vectors for simplicity.
- Add sparse vectors later.
- Reuse the coefficient-domain abstraction for row operations.

Tests:

- Small zero-dimensional systems with known lex bases.
- Round-trip validation: source basis and FGLM basis mutually reduce generators to zero.
- Non-zero-dimensional input must fail cleanly.
- Cancellation during FGLM must exit cleanly.

## Triangular Solving Plan

The first Groebner PR should expose `fglm_convert`; triangular solving can be a follow-up built on that.

### General Approach

1. Compute or accept a lex Groebner basis.
2. Identify univariate polynomials in the last variable.
3. Solve univariate polynomial using existing `solve_poly`.
4. Substitute each root into the remaining lex basis.
5. Continue one variable at a time.
6. Return:
   - exact finite solutions when all branches close exactly;
   - `ConditionSet` or a structured partial result when solving cannot finish.

Important limitations:

- Symbolic parameters can create branch conditions. A complete implementation needs to track denominators and exceptional cases.
- Multiple roots and positive-dimensional ideals need separate handling.

### First Practical Solver Shape

Implement the `symcse`-inspired shape first:

1. Add an auxiliary variable `zaug`.
2. Add the equation `zaug - selected_variable`.
3. Compute a `degrevlex` Groebner basis in the augmented variables.
4. Convert to `lex` with FGLM.
5. Find exactly one univariate polynomial in `zaug`.
6. Treat the remaining equations as linear in `variables - {selected_variable}` over the field extended by `zaug`.
7. Solve that linear system with existing fraction-free matrix methods where possible.
8. Return solutions expressed in terms of roots of the univariate polynomial, or pass the univariate polynomial to existing exact/numeric univariate root machinery.

This is narrower than a complete triangular decomposition, but it is a useful target because it matches the motivating workflow and the `symcse` chemical-equilibrium tests.

Add a helper analogous to `MPolySys::solve_in_terms_of_upoly`:

```cpp
struct TriangularUnivariateSolution {
    RCP<const Basic> univariate_polynomial;
    RCP<const Symbol> root_variable;
    vec_sym linear_variables;
    DenseMatrix linear_solution;
};

TriangularUnivariateSolution
extract_univariate_linear_shape(const GroebnerResult &lex_basis,
                                const RCP<const Symbol> &root_variable);
```

Validation:

- Verify exactly one usable univariate polynomial exists.
- Verify all non-univariate basis elements are linear in the remaining variables after treating `root_variable` as a coefficient.
- Fall back to `ConditionSet` or an unsupported-shape status when the lex basis does not have this form.

## Cancellation Design

Groebner computation can be extremely expensive. Cancellation must be cooperative and must leave SymEngine in a valid state.

### Do Not Install Global Signal Handlers by Default

SymEngine is a library. Installing a process-wide `SIGINT` handler inside `groebner_basis` would surprise embedding applications and is not thread-safe as a general policy.

### Native Cancellation Token

Implement:

```cpp
class GroebnerCancellationToken {
    std::atomic<bool> cancelled_;
public:
    void cancel();
    bool is_cancelled() const;
};
```

Check the token:

- before and after each S-pair reduction;
- every N reduction steps;
- during polynomial arithmetic loops;
- during FGLM monomial enumeration;
- during linear algebra row reduction.

On cancellation:

- throw a dedicated exception or return `GroebnerStatus::Cancelled`;
- rely on RAII containers for cleanup;
- do not return partial bases as successful results.

### Optional Signal Helper

Later, add an opt-in helper for command-line tools:

```cpp
class ScopedGroebnerInterruptHandler {
public:
    explicit ScopedGroebnerInterruptHandler(GroebnerCancellationToken &token);
    ~ScopedGroebnerInterruptHandler();
};
```

This helper may set the token from `SIGINT`, but it should not be used automatically by library calls.

### Time and Work Limits

Support optional limits in `GroebnerOptions`:

- `max_s_pairs`;
- `max_reduction_steps`;
- `max_milliseconds`.

These should use the same cooperative exit path as manual cancellation.

## File-Level Implementation Plan

### Phase 0: Scaffolding

Files:

- Add `symengine/polys/groebner.h`.
- Add `symengine/polys/groebner_internal.h`.
- Add `symengine/polys/groebner.cpp`.
- Update `symengine/CMakeLists.txt`.
- Add `symengine/tests/polynomial/test_groebner.cpp`.
- Update `symengine/tests/polynomial/CMakeLists.txt`.

Deliverables:

- Public API compiles.
- Monomial order comparators implemented and tested.
- Internal polynomial canonicalization implemented and tested.
- Cancellation token implemented and tested.

Use the "Monomial Ordering Tests" section above as the Phase 0 acceptance list.

### Phase 1: Rational Buchberger

Deliverables:

- Expression-to-internal conversion for integer/rational coefficients.
- Internal-to-expression conversion.
- Normal form reduction.
- Buchberger basis computation over rationals.
- Reduced basis output.

Tests:

- Principal ideal: `{x^2 - 1}`.
- Simple two-variable ideals:
  - `{x*y - 1, y^2 - 1}`;
  - `{x^2 + y^2 - 1, x - y}`;
  - `{x^2 - y, y^2 - x}`.
- `symcse` reference case:
  - `{x^2 - y, x^3 - x}` with lex output containing `x^2 - y`, `x*y - x`, and `y^2 - y`.
- Import the small exact-output cases from "Phase 1 Exact Groebner Basis Tests" above.
- Compare expected bases for `lex`, `grlex`, and `degrevlex` where small enough.
- Property tests for S-polynomial reduction.
- Import the normal-form and `is_groebner`/`is_reduced` cases from "Reduction and Basis-Property Tests".

### Phase 2: Symbolic Coefficients

Deliverables:

- `ExpressionCoeffDomain`.
- Conversion that treats non-variable symbols as coefficients.
- Formal coefficient division.
- Optional coefficient simplification policy.

Tests:

- Inputs containing parameters:
  - `{C1*x + y, x - 1}` over variables `{x, y}`;
  - `{x^2 + K1*x*y, x*y + 2*y^3 - K2}` over variables `{x, y}`;
  - `{C1*x*y + C2*x*z^2 - 42*y^3*z, ...}` once representative additional equations are chosen.
- Import the parameter-coefficient cases from "Symbolic-Parameter and `symcse` Acceptance Tests".
- Verify basis properties by reduction, not only exact expected string output.

### Phase 3: FGLM

Deliverables:

- Zero-dimensional detection.
- Standard monomial enumeration.
- Normal-form vector construction.
- Dense field linear algebra.
- `fglm_convert(degrevlex_basis, Lex)`.

Tests:

- Import all cases from "FGLM and Zero-Dimensional Tests".
- Degrevlex-to-lex conversion.
- Non-zero-dimensional rejection.
- Cancellation during conversion.

### Phase 4: Solver Integration

Deliverables:

- Add nonlinear polynomial-system solving entry point, likely in `solve.h`/`solve.cpp`.
- Keep Groebner basis computation separately usable.
- Implement `extract_univariate_linear_shape` for lex bases.
- Implement augmented-variable solving for a selected variable.
- Implement triangular solve for zero-dimensional lex bases where the shape is one univariate polynomial plus a linear system in the remaining variables.
- Return `ConditionSet` when exact solving cannot complete.

Tests:

- Import the finite exact, inconsistent, positive-dimensional, and high-degree cases from "Polynomial System Solving Tests".
- Systems requiring degrevlex plus FGLM.
- Systems with symbolic parameters where basis can be computed but solving is conditional.
- `symcse`-style augmented systems:
  - add `zaug - z[i]`;
  - compute `degrevlex`;
  - convert to `lex`;
  - extract a univariate polynomial in `zaug`;
  - solve the remaining variables linearly;
  - substitute back and verify residuals are zero or numerically small.

### Phase 5: Performance Improvements

Potential improvements:

- Gebauer-Moller pair criteria.
- Better pair priority queue.
- Sugar degree strategy.
- Packed monomial representation.
- Memory pooling for monomials and terms.
- F4 matrix reductions.
- Modular computation and rational reconstruction.
- Parallel reductions only after cancellation and determinism are well-defined.

## API and Behavior Decisions to Settle Before Coding

1. Cancellation result style:
   - throw `SymEngineException`, or
   - return `GroebnerStatus::Cancelled`.
2. Whether `GroebnerResult::basis` should be `vec_basic` only, or also expose typed polynomial data later.
3. Whether symbolic coefficient division should record denominator assumptions in the result.
4. Whether to add `MRatPoly` now or defer until after a working Groebner implementation.
5. Exact output ordering rules for basis polynomials and terms.
6. Whether memory budgeting should be exposed in the first native API or left to embedding applications. Unlike the `symcse` subprocess path, SymEngine cannot safely use `ulimit` internally.

## Resolved API and Behavior Decisions (as implemented)

These reflect choices made during the first implementation in `symengine/polys/groebner.{h,cpp}`. Revisit them if the rationale changes.

1. Cancellation/limit control flow:
   - The public API never throws for cancellation or resource limits. Internally, sentinel types `GroebnerCancelled`, `GroebnerLimitExceeded`, `GroebnerNotZeroDimensional` are used so status determination does not rely on string-matching exception messages. These are translated to `GroebnerStatus::{Cancelled, ResourceLimitExceeded, NotZeroDimensional}` on the result.
   - Other unexpected `SymEngineException`s (e.g. malformed input) propagate. This keeps user-facing errors distinguishable from cooperative stops.
2. `GroebnerResult::basis` exposes `vec_basic` only in this first PR. Typed polynomial overloads remain a follow-up.
3. Symbolic coefficient division proceeds formally in the field of fractions over the parameter symbols. Denominator/specialization conditions are not yet tracked, but the field is reserved for a future diagnostics extension.
4. `MRatPoly` is deferred. `ExpressionCoeffDomain` provides a workable parameter coefficient layer until performance demands it.
5. Output ordering: when `sort_output` is true, basis polynomials are sorted descending by leading monomial under the active order. Within each polynomial terms are stored descending in the active order. The reduced-basis canonical form (monic, no term divisible by another LM) plus this ordering yields a deterministic output.
6. Cooperative cancellation and `max_s_pairs` / `max_reduction_steps` / `max_milliseconds` are surfaced in `GroebnerOptions`. Memory budgeting is left to the embedder; the library does not call `ulimit` or signal handlers by default.

### Public introspection helpers (added beyond the initial sketch)

To keep property-based tests in-tree and avoid duplicating reduction logic, the following are exposed:

```cpp
RCP<const Basic> normal_form(const RCP<const Basic> &poly,
                             const vec_basic &G,
                             const vec_sym &variables,
                             MonomialOrder order = MonomialOrder::DegRevLex);

bool is_groebner(const vec_basic &G,
                 const vec_sym &variables,
                 MonomialOrder order = MonomialOrder::DegRevLex);

bool is_reduced_basis(const vec_basic &G,
                      const vec_sym &variables,
                      MonomialOrder order = MonomialOrder::DegRevLex);
```

- `normal_form` performs multivariate polynomial reduction modulo `G` under `order`.
- `is_groebner` checks Buchberger's criterion: every S-polynomial of pairs in `G` reduces to zero modulo `G`.
- `is_reduced_basis` additionally verifies that every leading coefficient is one and that no term of any basis element is divisible by another's leading monomial.

### Buchberger statistic accounting

- `stats.s_pairs_processed` counts only pairs whose S-polynomial was actually constructed and reduced.
- Pairs eliminated by the product criterion (leading-monomial coprimality) are not counted as processed and are not counted as reductions to zero — they are a separate criterion-hit class. If a more granular count is needed later, add a `criterion_hits` field.
- `stats.reductions_to_zero` counts only S-polynomials whose normal form was zero after reduction by the current basis.

### Augmented-variable convention

`augmented_groebner_basis(polys, variables, selected_variable, auxiliary_variable, options)` always appends the auxiliary variable last in the augmented variable list. Because lex ordering treats `variables[0] > variables[1] > … > variables[n-1]`, the auxiliary variable becomes the least under lex, and is therefore the variable that appears as the univariate polynomial in the lex Groebner basis. Callers should use the auxiliary variable as `root_variable` when calling `extract_univariate_linear_shape`.

## Recommended First Pull Request

Keep the first PR small enough to review:

1. Add monomial order and internal sparse polynomial infrastructure.
2. Add rational coefficient domain.
3. Add polynomial reduction and Buchberger over rationals.
4. Add public `groebner_basis`.
5. Add property tests and a few exact expected examples.
6. Add cancellation token checks in all long loops.

Do not include FGLM or solver integration in the first PR. Those should build on a verified basis engine.

## Risks

- Symbolic coefficient zero testing is undecidable in general. The implementation must document that formal symbolic coefficients are handled generically and may miss parameter-special cases.
- Groebner bases can grow explosively. Cancellation and work limits should be part of the first implementation, not an afterthought.
- Existing multivariate polynomial classes do not provide the ordered leading-term operations needed for efficient Groebner algorithms.
- FGLM requires zero-dimensional ideals. The solver workflow must detect and report positive-dimensional systems clearly.
- Exact expression simplification during symbolic coefficient arithmetic can dominate runtime. Keep coefficient normalization configurable.
- The augmented-variable workflow depends on variable ordering. The API must document where `auxiliary_variable` is appended and which variable becomes the univariate root variable after lex conversion.

## Acceptance Criteria

The feature should be considered initially complete when:

- `groebner_basis` works natively over exact rational coefficients for small multivariate systems.
- `lex`, `grlex`, and `degrevlex` are tested.
- Reduced basis correctness is verified by S-polynomial reductions.
- Cancellation can stop a long computation without leaking global state or returning a false successful result.
- FGLM converts small zero-dimensional degrevlex bases to lex order.
- The solver layer can use the lex basis for at least simple triangular exact systems.
- The solver layer supports the augmented-variable univariate-plus-linear shape demonstrated by `symcse`.
