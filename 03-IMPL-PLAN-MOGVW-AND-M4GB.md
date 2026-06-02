# Implementation Plan: MoGVW and M4GB for SymEngine Groebner Computation

## Summary

This plan extends the current standalone Groebner implementation with two faster native algorithms:

1. `GroebnerAlgorithm::MoGVW`: a monomial-oriented signature algorithm based on `mogvw.tex`.
2. `GroebnerAlgorithm::M4GB`: a tail-reduced, matrix-backed algorithm inspired by the `/opt/m4gb`.

The current checkout already has:

- `symengine/polys/groebner.h`
- `symengine/polys/groebner_internal.h`
- `symengine/polys/groebner.cpp`
- tests in `symengine/tests/polynomial/test_groebner.cpp`
- algorithms exposed in `GroebnerAlgorithm::{Buchberger, F5B, MoGVW, Auto}`.

The current implementation is a useful correctness baseline: rational and expression coefficient domains, sparse ordered `GPoly`, Buchberger, FGLM, cancellation, limits, basis-property helpers, and first-shape solving. MoGVW and M4GB should reuse the public API and verification helpers but should introduce their own internal engines rather than deforming Buchberger into unrelated code.

## Goals

- Add faster Groebner computation for hard finite-field and rational systems.
- Preserve deterministic reduced-basis output for a fixed input/order/domain.
- Keep cancellation and work limits as first-class behavior in every algorithm.
- Keep Singular out of the runtime and build graph.
- Make benchmark coverage broad enough to compare:
  - current Buchberger;
  - future F4;
  - MoGVW;
  - M4GB;
  - degrevlex plus FGLM versus direct lex.

## Non-Goals

- No bridge to Singular, Sage, FGb, OpenF4, or the `m4gb` executable.
- No external dependency on the `m4gb` build system or generated solver executables. Reuse source directly where helpful, under the separate author license agreement.
- No complete comprehensive Groebner systems with parameter case splitting in these algorithm PRs.
- No noncommutative Groebner bases.
- No attempt to make MoGVW handle every coefficient domain in its first PR. Start where the theory and data structures are cleanest.

## Current Implementation Constraints

### Existing `GPoly`

`GPoly<Coeff, Domain>` stores sorted sparse terms:

- leading term access is O(1);
- reductions erase from the front during normal form;
- coefficient arithmetic is dispatched through `RationalCoeffDomain` and `ExpressionCoeffDomain`;
- monomials are `std::vector<unsigned int>`.

This is correct but not ideal for fast algorithms:

- vector monomials are allocation-heavy;
- normal form repeatedly canonicalizes and sorts;
- pair queues are simple vectors resorted after each insert;
- there is no finite-field coefficient domain;
- no packed monomial id, basis-leading-monomial index, or reducer lookup table exists yet.

### Required Foundation Before MoGVW/M4GB

Add a shared internal layer before implementing either algorithm:

- `symengine/polys/groebner_monomial_index.h`
- `symengine/polys/groebner_matrix.h`
- `symengine/polys/groebner_signature.h`
- `symengine/polys/groebner_domains.h`

These can remain private headers.

The shared layer should provide:

- packed monomial ids for fixed variable count/order;
- monomial multiplication/division/lcm without heap allocation in hot loops;
- a leading-monomial divisor index;
- finite-field coefficient domains;
- sparse-to-dense row conversion for selected monomial supports;
- cancellation hooks in matrix and reduction loops;
- stats counters shared by all algorithms.

## Public API Changes

Extend `GroebnerAlgorithm`:

```cpp
enum class GroebnerAlgorithm {
    Buchberger,
    F4,
    F5B,
    MoGVW,
    M4GB,
    Auto
};
```

Extend `GroebnerOptions`:

```cpp
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

    // New algorithm controls.
    unsigned matrix_batch_size = 0;       // 0 = algorithm default
    unsigned max_matrix_rows = 0;         // 0 = unlimited
    unsigned max_matrix_columns = 0;      // 0 = unlimited
    unsigned max_degree = 0;              // 0 = unlimited
    bool verify_with_buchberger = false;  // expensive debug option
    bool allow_signature_algorithms = true;
};
```

Extend `GroebnerStats`:

```cpp
struct GroebnerStats {
    unsigned input_polys = 0;
    unsigned output_polys = 0;
    unsigned s_pairs_processed = 0;
    unsigned reductions_to_zero = 0;
    unsigned max_basis_size = 0;

    unsigned pairs_rejected_product = 0;
    unsigned pairs_rejected_chain = 0;
    unsigned pairs_rejected_signature = 0;
    unsigned pairs_rejected_rewritten = 0;
    unsigned pairs_rejected_syzygy = 0;

    unsigned matrices_built = 0;
    unsigned max_matrix_rows = 0;
    unsigned max_matrix_columns = 0;
    unsigned rows_reduced = 0;
    unsigned rows_zero = 0;

    unsigned monomials_indexed = 0;
    unsigned reducer_lookup_hits = 0;
    unsigned reducer_lookup_misses = 0;
};
```

Status handling remains unchanged: cancellation and limits return `GroebnerStatus::{Cancelled, ResourceLimitExceeded}` rather than throwing past the public API.

## Coefficient-Domain Roadmap

### Required for M4GB

M4GB is most valuable for finite-field, dense or overdefined systems. Add finite-field domains first:

- `GFpCoeffDomain` for prime fields.
- Later `GF2nCoeffDomain` if Boolean/MQ benchmarks justify it.

Proposed public entry point can be delayed. Internal support can initially be selected by options or by future overloads:

```cpp
GroebnerResult groebner_basis_mod(const vec_basic &polys,
                                  const vec_sym &variables,
                                  unsigned modulus,
                                  const GroebnerOptions &options);
```

Do not emulate finite fields through rational coefficients with modular-looking integers. The matrix algorithms need cheap packed field operations.

### Required for MoGVW

Start with `QQ` and `GF(p)`.

Avoid `ExpressionCoeffDomain` for the first MoGVW PR. Signature criteria require exact divisibility and comparisons; symbolic expression division can silently introduce branch-sensitive denominators. Add expression coefficients only after the rational/finite-field behavior is tested.

## Shared Monomial Infrastructure

### Packed Monomial Ids

Implement a packed representation alongside the current `std::vector<unsigned int>` representation:

```cpp
struct PackedMonomial {
    uint64_t id;
    unsigned total_degree;
};
```

For early support, use a fixed-width exponent pack when possible:

- up to 8 variables with 8-bit exponents;
- up to 16 variables with 4-bit exponents;
- fallback to vector monomials for large exponents or variable counts.

For M4GB-style dense monomial enumeration under degrevlex, also add an ordinal/index codec:

- map monomials to increasing degrevlex indices;
- enumerate all monomials up to a degree;
- decode an index back to exponent vector.

The `m4gb/lib/monomial_degrevlex.hpp` reference uses degree-layer indexing and lookup tables. With the separate author license agreement, SymEngine can adapt this implementation directly, but it should still be reshaped to fit SymEngine naming, allocation, cancellation, and test conventions.

### Divisor Index

Both MoGVW and M4GB need fast tests of the form:

- find whether any leading monomial divides `m`;
- find a preferred reducer for `m`;
- find whether a monomial belongs to an existing labeled-monomial set.

Add:

```cpp
class LeadingMonomialIndex {
public:
    void insert(MonomialId lm, size_t basis_index);
    void erase(MonomialId lm, size_t basis_index);
    bool find_divisor(MonomialId m, size_t *basis_index) const;
    std::vector<size_t> find_all_divisors(MonomialId m) const;
};
```

Implementation stages:

1. Simple scan, behind the interface.
2. Degree bucketed scan.
3. Divisibility trie or hashed quotient buckets.

MoGVW specifically benefits from membership tests in large sets; M4GB benefits from finding the first divisor while constructing matrix multiples.

## Matrix Infrastructure

Add a private row-reduction layer:

```cpp
template <typename Coeff, typename Domain>
struct MatrixRow {
    std::vector<unsigned> columns;
    std::vector<Coeff> values;
    Signature signature; // optional; used by signature algorithms
};

template <typename Coeff, typename Domain>
class GroebnerMatrixReducer {
public:
    MatrixReductionResult reduce(std::vector<MatrixRow<Coeff, Domain>> &rows,
                                 const std::vector<MonomialId> &columns,
                                 const GroebnerOptions &options);
};
```

Stages:

1. Sparse rows over `QQ` and `GF(p)`.
2. Dense rows over `GF(p)` for M4GB-style overdefined systems.
3. Optional block/simd rows for `GF(2)` or `GF(2^n)`.

One-sided elimination is required by MoGVW: rows with higher signatures may be reduced by lower-signature rows, but not the reverse.

Cancellation checks:

- before each pivot row;
- every N row operations;
- before growing row/column vectors;
- after symbolic preprocessing completes a batch.

## MoGVW Plan

### Algorithm Summary

`mogvw.tex` presents a monomial-oriented GVW algorithm. Instead of organizing work around J-pairs of labeled polynomials, it organizes work around labeled monomials:

```text
labeled monomial = (m, (u, f))
```

where:

- `m` is an available leading monomial;
- `(u, f)` is a labeled polynomial/vector relation;
- `lm(f)` divides `m`;
- signature is `(m / lm(f)) * lm(u)`.

The key design difference:

- GVW/F5 asks: for a given signature, what is the smallest polynomial leading monomial?
- MoGVW asks: for a given monomial, what is the smallest signature that can generate it?

This allows:

- avoiding explicit J-pair generation in the main loop;
- hash-table membership checks instead of expensive divisor searches in several criteria;
- degree-wise lifting of labeled monomials;
- matrix-style reductions grouped by degree/signature.

### MoGVW Data Structures

Add private types:

```cpp
struct Signature {
    MonomialId monomial;
    unsigned generator_index;
};

struct LabeledPolynomial {
    Signature signature;
    GPolyLike polynomial;
    unsigned id;
};

struct LabeledMonomial {
    MonomialId monomial;       // m
    unsigned labeled_poly_id;  // points to (u, f)
    Signature signature;       // (m / lm(f)) * lm(u)
    bool primitive;
    bool syzygy;
    bool lifted;
};
```

Indexes:

- `unordered_map<MonomialId, LabeledMonomialId> best_by_monomial`
- `unordered_map<Signature, LabeledMonomialId> best_by_signature`
- `unordered_set<Signature> principal_syzygy_signatures`
- degree buckets for unlifted labeled monomials
- `LeadingMonomialIndex` for primitive labeled polynomials

Ordering:

- polynomial monomial order is `options.order`, initially `DegRevLex`;
- signature order is position-over-term:
  - `x^alpha e_i < x^beta e_j` if `i > j`;
  - if `i == j`, compare monomials by the polynomial order.

### MoGVW Insertion Rule

When inserting a labeled monomial `(m, (u, f))`:

1. If it is syzygy, add it to the syzygy index.
2. If `m` is not present, insert it.
3. If `m` is present, keep the labeled monomial with smaller signature.
4. If replacement happens, enqueue the displaced labeled monomial for mutual reduction or mark it rejected depending on criteria.

This is central to the "smallest signature for a monomial" invariant.

### MoGVW Criteria

Implement in this order:

1. Rewritten criterion:
   - reject a labeled monomial when another labeled monomial with a smaller/equivalent signature already covers the same generated work.
   - In matrix mode, perform this during sorting of rows by signature, because equal signatures become adjacent.
2. Principal syzygy criterion:
   - for signature `x^alpha e_i`, reject if an existing labeled monomial implies a principal syzygy with a smaller generator position.
   - Start with principal syzygies only; full syzygy storage can be added later.
3. LCM criterion:
   - reject lifted multiples that are not the primitive J-pair and whose primitive J-pair is already covered.
   - Track this with a hash set keyed by `(primitive_id, other_primitive_id, lcm_monomial)` or an equivalent covered-J-pair key.

Stats must count each rejection type separately.

### MoGVW Scalar First, Matrix Second

Phase MGVW-1 should implement a non-matrix scalar version:

1. Convert input to internal labeled polynomials.
2. Initialize primitive labeled monomials.
3. Initialize principal syzygy index.
4. Process unlifted labeled monomials degree by degree.
5. For each variable, lift `x_i * labeled_monomial`.
6. Apply criteria.
7. Mutual-reduce colliding labeled monomials.
8. Insert primitive outputs into the basis.
9. At termination, extract primitive polynomials and reduce them to canonical basis form.

This establishes correctness with minimal linear algebra.

Phase MGVW-2 should implement the matrix version from the paper:

1. `lift(todo, G)` produces reduction candidates.
2. `append(H, G)` adds reducers needed for monomials in the candidate tails.
3. `eliminate(H)` performs one-sided row reduction sorted by signature.
4. `update(P, G)` inserts new primitive labeled monomials.

The matrix version should share the row-reduction layer with future F4/M4GB.

### MoGVW Integration

Dispatch from `groebner_basis`:

```cpp
switch (options.algorithm) {
case GroebnerAlgorithm::MoGVW:
    return groebner_basis_mogvw(...);
case GroebnerAlgorithm::Auto:
    if (can_use_mogvw(polys, variables, options)) ...
}
```

Initial `can_use_mogvw`:

- order must be `DegRevLex` or `GrLex`;
- coefficient domain must be `QQ` or `GF(p)`;
- no symbolic parameters;
- input count should be at least two;
- no denominator/negative exponent conversion issues.

Return reduced basis using existing interreduce/reduce helpers until MoGVW has native final reduction.

### MoGVW Testing

Correctness tests:

- every small exact basis from `01-IMPL-PLAN-GROEBNER.md` must match Buchberger by ideal equality;
- `{x^2 - y, x^3 - x}` must produce same reduced basis as Buchberger;
- toric/binomial ideal must pass `is_groebner` and input-reduction checks;
- Katsura-3 must match the Buchberger basis after final reduction;
- zero-dimensional systems must support FGLM after MoGVW result.

Signature-specific tests:

- collision replacement keeps smaller signature;
- principal syzygy criterion rejects known redundant lifted monomials;
- rewritten criterion rejects duplicate-signature rows;
- cancellation during lift, append, eliminate, and update returns `Cancelled`;
- `max_degree`, `max_matrix_rows`, and `max_matrix_columns` return `ResourceLimitExceeded`.

Benchmarks:

- compare `Buchberger` versus `MoGVW` on degrevlex cases from `benchmarks/groebner_bench.cpp`;
- add boolean/MQ benchmarks after finite fields exist.

## M4GB Plan

### Algorithm Summary

M4GB is designed around tail-reduced polynomials:

- every basis polynomial is monic;
- all tail terms are non-reducible by current basis leading monomials;
- multiples of basis polynomials are stored as dense coefficient vectors over non-reducible monomials;
- reducible monomials are represented by matrix entries keyed by their leading monomial;
- multiplication by a monomial performs full reduction immediately through the matrix/polymatrix.

Compared with F4:

- F4 processes batches of critical pairs into large matrices.
- M4GB still processes critical-pair work but maintains a persistent database of reduced basis multiples.
- This can avoid repeated construction/reduction of the same multiples and reduce the degree-regularity "staircase" behavior.

M4GB is especially attractive for dense overdefined finite-field systems, not for symbolic-parameter systems.

### M4GB Data Structures

Implement private classes:

```cpp
template <typename Coeff, typename Domain>
class TailReducedPolynomial {
public:
    MonomialId leading_monomial;
    std::vector<Coeff> tail; // dense over non-reducible monomial index
};

class NonReducibleMonomialIndex {
public:
    MonomialId upper_bound;
    std::vector<MonomialId> dense_index;
    unordered_map<MonomialId, unsigned> dense_inverse;
    void increase_upper_bound(MonomialId m);
    void decrease_upper_bound(MonomialId m);
};

template <typename Coeff, typename Domain>
class PolyMatrix {
public:
    struct Entry {
        MonomialId basis_lm;
        std::vector<Coeff> tail;
        int generation;
        bool computed;
    };
};
```

The `./m4gb` reference names the same concepts:

- `upper_bound`
- `dense_index`
- `dense_invindex`
- `matrix`: leading monomial -> reduced tail
- `basis`: leading monomial -> degree
- `basisitem`: per-basis metadata and known multiples
- lazy computation/generation tracking

Reimplement these ideas using SymEngine-owned code.

### M4GB Core Invariants

Maintain:

1. `basis` contains only monic tail-reduced polynomials.
2. `dense_index` contains monomials not divisible by any current basis leading monomial and below `upper_bound`.
3. `matrix[m]` exists for reducible monomials below `upper_bound`.
4. `matrix[m].tail` represents the fully reduced tail of the polynomial with leading monomial `m`.
5. A polynomial/vector has no trailing zero coefficient.
6. `upper_bound` changes only through controlled increase/decrease routines that keep all indexes consistent.

Every debug build should have an invariant checker:

```cpp
void assert_m4gb_invariants(const M4GBState &state);
```

Use it in tests behind a cheap flag.

### M4GB Operations

#### Increase Upper Bound

When `upper_bound` grows to include monomial `m`:

1. Enumerate monomials from current `upper_bound` through `m`.
2. For each monomial:
   - if divisible by a basis leading monomial, create a matrix entry for the preferred divisor;
   - otherwise append it to `dense_index`.
3. Update inverse indexes.
4. Check cancellation periodically.

Start with simple divisor scan. Add the M4GB-style basis sieve later.

#### Insert Basis Polynomial

When inserting a new basis polynomial with leading monomial `lm`:

1. Store its reduced tail in matrix form.
2. Add `lm` to the leading-monomial index.
3. Remove `lm` from `dense_index` if present.
4. For every existing non-reducible monomial now divisible by `lm`, move it from `dense_index` into `matrix`.
5. Generate lazy matrix entries for relevant multiples.
6. Shrink or reindex as needed.

#### Get Reduced Multiple

`get_u_g(u, g)` computes the fully reduced multiple `u * g`:

1. For each nonzero tail term `c*m` of `g`, compute `u*m`.
2. If `u*m` is non-reducible, put it in the dense tail.
3. If `u*m` is reducible, subtract the matrix entry for `u*m`.
4. Cleanup trailing zeros.

This is the core M4GB optimization.

#### Shrink

When new basis elements make dense monomials reducible:

1. Identify dense columns that are now reducible.
2. Use their matrix entries to eliminate them from existing rows.
3. Remove those columns from `dense_index`.
4. Update inverse indexes.

Implement eager shrink first. Add lazy generation/shrink later only if benchmarks show it matters.

#### Selection and Update

Maintain critical pairs ordered by lcm degree and monomial:

```cpp
struct CriticalPair {
    MonomialId lcm;
    MonomialId lm_a;
    MonomialId lm_b;
};
```

Process a bounded batch:

1. Select pairs of the current lcm degree.
2. Increase upper bound to cover selected lcms.
3. Build reduced rows from the two required multiples.
4. Row-reduce.
5. Convert nonzero rows to tail-reduced candidate polynomials.
6. Insert candidates into the basis.
7. Remove redundant basis elements whose leading monomials are divisible by new leading monomials.

This can start closer to F4-style batches and evolve toward the one-pair-at-a-time M4GB style once the polymatrix is reliable.

### M4GB Domain Scope

Phase M4GB-1:

- `GF(p)` only;
- `DegRevLex` only;
- fixed variable count with packed monomial ids;
- no symbolic parameters;
- no FGLM in the same PR, but result must be usable by existing FGLM when zero-dimensional and convertible to `QQ`-style basis is not required.

Phase M4GB-2:

- `QQ` via rational coefficients if row operations are efficient enough;
- modular-rational reconstruction can be a later separate project;
- optional `GF(2)`/`GF(2^n)` dense bit-matrix support.

M4GB over `ExpressionCoeffDomain` is not a near-term target.

### M4GB Integration

Dispatch:

```cpp
case GroebnerAlgorithm::M4GB:
    if (!can_use_m4gb(...)) {
        return unsupported_domain_or_fallback(...);
    }
    return groebner_basis_m4gb(...);
```

Do not silently fall back from explicitly requested `M4GB` to Buchberger unless the project decides that is acceptable. Prefer:

- `UnsupportedCoefficientDomain` if the domain/order is unsupported;
- fallback only for `Auto`.

`Auto` heuristic after M4GB is stable:

- use M4GB for finite-field, degrevlex, many quadratic equations, no parameters;
- use MoGVW for rational/finite-field degrevlex systems with meaningful signature filtering;
- use Buchberger for tiny systems and symbolic parameters.

### M4GB Testing

Unit tests:

- monomial id encoding/decoding under degrevlex;
- dense index insertion/removal;
- divisor detection;
- `get_u_g` reduces multiples to non-reducible tails;
- shrink removes reducible dense columns without changing represented polynomials;
- row reduction over `GF(p)`;
- basis insertion preserves invariants.

Correctness tests:

- all small rational tests after `QQ` support;
- finite-field case from SymPy:
  - `GF(7)[x, y, z]`;
  - `{3*x^2 + y*z - 5*x - 1, 2*x + 3*x*y + y^2, x - 3*y + x*z - 2*z^2}`;
  - compare reduced basis under canonical modular normalization.
- `m4gb/testdata/gf31/31_16_mq-test0.in` as a skipped/slow integration fixture until parser and field support exist.
- `m4gb/testdata/gf256/256_16_mq-test0.in` only after `GF(2^n)` exists.

Property tests:

- every input reduces to zero by M4GB output;
- `is_groebner` returns true;
- M4GB output and Buchberger output mutually reduce to zero on small cases;
- cancellation during upper-bound growth, matrix entry computation, row reduction, and update exits cleanly.

Benchmarks:

- dense random MQ over `GF(31)` with `m = n + 1` and `m = 2n`;
- n from 6 upward, stopping at a time limit;
- compare output stats: matrix rows/columns, basis size, critical pairs, memory proxy.

## Benchmark Plan

`benchmarks/groebner_bench.cpp` is the first native benchmark harness. It should remain lightweight and dependency-free.

Current benchmark categories:

- small SymPy elimination case in lex and degrevlex;
- SymPy minpoly case in direct lex and degrevlex plus FGLM;
- SymPy FGLM symmetric case in direct lex and degrevlex plus FGLM;
- Katsura-3 in lex and degrevlex;
- toric/binomial degrevlex case;
- algebraic-composition case inspired by `demo_groebner_bench_sympy.py`.

Next benchmark additions after finite fields:

- random dense quadratic systems over `GF(31)`;
- root-forced MQ systems over `GF(31)`;
- Boolean/HFE systems if `GF(2)`/`GF(2^n)` support is added;
- direct comparison of `Buchberger`, `F4`, `MoGVW`, and `M4GB` on the same inputs.

Benchmark output should always include:

- algorithm;
- order;
- domain;
- total time;
- basis size;
- S-pairs processed;
- rejected pairs by criterion;
- matrix count;
- max rows/columns;
- rows reduced to zero;
- cancellation/limit status.

## Implementation Phases

### Phase 0: Benchmark and Instrumentation

Deliverables:

- `benchmarks/groebner_bench.cpp`.
- Extended `GroebnerStats`.
- `GroebnerAlgorithm::F4` and `GroebnerAlgorithm::M4GB` enum entries.
- Explicit `UnsupportedCoefficientDomain` behavior for unsupported requested algorithms.

Tests:

- existing Groebner tests still pass;
- benchmark executable builds.

### Phase 1: Finite Fields

Deliverables:

- `GFpCoeffDomain`.
- modular expression conversion for integer/rational inputs.
- modular normal form and Buchberger path.
- finite-field test fixtures from SymPy.

Tests:

- arithmetic unit tests for `GF(p)`;
- compare small modular Groebner bases to expected output;
- verify canonical representative normalization.

### Phase 2: Shared Packed Monomial and Matrix Layer

Deliverables:

- monomial id codec for degrevlex;
- leading-monomial divisor index;
- sparse matrix row representation;
- row reduction over `QQ` and `GF(p)`;
- one-sided row reduction support.

Tests:

- encode/decode round trips;
- monomial ordering consistency with current `compare_monomial_less`;
- divisor lookup against brute force;
- row-reduction invariants.

### Phase 3: Scalar MoGVW

Deliverables:

- signature type and order;
- labeled polynomial/monomial storage;
- lift loop;
- collision handling;
- principal syzygy and rewritten criteria;
- final extraction and reduction.

Tests:

- small bases equal Buchberger by ideal equality;
- criterion unit tests;
- cancellation and limits.

### Phase 4: Matrix MoGVW

Deliverables:

- `lift`, `append`, `eliminate`, `update` functions.
- one-sided matrix elimination sorted by signature.
- matrix stats.

Tests:

- scalar and matrix MoGVW produce equivalent bases;
- matrix rows with identical signatures are rewritten correctly;
- benchmark against Buchberger.

### Phase 5: M4GB Polymatrix Core

Deliverables:

- tail-reduced polynomial representation.
- non-reducible monomial dense index.
- polymatrix entries.
- upper-bound increase/decrease.
- basis insertion and shrink.
- reduced multiple computation.

Tests:

- invariants after every operation;
- small hand-constructed polymatrix reductions.

### Phase 6: M4GB Solver Loop

Deliverables:

- critical pair selection;
- matrix row construction from reduced multiples;
- update loop;
- final reduced-basis conversion to public `vec_basic`.

Tests:

- small finite-field bases match Buchberger;
- dense MQ benchmarks complete;
- cancellation and limits.

### Phase 7: Auto Heuristics

Deliverables:

- `GroebnerAlgorithm::Auto` chooses between Buchberger, MoGVW, M4GB, and later F4.
- Add decision stats/logging in debug builds.

Initial heuristic:

- symbolic parameters -> Buchberger;
- requested lex -> Buchberger or degrevlex plus FGLM only when caller asks through solver/FGLM path;
- `GF(p)`, degrevlex, many dense quadratics -> M4GB;
- `QQ`/`GF(p)`, degrevlex, no parameters, medium system -> MoGVW;
- tiny systems -> Buchberger.

## Risk Register

- Signature algorithms are correctness-sensitive. Mitigation: keep scalar MoGVW first, add property checks, and verify against Buchberger for small cases.
- M4GB depends on a complex persistent matrix invariant. Mitigation: implement the polymatrix as an independently tested component before adding the solver loop.
- Direct source reuse from `./m4gb` still needs engineering adaptation. Mitigation: import only cohesive components, preserve provenance notes required by the separate agreement, and wrap them in SymEngine-native APIs/tests rather than carrying over the standalone solver architecture unchanged.
- Finite-field support is missing today. Mitigation: implement `GF(p)` before M4GB.
- Monomial packing can overflow. Mitigation: detect unsupported exponents and fall back to vector monomials or reject the optimized algorithm.
- Matrix algorithms can consume memory quickly. Mitigation: expose `max_matrix_rows`, `max_matrix_columns`, `max_degree`, and cooperative cancellation from the first PR.
- Reduced output may differ in ordering or normalization. Mitigation: final interreduce with existing canonical reducer and verify `is_reduced_basis`.

## Acceptance Criteria

MoGVW is acceptable when:

- explicitly requested `GroebnerAlgorithm::MoGVW` computes correct reduced bases for all small rational test cases without falling back;
- zero-dimensional MoGVW output can be passed to `fglm_convert`;
- signature criteria are unit-tested;
- cancellation and limits work in all major loops;
- benchmark output shows at least parity with Buchberger on selected degrevlex systems.

M4GB is acceptable when:

- `GF(p)` support exists;
- explicitly requested `GroebnerAlgorithm::M4GB` computes correct reduced bases for small finite-field systems;
- polymatrix invariants are tested independently;
- dense MQ benchmarks run through the native benchmark harness;
- cancellation and matrix limits work;
- reused M4GB source is covered by the separate author license agreement and is integrated as SymEngine-native code rather than as an external executable.
