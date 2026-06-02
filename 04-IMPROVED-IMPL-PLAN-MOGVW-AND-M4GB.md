# Implementation Plan: F5B, MoGVW and M4GB for SymEngine Groebner Computation

## Summary

This plan extends the current Groebner implementation in `symengine/polys/groebner.{h,cpp}` with three faster native algorithms:

1. `GroebnerAlgorithm::F5B`: the F5 algorithm in Buchberger's style — Sun & Wang (JSC), see `JSC-F5B.tex`. A signature-based Buchberger with two redundancy criteria (Syzygy + Rewritten) and F5-reduction. Reference implementation: `sympy/polys/groebnertools.py::_f5b`.
2. `GroebnerAlgorithm::MoGVW`: a monomial-oriented signature algorithm faithful to Algorithm 1 / Procedure 2 of `mogvw.tex`. Extends F5B's signature infrastructure to labeled monomials with a "smallest signature for each monomial" invariant.
3. `GroebnerAlgorithm::M4GB`: a tail-reduced, matrix-backed algorithm adapted from `/opt/m4gb` (Makarim & Stevens, ISSAC'17 — see `/opt/m4gb/ISSAC17-MS-M4GB.pdf`).

F5B is staged first because it is the simplest signature algorithm to verify (SymPy ships an extensively-tested reference), and the `Signature` / position-over-term ordering / F5-reduction / rewritten+syzygy criteria infrastructure it builds is reused verbatim by MoGVW.

The current checkout has:
- a Buchberger / Rational + Expression baseline in `groebner.cpp`,
- `GroebnerAlgorithm::{Buchberger, F5B, MoGVW, Auto}` already in the enum (`groebner.h:20-25`), but `groebner_basis` ignores `options.algorithm` and unconditionally dispatches to `buchberger<...>` at `groebner.cpp:388-425`,
- cancellation via two internal exception types (`GroebnerCancelled`, `GroebnerLimitExceeded`) translated to `GroebnerStatus::{Cancelled, ResourceLimitExceeded}` at `groebner.cpp:426-430`,
- a private `GPoly<Coeff,Domain>` with sparse sorted terms and `RationalCoeffDomain` / `ExpressionCoeffDomain` (`groebner_internal.h`).

MoGVW and M4GB must reuse the **public** API and the existing cancellation/result types but introduce their own private internal engines — they should not try to deform `buchberger` into a shared routine.

**Hard constraints carried from prior plans:**
- No bridge to Singular / Sage / FGb / OpenF4 / the `m4gb` executable; reuse `/opt/m4gb` source ideas only (the standalone author-license file in that tree covers excerpting algorithm logic; do not copy headers verbatim — re-author for SymEngine style).
- `/opt/m4gb` relies on **compile-time** template parameters (`MAXVARS`, `FIELDSIZE`, `MAX_INT_DEGREE`, …). SymEngine sees variable count and degree only at runtime. Every array / hash table / field characteristic must be runtime-sized (`std::vector`, runtime modulus), while keeping inner loops allocation-free.

---

## Public API Changes

### Enum

Add `M4GB`:

```cpp
enum class GroebnerAlgorithm {
    Buchberger,
    F5B,
    MoGVW,
    M4GB,
    Auto,
};
```

`F5B`, `MoGVW`, and `M4GB` all become real algorithms in this plan. Until each phase lands, an explicit request for that algorithm returns `GroebnerStatus::UnsupportedCoefficientDomain`. Document the staging order in the public header comment.

### `GroebnerOptions`

Extend with knobs the new algorithms need. All new fields default to `0` / `false`, preserving current behavior:

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

    // New: finite-field coefficient domain. 0 -> rational/Expression as today.
    // Any prime p with 2 <= p < 2^31 is valid. Composite/large p returns
    // UnsupportedCoefficientDomain.
    uint64_t modulus = 0;

    // New: matrix / degree controls.
    unsigned matrix_batch_size = 0;     // 0 = algorithm default
    unsigned max_matrix_rows = 0;       // 0 = unlimited
    unsigned max_matrix_columns = 0;    // 0 = unlimited
    unsigned max_degree = 0;            // 0 = unlimited

    // New: verification toggle for tests / debugging.
    bool verify_with_buchberger = false;
};
```

### `GroebnerStats`

Add rejection / matrix counters used by the new algorithms (keep existing fields):

```cpp
struct GroebnerStats {
    // existing
    unsigned input_polys = 0;
    unsigned output_polys = 0;
    unsigned s_pairs_processed = 0;
    unsigned reductions_to_zero = 0;
    unsigned max_basis_size = 0;

    // new (F5B + MoGVW — signature-based shared counters)
    unsigned rejected_by_syzygy = 0;     // Syzygy / "comparable" criterion
    unsigned rejected_by_rewritten = 0;  // Rewritten criterion
    unsigned f5_reductions = 0;          // F5-reduction steps (F5B only)
    unsigned labeled_monomials_lifted = 0; // MoGVW only
    unsigned rejected_by_lcm = 0;        // MoGVW LCM criterion
    unsigned collisions_resolved = 0;    // MoGVW only

    // new (M4GB / future F4)
    unsigned matrices_built = 0;
    unsigned max_matrix_rows_seen = 0;
    unsigned max_matrix_columns_seen = 0;
    unsigned rows_reduced_to_zero = 0;
    unsigned polymatrix_entries_built = 0;
    unsigned polymatrix_generation = 0;
};
```

### Dispatch in `groebner_basis`

Replace the unconditional `buchberger` call (`groebner.cpp:388-425`) with:

```cpp
switch (resolve_algorithm(polys, variables, options)) {
case GroebnerAlgorithm::Buchberger: return groebner_basis_buchberger(...);
case GroebnerAlgorithm::F5B:        return groebner_basis_f5b(...);
case GroebnerAlgorithm::MoGVW:      return groebner_basis_mogvw(...);
case GroebnerAlgorithm::M4GB:       return groebner_basis_m4gb(...);
}
```

`resolve_algorithm` materializes `Auto` to a concrete choice (heuristic table in §"Auto Heuristic" below).

For an explicitly requested algorithm whose preconditions fail (e.g. `M4GB` with `modulus == 0`), return `GroebnerStatus::UnsupportedCoefficientDomain` rather than silently falling back. `Auto` is permitted to fall back without surfacing an error.

---

## Shared Infrastructure (Phase 1 prerequisites)

### 1. Finite-field coefficient domain `GFpCoeffDomain`

Both algorithms need cheap modular arithmetic.

```cpp
class GFpCoeffDomain {
    uint64_t p_;        // prime modulus, 2 <= p < 2^31
public:
    explicit GFpCoeffDomain(uint64_t p) : p_(p) {
        // Constructor caller guarantees primality + range; tested with
        // GroebnerOptions::modulus validation at the API boundary.
    }
    uint64_t modulus() const { return p_; }
    uint64_t zero() const { return 0; }
    uint64_t one()  const { return 1; }
    bool is_zero(uint64_t c) const { return c == 0; }
    bool is_one (uint64_t c) const { return c == 1; }

    // Invariant on inputs to add/sub/mul/neg: 0 <= a, b < p_.
    // With p_ < 2^31, both a+b and a*b fit in uint64_t with room to spare.
    uint64_t add(uint64_t a, uint64_t b) const {
        uint64_t s = a + b;
        return s >= p_ ? s - p_ : s;
    }
    uint64_t sub(uint64_t a, uint64_t b) const {
        return a >= b ? a - b : a + p_ - b;
    }
    uint64_t mul(uint64_t a, uint64_t b) const {
        return (a * b) % p_;
    }
    uint64_t neg(uint64_t a) const { return a == 0 ? 0 : p_ - a; }

    // Extended Euclidean inverse; precondition gcd(a, p_) == 1.
    uint64_t inv(uint64_t a) const;
    uint64_t div(uint64_t a, uint64_t b) const { return mul(a, inv(b)); }
};
```

Notes for implementers:
- Do **not** use `std::pow` or floating-point for modular inverse. Use either extended Euclidean or Fermat (`pow_mod(a, p-2, p)`).
- The `p < 2^31` cap bounds *coefficient* magnitude only. It says nothing about exponents — those are bounded separately by `PackedMonomial`'s bit-schedule and `RuntimeDegRevLexCodec`'s `max_degree`. With `p < 2^31` we get `a + b < 2^32` (safe in `uint64_t`) and `a * b < 2^62` (safe in `uint64_t` with two bits of headroom). Lifting the cap means `__int128_t` or Barrett/Montgomery; for one-prime modular Groebner you almost never want anything larger than the Mersenne prime `2^31 - 1`. Document this in the public header next to `GroebnerOptions::modulus`.
- `normalize(c)` is a no-op for GF(p) but should exist for parity with `RationalCoeffDomain` so generic code compiles.

Conversion from `RCP<const Basic>` to `uint64_t` mod `p`: integer terms are reduced; rational terms must have denominators coprime to `p` (otherwise `UnsupportedCoefficientDomain`).

### 2. Monomial representations

Two distinct representations are needed; the plan must **not** conflate them:

**(a) `PackedMonomial`** — used by MoGVW and as the "natural" key for the leading-monomial index.

```cpp
struct PackedMonomial {
    uint64_t bits;       // bit-packed exponents
    uint32_t total_degree;
    uint32_t num_vars;   // for the comparator + decoding

    bool operator==(const PackedMonomial&) const;
};
struct PackedMonomialHash {
    size_t operator()(const PackedMonomial&) const;
};
```

Choose a packing schedule at construction based on `num_vars` and observed max exponent: 8 vars × 8 bits, 16 vars × 4 bits, 32 vars × 2 bits, or fall back to a `std::vector<exponent_t>`-backed variant for systems outside those ranges (degree overflow / variable overflow returns `UnsupportedCoefficientDomain` for explicit M4GB; `Auto` falls back to Buchberger).

DegRevLex comparison must be implemented explicitly — bit-packed comparison alone does **not** give degrevlex. Implement as: first compare `total_degree`; if equal, compare reverse-lex on the unpacked exponents. Provide an `operator<` for the chosen `MonomialOrder` only after validating against the existing `compare_monomial_less` in `groebner_internal.h:64` via property tests.

**(b) `OrdinalMonomialId`** — a `uint64_t` that *is* the degrevlex index produced by the combinatorial codec. This is what M4GB uses for `dense_index`, `matrix`, `basis`, and `upper_bound`.

```cpp
using OrdinalMonomialId = uint64_t;

class RuntimeDegRevLexCodec {
public:
    RuntimeDegRevLexCodec(unsigned num_vars, unsigned max_degree);

    OrdinalMonomialId to_index(const std::vector<exponent_t>& exps) const;
    void from_index(OrdinalMonomialId id, std::vector<exponent_t>& out) const;
    unsigned degree(OrdinalMonomialId id) const;

    // Monomial arithmetic on ordinals (uses unpack/repack internally):
    OrdinalMonomialId multiply(OrdinalMonomialId a, OrdinalMonomialId b) const;
    bool divides(OrdinalMonomialId a, OrdinalMonomialId b) const;
    OrdinalMonomialId divide (OrdinalMonomialId a, OrdinalMonomialId b) const;
    OrdinalMonomialId lcm    (OrdinalMonomialId a, OrdinalMonomialId b) const;

    OrdinalMonomialId next(OrdinalMonomialId id) const { return id + 1; }
    OrdinalMonomialId max_id() const;   // L[max_degree + 1] - 1
};
```

This is the runtime port of `/opt/m4gb/lib/monomial_degrevlex.hpp::monomial_degrevlex_intcodec::static_data_t`. Replace the compile-time `T2`, `T3`, `L` arrays with `std::vector<std::vector<std::vector<uint64_t>>>` allocated in the constructor and populated from `detail::multiset_coefficient` (port of `/opt/m4gb/lib/detail.hpp:44-50`). Check for `uint64_t` overflow during construction; if `multiset_coefficient(num_vars + 1, max_degree)` overflows, return `UnsupportedCoefficientDomain` for explicit M4GB.

The plan's two monomial representations are bridged by an unpack/repack at the codec boundary. MoGVW data structures stay in `PackedMonomial`; M4GB structures stay in `OrdinalMonomialId`. **Never** type a structure as `PackedMonomial` and then key it by an ordinal index, or vice versa.

### 3. Sparse matrix row & reducer

```cpp
template <typename Coeff>
struct MatrixRow {
    std::vector<uint64_t> columns;  // monomial ids (ordinal *or* packed; one or the other per matrix instance)
    std::vector<Coeff>    values;   // parallel to columns
    Signature             signature; // ignored when the matrix is not signature-driven
};

template <typename Coeff, typename Domain>
class GroebnerMatrixReducer {
public:
    struct Result {
        unsigned rows_in;
        unsigned rows_out;
        unsigned rows_zero;
    };
    Result reduce(std::vector<MatrixRow<Coeff>>& rows,
                  const GroebnerOptions& options,
                  GroebnerStats& stats);
};
```

Implementation contract:
- Cancellation check every `options.cancellation_check_interval` pivot steps; throw `GroebnerCancelled` on `is_cancelled()`.
- Throw `GroebnerLimitExceeded` when `rows.size() > options.max_matrix_rows` (when non-zero) before growth, and similarly for columns.
- For MoGVW: one-sided elimination (a row may only be reduced by rows with **strictly smaller** signature).
- For M4GB / F4: classical Gaussian elimination, monic on output.

---

## Implementing F5B (Phase 2)

F5B is "F5 in Buchberger's style" (Sun & Wang, `JSC-F5B.tex`). It is the same overall loop as Buchberger — pop a critical pair, compute its S-polynomial, reduce, insert if nonzero — with two additions:

1. Every polynomial carries a **signature** (a labeled-polynomial structure) tracking its derivation from the input.
2. Critical pairs are rejected up-front by two criteria (**Syzygy** + **Rewritten**) based on signatures, and reduction is replaced by **F5-reduction**, which only allows reductions whose multiplier raises the signature.

These additions cut redundant zero-reductions dramatically on systems where Buchberger spends most of its time reducing S-polys to zero (the Katsura / cyclic / Czichowski family, among others).

**Reference implementation to mirror:** `sympy/polys/groebnertools.py::_f5b` lines 567-688 plus its helpers `sig`, `lbp`, `sig_cmp`, `sig_key`, `sig_mult`, `lbp_sub`, `lbp_mul_term`, `lbp_cmp`, `lbp_key`, `critical_pair`, `cp_cmp`, `cp_key`, `s_poly`, `is_rewritable_or_comparable`, `f5_reduce`, `red_groebner`. Tested against SymPy's `test_groebnertools.py` fixtures (port them — see Tasks below).

### Data structures

Reuse the same `Signature` and signature-ordering machinery defined for MoGVW (§"Implementing MoGVW" below) — F5B *is* the first consumer of that infrastructure. The MoGVW section's `SignatureLess` lives in a shared header `groebner_signature.h` and is used by both algorithms.

```cpp
// signature, polynomial, and the serial 'k' assigned when added to the basis.
// (Equivalent to SymPy's lbp(signature, polynomial, number).)
template <typename Coeff, typename Domain>
struct LabeledPoly {
    Signature   signature;
    GPoly<Coeff, Domain> polynomial;
    unsigned    number;     // monotonically increasing, used by the rewritten criterion
};

// Critical pair, packed as in SymPy (groebnertools.py:402-434): the larger-
// signature side comes first; S-polynomial of cp is um*f - vm*g.
template <typename Coeff, typename Domain>
struct CriticalPair {
    Signature                 sig_um_f;     // signature of um * f
    PackedMonomial            um;
    const LabeledPoly<Coeff, Domain>* f;
    Signature                 sig_vm_g;     // signature of vm * g
    PackedMonomial            vm;
    const LabeledPoly<Coeff, Domain>* g;
};
```

### Signature ordering — match SymPy's `sig_cmp`

```text
u < v  iff  u.index > v.index
         OR (u.index == v.index AND u.monomial < v.monomial)
```

i.e. **the higher generator index gives the lower signature**. The MoGVW plan adopts the same convention; both algorithms link against the same `SignatureLess`. Critical-pair ordering is then `cp_key` (a lex pair of the two side signatures); the main loop pops the *smallest* CP.

This convention makes the principal-syzygy of `(f_i, f_j)` with `i < j` always reachable before its dependents and matches the JSC paper.

### The two criteria — `is_rewritable_or_comparable`

Following `sympy/polys/groebnertools.py:489-513`, both criteria are checked in one routine:

```cpp
template <typename Coeff, typename Domain>
bool is_rewritable_or_comparable(const Signature& s,
                                 unsigned num,
                                 const std::vector<LabeledPoly<Coeff, Domain>>& B,
                                 GroebnerStats& stats)
{
    for (const auto& h : B) {
        // Comparable (Syzygy criterion): h has a *higher* index and its
        // leading monomial divides s.monomial.
        if (h.signature.generator_index > s.generator_index &&
            packed_divides(h.polynomial.leading_monomial(), s.monomial)) {
            ++stats.rejected_by_syzygy;
            return true;
        }
        // Rewritable: same index, h has a *later* number, and h.signature.monomial
        // divides s.monomial.
        if (h.signature.generator_index == s.generator_index &&
            h.number > num &&
            packed_divides(h.signature.monomial, s.monomial)) {
            ++stats.rejected_by_rewritten;
            return true;
        }
    }
    return false;
}
```

Note the convention flip: SymPy stores signatures with `index = i + 1` and "smaller index = larger signature", so its `sign[1] < Sign(h)[1]` test compares indices the same way ours does. Match the SymPy logic line-for-line; do not "fix" the direction without re-running the SymPy fixtures.

### F5-reduction

JSC paper §4.1 ("F5-Reduction") and `sympy/polys/groebnertools.py:516-564`. A labeled polynomial `f` is F5-reduced by `B`: while there exists `h ∈ B` with non-zero `polynomial`, such that `LM(h) | LM(f)` and `sig(t·h) < sig(f)` where `t = LM(f)/LM(h)`, replace `f := f - (LT(f)/LT(h))·h`. Stop when no such `h` exists.

```cpp
template <typename Coeff, typename Domain>
LabeledPoly<Coeff, Domain> f5_reduce(LabeledPoly<Coeff, Domain> f,
                                     const std::vector<LabeledPoly<Coeff, Domain>>& B,
                                     const GroebnerOptions& options,
                                     GroebnerStats& stats)
{
    if (f.polynomial.is_zero()) return f;
    for (;;) {
        check_cancellation(options);
        bool reduced = false;
        for (const auto& h : B) {
            if (h.polynomial.is_zero()) continue;
            if (!packed_divides(h.polynomial.leading_monomial(), f.polynomial.leading_monomial())) continue;
            PackedMonomial t = packed_divide(f.polynomial.leading_monomial(), h.polynomial.leading_monomial());
            Signature th_sig = signature_mul(h.signature, t);
            if (!SignatureLess{}(th_sig, f.signature)) continue;
            // Note: SymPy explicitly does NOT call is_rewritable_or_comparable here
            // ("is in general slower than without"); follow that choice.
            f = lbp_sub(f, lbp_mul_term(h, t));
            ++stats.f5_reductions;
            reduced = true;
            break;
        }
        if (!reduced || f.polynomial.is_zero()) return f;
    }
}
```

### Main loop

Mirrors `sympy/polys/groebnertools.py:600-688` step-for-step:

```cpp
template <typename Coeff, typename Domain>
std::vector<GPoly<Coeff, Domain>>
f5b_impl(std::vector<GPoly<Coeff, Domain>> F,
         const GroebnerOptions& options,
         GroebnerStats& stats)
{
    // 1. Pre-reduce inputs (Pernici / Becker-Weispfenning p.203).
    //    Loop: B := each F[i] reduced by F[0..i-1]; repeat until fixpoint.
    F = pernici_pre_reduce(F, options);

    // 2. Initialize labeled polynomials with sig (1, i+1) and number i+1.
    std::vector<LabeledPoly<Coeff, Domain>> B;
    B.reserve(F.size());
    for (size_t i = 0; i < F.size(); ++i)
        B.push_back({ make_signature(one_monomial(), i + 1u), F[i], (unsigned)(i + 1) });

    // 3. Sort B by leading monomial, descending (largest LM first).
    std::sort(B.begin(), B.end(), [](const auto& a, const auto& b) {
        return PackedMonomialLess{}(b.polynomial.leading_monomial(),
                                    a.polynomial.leading_monomial());
    });

    // 4. Build initial critical-pair set, sorted by cp_key, descending.
    auto CP = build_initial_critical_pairs(B);
    sort_critical_pairs_desc(CP);

    unsigned k = (unsigned)B.size();

    while (!CP.empty()) {
        check_cancellation(options);
        // Pop the *smallest* critical pair.
        CriticalPair<Coeff, Domain> cp = std::move(CP.back());
        CP.pop_back();
        ++stats.s_pairs_processed;

        // Reject by either criterion on either side.
        if (is_rewritable_or_comparable(cp.sig_um_f, cp.f->number, B, stats)) continue;
        if (is_rewritable_or_comparable(cp.sig_vm_g, cp.g->number, B, stats)) continue;

        // Build the S-poly via lbp_sub(lbp_mul_term(f, cp.um), lbp_mul_term(g, cp.vm)).
        LabeledPoly<Coeff, Domain> s = s_polynomial(cp);

        // F5-reduce, then make monic.
        LabeledPoly<Coeff, Domain> p = f5_reduce(std::move(s), B, options, stats);
        p.polynomial.make_monic();
        p.number = ++k;

        if (p.polynomial.is_zero()) {
            ++stats.reductions_to_zero;
            continue;
        }

        // Prune old CPs that become redundant against the new p.
        CP.erase(std::remove_if(CP.begin(), CP.end(), [&](const auto& cpx) {
            return is_rewritable_or_comparable(cpx.sig_um_f, cpx.f->number, {p}, stats) ||
                   is_rewritable_or_comparable(cpx.sig_vm_g, cpx.g->number, {p}, stats);
        }), CP.end());

        // Add new critical pairs (p, g) for each g in B, filtered by redundancy
        // against {p} (same filter SymPy uses).
        for (const auto& g : B) {
            if (g.polynomial.is_zero()) continue;
            CriticalPair<Coeff, Domain> cp_new = make_critical_pair(p, g);
            if (is_rewritable_or_comparable(cp_new.sig_um_f, cp_new.f->number, {p}, stats)) continue;
            if (is_rewritable_or_comparable(cp_new.sig_vm_g, cp_new.g->number, {p}, stats)) continue;
            CP.push_back(std::move(cp_new));
        }
        sort_critical_pairs_desc(CP);

        // Insert p into B in LM-descending position (sympy uses linear scan).
        insert_in_lm_order(B, std::move(p));

        if (options.max_basis_size && B.size() > options.max_basis_size)
            throw GroebnerLimitExceeded{};
    }

    // 5. Reduce the basis (red_groebner).
    return red_groebner(extract_polynomials(B), options);
}
```

The `red_groebner` step (sympy:691-720) selects the minimal-LM subset and inter-reduces it. Reuse the existing `interreduce` / `reduce` helpers in `groebner.cpp` rather than reimplementing — that's the same canonical-form code path that Buchberger uses today, so F5B output is bit-identical to Buchberger output (modulo input order).

### Domain support

F5B's signature machinery is **purely monomial** — no coefficient operations enter the criteria or the signature comparisons. The only places where the coefficient domain matters are the polynomial arithmetic in `f5_reduce`, `s_polynomial`, and `make_monic`. Therefore F5B works over any `Coeff/Domain` pair that `GPoly` already supports:

- **`RationalCoeffDomain`** — the primary target; produces bit-identical bases to Buchberger.
- **`GFpCoeffDomain`** — once Phase 1 ships, free of charge.
- **`ExpressionCoeffDomain`** — works "out of the box" modulo the genericity caveat from Phase 6 (the rewritten criterion's `Num(h) > num` check is exact regardless of coefficient ring). This makes F5B the *first* algorithm in the plan that gives signature-criterion pruning to symbolic-coefficient inputs — a substantial improvement for the parametric workloads exercised by the symcse corpus.

`can_use_f5b(polys, variables, options)` therefore has weaker preconditions than MoGVW or M4GB: any `MonomialOrder`, any supported coefficient domain, no minimum input size. Explicit `algorithm == F5B` always proceeds (no `UnsupportedCoefficientDomain` fallback) unless the input cannot be parsed into `GPoly`.

### F5B precondition: `ExpressionCoeffDomain` zero-test

The F5-reduction loop and `make_monic` both call `is_zero`. Over `ExpressionCoeffDomain`, that predicate uses post-`expand` structural equality (`groebner_internal.h:126`), which is correct for *generically* nonzero coefficients but flags `K1**2 - K2*K3` as nonzero even on the subvariety `K1**2 = K2*K3`. Document this clearly in the F5B header comment. Symcse-style inputs work fine with this assumption — see the test plan below.

### Tasks

1. Implement `Signature`, `SignatureLess`, signature arithmetic (`signature_mul`, etc.) in `symengine/polys/groebner_signature.h`. This is the **shared** signature header consumed by both F5B (Phase 2) and MoGVW (Phase 3); design it for both consumers up front.
2. Implement `LabeledPoly`, `CriticalPair`, `lbp_sub`, `lbp_mul_term`, `make_critical_pair`, `s_polynomial`, `cp_key`-style ordering.
3. Implement `is_rewritable_or_comparable` and `f5_reduce` exactly mirroring SymPy's logic.
4. Implement `pernici_pre_reduce` (the iterative input pre-reduction at the top of `_f5b`).
5. Implement the main loop and wire `GroebnerAlgorithm::F5B` to `groebner_basis_f5b`.
6. **Port the SymPy F5B test corpus** from `/opt-3/cpython-v3.13-apt-deb/lib/python3.13/site-packages/sympy/polys/tests/test_groebnertools.py` into `symengine/tests/polynomial/test_groebner_f5b.cpp`:
   - `test_groebner_f5b` (small QQ basis, lex+grlex).
   - `test_benchmark_minpoly_f5b` (3-variable QQ minimum-polynomial system in lex).
   - `test_benchmark_katsura3_f5b` (Katsura-3 over `ZZ`/QQ in lex and grlex).
   - `test_benchmark_kastura_4_f5b` (Katsura-4 — tag as slow benchmark).
   - `test_benchmark_czichowski_f5b` (Czichowski — tag as slow benchmark).
   - `test_benchmark_cyclic_4_f5b` (cyclic-4 — tag as slow benchmark).
   - Hand-translate inputs and expected outputs; do not pull in a Python runtime.
7. Run all existing Buchberger tests in `symengine/tests/polynomial/test_groebner.cpp` *also* through F5B (via a `TEST_F5B` macro or parametrized fixture) and assert ideal equality with the Buchberger output.
8. Run the symcse parametric corpus (defined in Phase 6) *also* through F5B once `ExpressionCoeffDomain` support is wired in. Expectation: F5B over parametric inputs produces the same basis as Buchberger but with fewer reductions to zero. Record `f5_reductions`, `rejected_by_syzygy`, `rejected_by_rewritten`, `reductions_to_zero` in the benchmark output for direct comparison.

### Acceptance criteria

- All ported SymPy F5B fixtures pass.
- Every Buchberger-passing test in `test_groebner.cpp` also passes for F5B (with byte-identical reduced bases under `sort_output`).
- `rejected_by_syzygy + rejected_by_rewritten > 0` on Katsura-3, Katsura-4, and cyclic-4 (otherwise F5B is silently degenerating into Buchberger).
- Cancellation triggers `Cancelled` from inside `f5_reduce` and from inside the main loop.
- F5B over `ExpressionCoeffDomain` returns the same basis as Buchberger on every symcse test (with the genericity caveat documented).

---

## Implementing MoGVW (Phase 3)

### Concept mapping (`mogvw.tex` §3)

A *labeled monomial* is `(m, (u, f))`:
- `m` is the *available* monomial being tracked.
- `(u, f)` is a labeled polynomial — `u` is the signature carrier in a free module, `f` is its image in `R`.
- The induced signature is `S = (m / lm(f)) * lm(u)`.

The algorithm tracks, for each monomial `m`, the labeled polynomial `(u, f)` with the **smallest** signature whose `lm(f)` divides `m`.

### Data structures

```cpp
struct Signature {
    PackedMonomial monomial;
    unsigned generator_index;   // 1 .. input_polys
};

// Position-over-term ordering. The convention chosen here matches mogvw.tex
// Lemma 3.1 (prec_s compatible with prec_p): for the same position, compare
// monomials by the polynomial order; across positions, the smaller index has
// the *larger* signature (so principal syzygies of (f_i, f_j) with i < j are
// available when their J-pair is examined).
//
// Tests in Phase 3 must verify this convention against worked examples in
// mogvw.tex §3.4 (Example 3.4: r_3, r_4, r_5 ordering). The same SignatureLess
// is exercised by F5B in Phase 2; both must agree.
struct SignatureLess {
    bool operator()(const Signature& a, const Signature& b) const {
        if (a.generator_index != b.generator_index)
            return a.generator_index > b.generator_index;
        return PackedMonomialLess{}(a.monomial, b.monomial);
    }
};

struct LabeledPolynomial {
    Signature signature;
    GPoly<uint64_t, GFpCoeffDomain> polynomial;     // or rational variant
    unsigned id;
};

struct LabeledMonomial {
    PackedMonomial m;                    // the "available" monomial
    unsigned       labeled_poly_id;      // index into ctx.labeled_polys
    Signature      signature;            // (m / lm(f)) * lm(u)
    bool           is_syzygy;            // polynomial component reduces to 0
    bool           lifted;               // x_i * this has been queued for each i
};
```

`MoGVWContext` holds:

```cpp
struct MoGVWContext {
    std::vector<LabeledPolynomial> labeled_polys;
    std::vector<LabeledMonomial>   G;                  // the full set
    std::unordered_map<PackedMonomial, size_t, PackedMonomialHash>
        best_by_monomial;                              // m -> index into G
    std::set<Signature, SignatureLess>
        principal_syzygies;                            // sig where polynomial = 0
    LeadingMonomialIndex lm_index;                     // for divisor / cover tests
    GroebnerOptions options;
    GroebnerStats   stats;
};
```

### `mutualreduce` — faithful to `mogvw.tex` Procedure 2

The original (Procedure 2 at `mogvw.tex:672-724`) is:

```text
mutualreduce(brm, G):
    if brm is REDUCIBLE by G and (LCM-rejected or Syzygy-rejected or Rewritten-rejected):
        return
    reduce brm to brm''  by G                          // top-reduce; brm'' may equal brm
    if m'' != 0 and brm'' collides with brm' = (m'', (v, g)) in G:
        if sig(brm'') strictly precedes sig(brm'):
            G := (G \ {brm'}) ∪ {brm''}
            mutualreduce(brm', G)                       // recurse on the *displaced* one
        // else: no-op, brm'' is silently dropped
    else:
        G := G ∪ {brm''}
```

Implementation, with the bug fixed from the prior draft:

```cpp
void mutual_reduce(LabeledMonomial brm, MoGVWContext& ctx) {
    check_cancellation(ctx);

    // 1. Combined gate: only reject when the labeled monomial is REDUCIBLE.
    bool reducible = ctx.lm_index.has_divisor(brm.m);
    if (reducible) {
        if (is_rejected_by_lcm(brm, ctx))       { ++ctx.stats.rejected_by_lcm;       return; }
        if (is_rejected_by_syzygy(brm, ctx))    { ++ctx.stats.rejected_by_syzygy;    return; }
        if (is_rejected_by_rewritten(brm, ctx)) { ++ctx.stats.rejected_by_rewritten; return; }
    }

    // 2. Top-reduce in place. May leave brm unchanged when not reducible.
    LabeledMonomial brm_pp = top_reduce(brm, ctx);

    // 3. Collision resolution.
    if (!is_zero_monomial(brm_pp.m)) {
        auto it = ctx.best_by_monomial.find(brm_pp.m);
        if (it != ctx.best_by_monomial.end()) {
            LabeledMonomial& existing = ctx.G[it->second];
            if (SignatureLess{}(brm_pp.signature, existing.signature)) {
                // Strictly smaller signature: displace and recurse on the
                // *old* one (it might still produce something useful).
                LabeledMonomial displaced = existing;     // copy before mutation
                existing = brm_pp;                         // replace in G in place
                ++ctx.stats.collisions_resolved;
                mutual_reduce(displaced, ctx);             // recurse with displaced
            }
            // Otherwise: drop brm_pp silently (no insert, no recurse).
            return;
        }
    }

    // 4. No collision: insert.
    size_t new_id = ctx.G.size();
    ctx.G.push_back(brm_pp);
    if (is_zero_monomial(brm_pp.m)) {
        // Syzygy labeled monomial. Track its signature for the criterion.
        ctx.principal_syzygies.insert(brm_pp.signature);
    } else {
        ctx.best_by_monomial[brm_pp.m] = new_id;
        ctx.lm_index.insert(brm_pp.m, new_id);
    }
}
```

The bug in the prior draft: it fell through to step 4 (unconditional insert into `G`) even when a collision resolution happened, and only wrote the new value into `best_by_monomial`, leaving `G` inconsistent. The version above keeps `G` as the canonical set, with `best_by_monomial` as an index into it.

### Main loop (faithful to `mogvw.tex` Algorithm 1, lines 1-11)

```cpp
void groebner_basis_mogvw_impl(MoGVWContext& ctx) {
    initialize_G_from_input(ctx);  // primitive labeled monomials + initial principal syzygies
    unsigned liftdeg = max_primitive_degree(ctx);

    for (;;) {
        bool progressed = false;
        for (size_t i = 0; i < ctx.G.size(); ++i) {
            LabeledMonomial& brm = ctx.G[i];
            if (brm.is_syzygy || brm.lifted) continue;
            if (degree(brm.m) > liftdeg) continue;

            for (unsigned v = 0; v < ctx.options.num_variables; ++v) {
                LabeledMonomial lifted = multiply_by_variable(brm, v);
                mutual_reduce(lifted, ctx);
                ++ctx.stats.labeled_monomials_lifted;
            }
            brm.lifted = true;
            progressed = true;
            liftdeg = std::max(liftdeg, max_primitive_degree(ctx));
        }

        unsigned maxcpdeg = max_cp_lcm_degree(ctx);
        if (maxcpdeg > liftdeg + 1) {
            liftdeg = maxcpdeg - 1;     // goto step 4 in the paper
            continue;
        }
        if (!progressed) break;
    }

    extract_reduced_basis(ctx);
}
```

Note the "goto step 4" jump in the paper. Implement it as the `liftdeg` reset above; do not literal-`goto`.

### MoGVW criteria

Implement strictly in the order the paper uses them inside `mutualreduce`:

1. **LCM criterion** (Corollary 3.3 in `mogvw.tex`): reject `brm = t(lm(f), (u,f))` when `t·lm(f) ≠ lcm(lm(f), lm(g))` for some `(lm(g), (v,g)) ∈ G` and the J-pair `(lcm/lm(f))·(u,f)` is already covered.
2. **Principal syzygy criterion**: reject when `brm` has a signature that is divisible by an already-recorded principal-syzygy signature with smaller generator index. Start with **only principal syzygies**; full syzygy tracking is out of scope for the first PR.
3. **Rewritten criterion** (Corollary 3.2): reject when there exists `(m', (u', f')) ∈ G` with the same generator index whose signature divides `brm.signature` and whose monomial does *not* equal `brm.m`.

Each criterion has its own stat counter; tests must verify the counters separately on hand-constructed inputs.

### MoGVW dispatch & preconditions

`can_use_mogvw(polys, variables, options)`:

- `options.order ∈ {DegRevLex, GrLex}`,
- `options.modulus == 0` (rational) **or** prime and < 2^31 (GF(p)),
- no symbolic parameters (`!has_symbolic_parameters(polys, variables)`),
- input count ≥ 2,
- packed-monomial schedule for `variables.size()` is feasible.

For `algorithm == MoGVW` and any check failing, return `UnsupportedCoefficientDomain`. For `Auto`, fall back silently to Buchberger.

After MoGVW computes its internal `G`, **always** post-process through the existing `interreduce` / `reduce` helpers from `groebner.cpp` to get a canonical reduced basis. This both protects against MoGVW-specific normalization drift and lets the existing tests pass without modification.

### MoGVW scope

This plan ships **scalar MoGVW only** (no matrix elimination of grouped signatures). Matrix-MoGVW is deferred — its design depends on a sound `GroebnerMatrixReducer` and on the M4GB row infrastructure, both of which mature in Phase 4.

---

## Implementing M4GB (Phase 4 + 5)

M4GB's defining trick is a persistent, fully tail-reduced cache of basis multiples — the *polymatrix*. Every term operation (multiply by monomial, reduce) is answered by table lookup into the polymatrix.

### Data structures

Note: all monomial IDs in M4GB are **ordinal IDs** produced by `RuntimeDegRevLexCodec`. Ordering and successor (`+1`) operate on these IDs directly.

```cpp
using DensePoly = std::vector<uint64_t>;   // GFp coefficients; index i corresponds to dense_index[i]

struct PolyMatrixEntry {
    DensePoly tail;            // reduced tail (no trailing zeros)
    int generation;            // matches the polymatrix's current generation if up-to-date; -1 = not computed
    OrdinalMonomialId basis_lm; // the basis polynomial this is a multiple of
};

struct BasisRep {
    std::unordered_set<OrdinalMonomialId> matrix_multiples;
    OrdinalMonomialId next_mul; // BASISSIEVE: next multiplier to try
    OrdinalMonomialId next_res; // BASISSIEVE: lm * next_mul
};

class PolyMatrix {
public:
    GFpCoeffDomain field;
    RuntimeDegRevLexCodec codec;

    OrdinalMonomialId upper_bound = 0;        // strictly larger than every monomial in the matrix

    std::vector<OrdinalMonomialId>            dense_index;     // sorted ascending
    std::unordered_map<OrdinalMonomialId, size_t> dense_invindex;

    // Reduced multiples of basis polys: lm (reducible) -> entry
    std::unordered_map<OrdinalMonomialId, PolyMatrixEntry> matrix;

    // Basis: lm -> degree(lm)
    std::map<OrdinalMonomialId, unsigned> basis;
    std::unordered_map<OrdinalMonomialId, BasisRep> basis_info;

    int generation = 0;
};
```

`OrdinalMonomialId` is `uint64_t`, so `std::unordered_map` works without a custom hash.

### Invariants (assert in debug builds)

After every public operation:

1. Every `basis[lm]` polynomial is monic and tail-reduced by the rest of the basis.
2. Every `id ∈ dense_index` is **not** divisible by any current basis leading monomial and satisfies `id < upper_bound`.
3. Every `id ∈ matrix` is divisible by exactly one basis leading monomial and satisfies `id < upper_bound`.
4. No `DensePoly` has a trailing zero.
5. `dense_invindex[dense_index[i]] == i` for all `i`.

Provide a `void assert_invariants(const PolyMatrix&)` callable from tests; gate routine calls behind a cheap toggle so it isn't compiled out by accident.

### Core operations

#### `increase_upper_bound(m)`

Enumerate ordinals from the current `upper_bound` through `m` in ascending degrevlex order. For each new ordinal:
- if some basis lm divides it, register a matrix entry pointing at the chosen divisor (eager or lazy per the lazy flag);
- otherwise, append to `dense_index` and `dense_invindex`.

Phase 4 ships the simple "scan basis polynomials looking for a divisor" version. The BASISSIEVE variant from `/opt/m4gb/src/m4gb.hpp:355-394` (each basis poly carries `next_mul` / `next_res` cursors) is an optional follow-up if benchmarks demand it.

#### `decrease_upper_bound(m)`

Remove all `dense_index` entries `> m`, prune `matrix` keys `> m`, and (with BASISSIEVE) reset each basis cursor to `lm * minimum_of_degree(...)`. Mirrors `m4gb.hpp:325-349`.

#### `insert_basis_polynomial(lm, tail)`

1. Insert `(lm, PolyMatrixEntry{tail, generation, lm})` into `matrix`.
2. Add `lm → degree(lm)` to `basis`; initialize `basis_info[lm]`.
3. Remove `lm` from `dense_index` if present.
4. For every `id ∈ dense_index` with `lm | id`: create a matrix entry for `id` (lazy or eager) and remove `id` from `dense_index`.

Steps 3-4 may shift indices in `dense_index` — recompute `dense_invindex` accordingly (incremental update if the changed prefix is short, otherwise full rebuild — mirror `m4gb.hpp:868-897`).

#### `get_u_g(u, g_tail)` — the heart of M4GB

```cpp
DensePoly get_u_g(OrdinalMonomialId u,
                  const DensePoly& g_tail,
                  PolyMatrix& pm)
{
    DensePoly ret;
    if (g_tail.empty()) return ret;

    // Pre-size based on the largest possible target dense index (= dense
    // position of u * top monomial of g_tail). Avoids reallocation but does
    // not allocate the full dense_index.size().
    OrdinalMonomialId top = pm.codec.multiply(u, pm.dense_index[g_tail.size() - 1]);
    auto ub_it = std::lower_bound(pm.dense_index.begin(), pm.dense_index.end(), top);
    ret.reserve(ub_it - pm.dense_index.begin() + 1);

    for (size_t i = g_tail.size(); i-- > 0; ) {
        if (pm.field.is_zero(g_tail[i])) continue;
        OrdinalMonomialId m = pm.codec.multiply(u, pm.dense_index[i]);

        auto dit = pm.dense_invindex.find(m);
        if (dit != pm.dense_invindex.end()) {
            // Non-reducible target: accumulate in dense slot.
            if (dit->second >= ret.size()) ret.resize(dit->second + 1, 0);
            ret[dit->second] = pm.field.add(ret[dit->second], g_tail[i]);
        } else {
            // Reducible target: subtract c * reducer_tail.
            auto mit = pm.matrix.find(m);
            if (mit == pm.matrix.end())
                throw std::runtime_error("M4GB invariant violated: m neither dense nor reducible");
            const DensePoly& reducer = ensure_entry_current(mit, pm);
            for (size_t j = 0; j < reducer.size(); ++j) {
                if (pm.field.is_zero(reducer[j])) continue;
                uint64_t prod = pm.field.mul(g_tail[i], reducer[j]);
                if (j >= ret.size()) ret.resize(j + 1, 0);
                ret[j] = pm.field.sub(ret[j], prod);
            }
        }

        // Cooperative cancellation inside the hot loop.
        if ((i & (pm.options.cancellation_check_interval - 1)) == 0)
            check_cancellation(pm);
    }
    cleanup_trailing_zeros(ret, pm.field);
    return ret;
}
```

Differences from the prior draft: incremental `resize` (not `vector(dense_index.size(), 0)`), descending iteration matching `m4gb.hpp:703-781`, cancellation check inside the loop.

#### `ensure_entry_current(it, pm)` — lazy generation handling

Phase 4 ships **eager** computation: every matrix entry is computed at insertion time and `generation` is never used (always equal to `pm.generation`). Phase 5 adds the lazy / shrink-on-demand machinery from `m4gb.hpp:497-562` if benchmarks justify it; it is not on the critical path for correctness.

#### `shrink()`

After a batch of new basis elements lands, some columns in `dense_index` may have become reducible. Walk `dense_index`, move any entry that now has a corresponding `matrix` key out, eliminate that column from every active row (or bump generation under lazy shrink), and rebuild `dense_invindex`. Mirrors `m4gb.hpp:822-898`.

For Phase 4 ship the eager flavor: synchronously reduce every active row by the eliminated column. Lazy shrink is a Phase-5 optional.

### M4GB solver loop

Sketch (mirrors `m4gb.hpp:1378-1425`):

```cpp
void groebner_basis_m4gb_impl(M4GBContext& ctx) {
    initialize_from_input(ctx);            // converts input to DensePoly, runs pre-rowreduce
    while (!ctx.critical_pairs.empty() || !ctx.to_add.empty()) {
        update(ctx);                       // drain to_add into basis, update CP
        ctx.matrix.shrink();
        if (ctx.critical_pairs.empty() && ctx.to_add.empty()) break;
        selection(ctx);                    // build matrix rows, row-reduce, push results to to_add
    }
    extract_basis(ctx);
}
```

`selection`:
1. Take the batch of critical pairs with the smallest LCM degree (capped by `options.matrix_batch_size`, default 512).
2. `ctx.matrix.increase_upper_bound(largest LCM in batch)`.
3. For each pair `(f, g)` with LCM `L`: produce two rows, `get_u_g(L/lm(f), lm(f))` and `get_u_g(L/lm(g), lm(g))`, take their difference and put both candidate rows in the submatrix.
4. `reducer.reduce(submatrix, ...)` (cancellation + row/column limits enforced inside).
5. Nonzero rows go to `to_add`.

`update`:
1. For each `p ∈ to_add` (smallest first), make `p` monic, split into LM + tail.
2. `pm.insert_basis_polynomial(lm, tail)`.
3. Lead-reduce all still-pending `to_add` polynomials against the new basis element.
4. Remove critical pairs invalidated by the new lm (chain criterion + product criterion).
5. Form new critical pairs with every other basis element; prune by Gebauer–Möller (mirrors `m4gb.hpp:1591-1624`).
6. Detect and process any basis element whose lm is now divisible by the new lm (the `IMMEDIATE_BASIS_REDUCE` path in `m4gb.hpp:1626-1663`).

### M4GB dispatch & preconditions

`can_use_m4gb(polys, variables, options)`:

- `options.modulus` is prime and `< 2^31`,
- `options.order == DegRevLex`,
- no symbolic parameters,
- the codec for `(variables.size(), max_input_degree + selection_headroom)` constructs without overflow.

Phase 4 supports **GF(p) + DegRevLex only**. For `algorithm == M4GB` with `modulus == 0`, return `UnsupportedCoefficientDomain`. For `Auto`, fall back silently.

---

## Auto Heuristic

After Phase 7 (all four algorithms available):

| Conditions | Choice |
|---|---|
| input ≤ 3 polynomials or ≤ 2 variables | Buchberger |
| `has_symbolic_parameters` | F5B (criteria-aware, supports `ExpressionCoeffDomain`) |
| `modulus != 0`, DegRevLex, ≥ 6 vars | M4GB |
| `modulus == 0` or `modulus != 0` with non-DegRevLex order | MoGVW |
| no signature benefit expected (tiny systems, lex on big systems) | Buchberger |

Staging:

- After Phase 2 (F5B only): `Auto` chooses Buchberger or F5B.
- After Phase 3 (F5B + MoGVW): add MoGVW for DegRevLex+QQ/GF(p) medium systems.
- After Phase 5 (all): full table above.

---

## Implementation Phases

### Phase 0 — Benchmark + dispatch scaffolding

**Tasks:**
1. Add `M4GB` to `GroebnerAlgorithm` (`F5B` and `MoGVW` are already in the enum).
2. Extend `GroebnerOptions` and `GroebnerStats` with the fields above.
3. Replace `groebner.cpp:388-425` with a `switch` on `options.algorithm`. Buchberger remains the only real implementation; new cases route to `unsupported(...)` for now.
4. Extend `benchmarks/groebner_bench.cpp` to print all new stat counters (zeroed by Buchberger so it just continues to work).

**Acceptance:**
- All existing tests in `symengine/tests/polynomial/test_groebner.cpp` pass unchanged.
- `benchmarks/groebner_bench.cpp` builds.
- Requesting `F5B`, `M4GB`, or `MoGVW` explicitly returns `UnsupportedCoefficientDomain` (until the corresponding phase wires it).

### Phase 1 — Shared math & codecs

**Tasks:**
1. Implement `GFpCoeffDomain` with extended-Euclidean `inv`.
2. Implement `RuntimeDegRevLexCodec` (runtime port of `monomial_degrevlex_intcodec::static_data_t`). Detect `multiset_coefficient` overflow at construction.
3. Implement `PackedMonomial` + `PackedMonomialHash` + `PackedMonomialLess` for each supported variable/exponent schedule.
4. Implement `LeadingMonomialIndex` (start: linear scan with degree filter).
5. Implement `MatrixRow<Coeff>` + `GroebnerMatrixReducer<Coeff, Domain>` with cancellation/limits.
6. Create `symengine/polys/groebner_signature.h` with `Signature`, `SignatureLess`, `signature_mul`, and a `static_assert`-style test that ordering is total and matches both F5B and MoGVW conventions.

**Acceptance:**
- Codec round-trip property tests: `from_index(to_index(m, d)) == m` for every monomial up to the chosen `(num_vars, max_degree)`.
- DegRevLex ordering produced by `PackedMonomialLess` matches `compare_monomial_less(..., DegRevLex)` on randomized inputs.
- `GFpCoeffDomain`: `inv(a) * a ≡ 1` for all `a ∈ [1, p)` on `p ∈ {2, 3, 7, 31, 2^31 - 1}`.
- `GroebnerMatrixReducer` over GF(31): reduces a hand-constructed full-rank 5×5 system to RREF.
- Worked-example signature comparisons from `JSC-F5B.tex` §3 and `mogvw.tex` §3.4 agree with `SignatureLess`.

### Phase 2 — F5B (signature-based Buchberger)

**Tasks:**
1. Implement `LabeledPoly`, `CriticalPair`, `lbp_sub`, `lbp_mul_term`, `make_critical_pair`, `s_polynomial`, and CP ordering — as per §"Implementing F5B" above.
2. Implement `is_rewritable_or_comparable` and `f5_reduce` mirroring `sympy/polys/groebnertools.py`.
3. Implement `pernici_pre_reduce` and the main F5B loop.
4. Wire `GroebnerAlgorithm::F5B` to `groebner_basis_f5b` for all three coefficient domains (`RationalCoeffDomain`, `GFpCoeffDomain`, `ExpressionCoeffDomain`).
5. Port the SymPy F5B fixtures (`test_groebner_f5b`, `test_benchmark_minpoly_f5b`, `test_benchmark_katsura3_f5b`, `test_benchmark_kastura_4_f5b` [slow], `test_benchmark_czichowski_f5b` [slow], `test_benchmark_cyclic_4_f5b` [slow]) into `symengine/tests/polynomial/test_groebner_f5b.cpp`.
6. Parametrize the existing `test_groebner.cpp` corpus to run through F5B as well as Buchberger and assert ideal equality of the reduced bases.
7. Add F5B to `benchmarks/groebner_bench.cpp` for every benchmark that currently runs Buchberger.

**Acceptance:**
- All ported SymPy fixtures pass.
- Every Buchberger-passing test in `test_groebner.cpp` also passes for F5B.
- Stat counters `rejected_by_syzygy + rejected_by_rewritten > 0` on Katsura-3, Katsura-4, and cyclic-4.
- F5B benchmarks on Katsura-3 / Katsura-4 / cyclic-4 show fewer `reductions_to_zero` than Buchberger on the same input.
- Cancellation triggers `Cancelled` from `f5_reduce` and the main loop.

### Phase 3 — Scalar MoGVW

**Tasks:**
1. Reuse the `Signature` infrastructure from Phase 1 / Phase 2. Add `LabeledMonomial`, `MoGVWContext`.
2. Implement the `mutualreduce` procedure faithfully (use the code in §"Implementing MoGVW" above).
3. Implement the main loop with `liftdeg` reset.
4. Implement LCM / principal-syzygy / rewritten criteria with separate stat counters.
5. Wire `GroebnerAlgorithm::MoGVW` to `groebner_basis_mogvw_impl`.
6. Run final reduced-basis extraction through the existing `interreduce`/`reduce` helpers in `groebner.cpp`.

**Acceptance:**
- Every test in `test_groebner.cpp` that currently passes for Buchberger also passes for MoGVW (rational coefficients only — the test corpus does not yet exercise GF(p)).
- Worked Example 3.4 in `mogvw.tex` produces the same `r_1, ..., r_5` sequence (drives a criterion-counter assertion).
- Property test: `is_groebner(MoGVW output, variables, order) == true` on Katsura-3, Katsura-4, and the toric/binomial fixture from `benchmarks/groebner_bench.cpp`.
- MoGVW and F5B agree on every shared test case (ideal equality of bases — implementation cross-check).
- Cancellation triggers `Cancelled` from inside `mutual_reduce`, `top_reduce`, and the main lift loop.
- `verify_with_buchberger == true` runs both algorithms and asserts ideal equality.

### Phase 4 — M4GB polymatrix core

**Tasks:**
1. Implement `PolyMatrix` with eager generation.
2. Implement `increase_upper_bound` (simple scan), `insert_basis_polynomial`, `get_u_g`, `shrink` (eager).
3. Implement `assert_invariants` and call it in debug builds after every mutating method.

**Acceptance:**
- Unit tests on hand-constructed polymatrices over GF(7):
  - Inserting `{x² - y, y² - x}`, increasing upper bound to `x³`, and reading every entry produces tail-reduced multiples.
  - `get_u_g(x, x² - y)` returns `x³ - xy` reduced through `matrix[x³]`.
  - `shrink` after inserting a new basis poly correctly removes columns and rewrites rows.
- Invariants hold across 1000 randomized insert/increase/shrink sequences.

### Phase 5 — M4GB solver loop

**Tasks:**
1. Implement `selection`, `update`, Gebauer–Möller CP pruning, and the IMMEDIATE_BASIS_REDUCE path.
2. Wire `GroebnerAlgorithm::M4GB` to `groebner_basis_m4gb_impl` (GF(p) + DegRevLex only).
3. Extend `benchmarks/groebner_bench.cpp` with the `GF(31)` MQ inputs from `/opt/m4gb/testdata/gf31/` (parse on the fly — do not vendor the files).

**Acceptance:**
- On every small GF(p) system added to the test corpus, M4GB output is a reduced Groebner basis (`is_reduced_basis == true`) and ideal-equivalent to the Buchberger / F5B output reduced modulo `p`.
- `FGLM` can convert an M4GB DegRevLex output to Lex on a zero-dimensional input.
- Dense MQ over GF(31) at `n = 6` completes within the benchmark time budget.
- Cancellation triggers `Cancelled` from inside `get_u_g`, `selection`, `update`, and the reducer.
- `max_matrix_rows`, `max_matrix_columns`, `max_degree` trigger `ResourceLimitExceeded`.

### Phase 6 — ExpressionCoeffDomain hardening + parametric test corpus

The `ExpressionCoeffDomain` path is already wired up in `groebner.cpp:386-407` for the symbolic-coefficient case; this phase is about making it robust and bringing in real-world parametric test cases — not adding a new domain.

**Tasks:**

1. **Coefficient-domain hardening.** `ExpressionCoeffDomain` currently calls `expand(...)` on every `add`/`sub`/`mul`/`div` (`groebner_internal.h:121-142`), with no normalization between operations. For parametric inputs whose intermediate coefficients contain common factors, this leads to monotonic coefficient bloat. Improvements:
   - Add an internal `normalize` step (call `expand` + optional `simplify`/`gcd` factoring on the leading coefficient before scaling for monic conversion).
   - Detect coefficients that are *exactly zero* vs. *generically zero* (e.g. `K1 - K1` after `expand` is zero; `K1*K2 - K2*K1` after `expand` is zero; but `K1**2 - K2*K3` is generically nonzero). The existing `is_zero` predicate uses `eq(*c.get_basic(), *SymEngine::zero)` which is sufficient *only* after `expand`; document this and verify all `is_zero` callers post-`expand`.
   - Wire `options.cancellation_token` into the inner `expand` chain so cancellation propagates during expensive symbolic operations. Today the token is only checked at the Buchberger loop level.

2. **Genericity policy.** Document explicitly in the header that for `ExpressionCoeffDomain`:
   - The returned basis is correct **under the assumption that no parameter satisfies a leading-coefficient relation that would have triggered an alternate branch**. This matches what Maple/Mathematica do by default.
   - Side-conditions (the non-vanishing leading-coefficient products encountered during the run) should be exposed via a new `GroebnerStats::genericity_assumptions: std::vector<RCP<const Basic>>` field (each entry is a polynomial that was *assumed* nonzero). Tests assert these are recorded.
   - Building a Comprehensive Groebner System (case-splitting on each branch) is **out of scope**; recording the assumptions is the cheap, useful subset.

3. **Port the symcse parametric test corpus** into `symengine/tests/polynomial/test_groebner_parametric.cpp`. Each test asserts: basis size, ideal-equivalence to the input (every input poly reduces to zero modulo the basis), and reduced-basis canonical form. Sources, in order of increasing system size:

   - From `symcse/tests/polynomials/test__groebner_sage.py`:
     - `test_call_groebner_sage__1`: `[x**2 - y, x**3 - x]` with `vars = [x, y]`. Pure-rational, no parameters; useful as a sanity baseline that also matches Buchberger in `lex`, `degrevlex + fglm`, and `revlex + fglm`. Three test variants (one per ordering).
     - `test_call_groebner_sage__2`: `[x**2 + K1*x*y, x*y + 2*y**3 - K2]` with `vars = [x, y]`, `params = [K1, K2]`. **Smallest parametric case** — good first test for `ExpressionCoeffDomain` correctness, and produces a basis where the `K1*K2` cross-product terms exercise coefficient expansion.

   - From `symcse/tests/polynomials/test__singular_groebner.py`:
     - `TestCase1`: 3 vars `[z, y, x]`, 2 params `[K1, ux]`, 3 polynomials. Includes a `K1*(ux - x)` factor that becomes a stress test for the `expand` step.
     - `TestCase2`: 3 vars, 4 params `[K1, b0, b1, ux]`, 3 polynomials. Asserts `len(basis_raw) == 3`.
     - `TestCase3`: 9 vars, 8 params, 9 polynomials (chemical equilibrium with augmentation variable `zaug`). Asserts `len(degrevlex_basis) >= 5`. This is the largest parametric system in the corpus and the most useful **benchmark**.

   - From `symcse/tests/polynomials/test__groebner.py`:
     - `MyPolySys1` (ammonia equilibrium): 5 vars `[NH4p, NH3, Hp, H2O, OHm]`, 5 params `[Ka, Kw, iNH3, iOHm, iHp]`, 5 polynomials. Tests the full `mk_aug_gbasis_zi(zi)` flow for `zi ∈ {0, 1, 2, 3, 4}` — five distinct augmented systems. Worth porting because it exercises `augmented_groebner_basis` (already in `groebner.cpp:834`) with symbolic coefficients.

   Port format: hand-translate each system into a C++ test using SymEngine `symbols(...)` + `parse(...)`. Do not pull in a Python runtime. Where the Python tests assert basis *size* only, the C++ tests should additionally assert reducedness and ideal-equivalence.

4. **Add the parametric corpus to `benchmarks/groebner_bench.cpp`** as a new `parametric_benchmarks()` group. For each system, time:
   - Buchberger,
   - F5B (Phase 2 wires this domain in; expect criterion-based reduction count savings even when end-to-end wall time is similar),
   - MoGVW (after Phase 3, falling back to Buchberger because MoGVW does not support `ExpressionCoeffDomain` in this plan — verify the fallback path takes the same time, modulo dispatch overhead).

   Emit per-benchmark: total time, basis size, terms per polynomial (proxy for coefficient bloat), `genericity_assumptions.size()`, and the rejection / `f5_reductions` counters.

**Acceptance:**
- All ported tests pass against the current Buchberger + `ExpressionCoeffDomain` path. If any do not, that's an existing correctness bug worth fixing in this phase.
- `TestCase3` (the largest) completes in under a generous time budget (e.g. 30 s on the CI machine). If it exceeds that, document it as a known limitation and tag it `[!benchmark]` so it doesn't block the test suite.
- `genericity_assumptions` is non-empty for `test_call_groebner_sage__2` and `TestCase1`/`TestCase2`/`TestCase3` (each branches on at least one leading-coefficient assumption).

### Phase 7 — Auto + cleanup

**Tasks:**
1. Implement `resolve_algorithm` per the heuristic table.
2. Expose `verify_with_buchberger` in the benchmark output (asserts on mismatch).
3. Document the public API additions in the header.
4. Add a cross-algorithm regression sweep: for every test in `test_groebner.cpp` and `test_groebner_f5b.cpp`, run all four algorithms whose preconditions are satisfied and assert pairwise ideal equality.

---

## Extension Paths

Concrete follow-ups after this plan lands, in roughly increasing engineering cost. None of these are committed work; they're the directions the architecture is shaped to allow.

### A. M4GB over Q (rational, no parameters)

**Effort:** medium. **Value:** medium for users without GF(p) inputs.

Reuse the entire `PolyMatrix` + `get_u_g` + selection / update architecture from Phases 3-4, swapping `uint64_t` + `GFpCoeffDomain` for `rational_class` + `RationalCoeffDomain`. The polymatrix invariants are field-agnostic. The performance benefit vs. Buchberger shrinks because rational arithmetic is no longer O(1), but the cached-multiples idea still amortizes redundant reductions. Main risk: rational coefficient growth in the dense vectors. Mitigation: normalize denominators eagerly and shrink trailing zeros aggressively.

### B. MoGVW over Q(C₁, …, Cₖ) — signature criteria over an expensive field

**Effort:** medium-high. **Value:** high for parametric workloads.

This is the *right* match for symbolic-coefficient Groebner. Signature criteria (LCM, principal syzygy, rewritten) prune *reductions* — and over `Q(C₁,…,Cₖ)` each reduction is expensive precisely because of the symbolic coefficient. Pruning earlier is more valuable than amortizing later. The matrix tricks of M4GB are the wrong investment here: caching dense vectors of rational functions is more expensive than recomputing them sparsely.

To make this work:
- Extend `MoGVWContext` to instantiate over `(Expression, ExpressionCoeffDomain)` instead of `(uint64_t, GFpCoeffDomain)`.
- Tighten the genericity-assumption tracking (see Phase 6) so the signature comparisons and zero-tests have a documented contract.
- The "smaller signature is preferred" invariant only works when we can decide signature equality / divisibility cheaply. Signatures live in the monomial part, not the coefficient part, so this is fine even with symbolic coefficients.

This is the path I recommend if the medium-term goal is parametric Groebner performance.

### C. Comprehensive Groebner System (CGS) on top of (B)

**Effort:** high. **Value:** high for users who need correctness on parameter-space boundaries.

Instead of recording the genericity assumptions as side conditions, branch on each one and compute a separate basis for each parameter-space cell. Produces a *system* of bases with associated conditions. This is what Magma / Maple's `Groebner[Solve]` do for "really" parametric inputs. It is a project in its own right and should be designed after (B) is in.

### D. Multi-prime CRT for Q from GF(p) M4GB

**Effort:** medium-high. **Value:** medium — competes with B for rational inputs without parameters.

Run M4GB over several primes `p₁, p₂, …`, reconstruct the rational basis via CRT + rational reconstruction. The performance argument is that several GF(p) M4GB runs are cheaper than one Q-Buchberger run for medium-size systems. Requires bad-prime detection (a prime is "bad" if some leading coefficient vanishes mod p) and a coefficient-stable basis comparator across primes.

### E. GF(2) and GF(2^n) support

**Effort:** medium. **Value:** depends on workload (high for cryptographic MQ, low otherwise).

For GF(2), exponents satisfy `x^2 = x` so the polynomial ring is effectively `F_2[x_1, ..., x_n] / (x_i^2 - x_i)`; this lets you bit-pack entire rows. Worth doing if Boolean MQ benchmarks (HFE, Rain) become a target. For GF(2^n), the field arithmetic switches to polynomial-basis multiplication tables, but the polymatrix architecture is otherwise unchanged.

### F. M4GB over an arbitrary fraction field — explicitly **not** recommended

The polymatrix's value proposition is amortizing field operations that are individually cheap. Over `Q(C₁, …, Cₖ)`, each field operation is itself a polynomial computation, and the cached dense vectors carry that cost forward instead of saving it. Even where M4GB beats Buchberger by an order of magnitude over GF(p), it would likely *lose* to Buchberger over a fraction field because the dense representation wastes memory on slots that hold large rational functions. If you need symbolic-coefficient performance, do (B) instead.

---

## Risk Register & Mitigations

- **Signature correctness (F5B / MoGVW).** Mitigation: F5B Phase 2 is a near-line-for-line port of SymPy's well-tested `_f5b` and is validated against SymPy fixtures; MoGVW Phase 3 reuses the same Signature/ordering machinery and is cross-checked against F5B for ideal equality. Per-criterion stat counters force regressions to be visible; `verify_with_buchberger` toggle covers the entire test corpus.
- **Polymatrix invariants (M4GB).** Mitigation: `assert_invariants` after every mutator in debug builds; randomized sequence test; eager generation in Phase 4 (avoid the lazy-shrink edge cases until Phase 5).
- **Monomial overflow.** Mitigation: codec constructor checks `multiset_coefficient` overflow; `PackedMonomial` schedule chosen at runtime and rejects out-of-range exponents (explicit M4GB returns `UnsupportedCoefficientDomain`; `Auto` falls back).
- **Memory blowout.** Mitigation: `max_matrix_rows` / `max_matrix_columns` / `max_degree` honored inside `GroebnerMatrixReducer::reduce`, `get_u_g`, `increase_upper_bound`, and `selection`; check `cancellation_token` every `cancellation_check_interval` iterations.
- **Output normalization drift.** Mitigation: every algorithm's output flows through the existing `interreduce`/`reduce` helpers before `gpoly_to_basic` so canonicalization is shared. Verify with `is_reduced_basis`.
- **Reference-code license / provenance.** Mitigation: re-author every component for SymEngine style; preserve a header-level provenance note in `groebner_m4gb.cpp` pointing at `ISSAC17-MS-M4GB.pdf` and the `/opt/m4gb` algorithm logic, but do not paste source verbatim.

---

## Out of scope (explicit) for this plan

- Matrix-MoGVW (one-sided F4-style elimination of signature-grouped rows). See Extension Path B.
- Lazy generation / lazy shrink in M4GB beyond the optional Phase-5 follow-up.
- GF(2^n), GF(2), GF(p) for `p ≥ 2^31`. See Extension Path E.
- Rational coefficients in M4GB. See Extension Path A.
- Modular reconstruction for QQ from multiple GF(p) runs. See Extension Path D.
- MoGVW over symbolic/parametric coefficients. See Extension Path B.
- Comprehensive Groebner Systems with explicit case-splitting on parameters. See Extension Path C (Phase 6 records side-conditions but does not branch on them).
- Noncommutative bases.
