# Plan: SymEngine fork changes needed by symcsepp (sep24 branch)

Decision (2026-07-06, Björn): no upstreaming to `symengine/symengine`. The
vendored fork (`bjodah/symengine`, branch `sep24`) is fixed directly for our
purposes. Changes are still written in an upstreamable STYLE (generic,
self-contained, formatted, tested) — not in order to submit them, but so the
fork's divergence stays small and mechanical to carry across a future rebase
onto upstream master.

## Objective

Implement, directly on the fork, only the gaps from symcsepp's symbolic-core
port that cannot be implemented correctly as consumer-side utilities. Keep
the changes generic to SymEngine and split them into independently
reviewable commits (one merge to `sep24` per change).

The
[companion symcsepp plan](../../../../symcsepp/.meta-docs/plans/symengine-gap-utilities.md)
owns expression decomposition, occurrence counting, CSE policy/naming, the
required `powsimp` subset, and symbolic `cosm1`. Those APIs must NOT be added
to the fork: every SymPy-shaped convenience kept out of SymEngine is
divergence we do not have to carry across a future rebase.

## Baseline audit

The upstream log was written against an older understanding of the API. In
the current vendored tree:

- `SreprPrinter` exists in `symengine/printers/srepr.h`, but its generic
  representation of a `FunctionWrapper` omits `FunctionSymbol::get_name()`;
- `UnevaluatedExpr` exists in `functions.h`; there is no general `doit()`
  protocol, and symcsepp does not require one;
- `linsolve` returns only `vec_basic`; a failed fraction-free solve becomes an
  empty vector, conflating inconsistent and underdetermined systems;
- all external `FunctionWrapper` subclasses inherit the same SymEngine type
  code. `dynamic_cast`/`is_a_sub` can discriminate them when RTTI is enabled,
  while `is_a<ExternalWrapper>` is unsafe because `is_a` compares the inherited
  type code.

Linear solving and wrapper introspection require SymEngine-owned API design.
The printer defect requires a core correction because a structural
representation must know about core node interfaces. `UnevaluatedExpr`/`doit()`
is explicitly out of scope.

## Step 0: record the divergence being added (divergence ledger)

The fork's history was squash-imported (`git log` attributes essentially the
whole tree to one commit, `42d05d2`, plus two LLVM-related commits on top),
so `git log`/`git blame` inside the submodule cannot separate existing
fork-local features from upstream ones. Snapshot the current delta once so
the new work is added to a known baseline:

```bash
cd /work/external/symengine
git remote add upstream https://github.com/symengine/symengine.git
git fetch upstream
git diff --stat upstream/master...HEAD -- symengine/
```

Start a short ledger at `.meta-docs/fork-divergence.md`: one line per
fork-local feature (known so far: the srepr printer, the `int`-returning
`fraction_free_gauss_jordan_solve` plus `linsolve`'s `{}`-on-failure, the
`#if HAVE_SYMENGINE_RTTI` guard around `is_a_sub`, the LLVM customization
points, and `UnevaluatedExpr` if the diff shows it is ours), plus one line
per change below as it lands. This ledger is what turns a future rebase onto
upstream master into a checklist instead of an archaeology project.

Compatibility baseline for all three changes is the CURRENT FORK behavior —
that is what symcsepp consumes. Whether upstream behaves differently is not
a delivery question anymore; it is only a note for the ledger.

## Change 1: diagnostic linear-system results

### Public contract

Add a non-throwing detailed API while preserving the existing `linsolve`
overloads:

```cpp
enum class LinearSolveStatus {
    unique,
    underdetermined,
    inconsistent
};

struct LinearSolveResult {
    LinearSolveStatus status;
    vec_basic solution;       // populated only for `unique`
    unsigned coefficient_rank;
    unsigned augmented_rank;
    std::vector<unsigned> pivot_columns;
};

LinearSolveResult linsolve_detailed(const DenseMatrix &system,
                                    const vec_sym &syms);
LinearSolveResult linsolve_detailed(const vec_basic &system,
                                    const vec_sym &syms);
```

The status must be explicit and machine-readable. Do not use an empty
solution to encode three outcomes. The existing `linsolve` functions
delegate to the detailed implementation and retain their current fork
behavior — return the solution for `unique`, and `{}` otherwise — because
symcsepp's `symb_solve_fully` consumes exactly that today.

### Current state (verified in the vendored tree, 2026-07-06)

- `solve.cpp`: `linsolve_helper` treats any nonzero return of
  `fraction_free_gauss_jordan_solve` as "return `{}`"; the fail code is
  `i + 1` for the first column without a structural pivot — exactly the
  underdetermined/inconsistent conflation this change removes.
- `dense_matrix.cpp`: `fraction_free_gauss_jordan_solve` is square-only
  (`SYMENGINE_ASSERT(A.row_ == A.col_)`, compiled out without
  `WITH_SYMENGINE_ASSERT`, so rectangular input in a release build indexes
  out of bounds instead of failing cleanly). `DenseMatrix::rank()` throws
  `NotImplementedError`, confirming the "do not call `rank()`" instruction
  below.
- Candidate building blocks already declared in `matrix.h`:
  `pivoted_fraction_free_gauss_jordan_elimination` and
  `reduced_row_echelon_form`. Prefer factoring pivot/rank extraction on top
  of `reduced_row_echelon_form` of the augmented matrix over writing a new
  eliminator — but add rectangular-input tests for the chosen routine first,
  since its existing callers are square-only.
- symcsepp integration point: `symcsepp/symb_eqsys/lin_eq.cpp`,
  `LinEqGroup::symb_solve_fully`, currently treats an empty `linsolve`
  result as the "fall back to matchprob + pivoted LU" signal. After this
  change: `underdetermined` -> the existing partial-solve fallback,
  `inconsistent` -> the error path; that mapping is the integration test
  described at the end of this plan.

### Elimination work

1. Implement fraction-free row reduction for rectangular coefficient and
   augmented matrices, or factor the current square-only routine into an
   internal elimination result that exposes pivots and inconsistent rows.
2. Compute structural ranks from pivots; do not call `DenseMatrix::rank()`,
   which is currently unimplemented.
3. Classify a row as inconsistent only when every coefficient is structurally
   zero and its right-hand side is structurally nonzero. A consistent system
   with fewer pivots than variables is underdetermined.
4. Populate `solution` only when rank equals the number of unknowns. Support
   overdetermined-but-consistent unique systems as well as square systems.
5. Document symbolic-pivot semantics: an expression not structurally equal to
   zero is treated as a generic nonzero pivot, so parameter-specialized rank
   drops are not inferred without an assumptions system.
6. Preserve rejection of nonlinear input in `linear_eqns_to_matrix` and
   validate matrix dimensions in release builds rather than relying only on
   assertions.

### Tests

Extend `symengine/tests/basic/test_solve.cpp` (and matrix tests if elimination
is factored there) with:

- square and overdetermined unique systems;
- a consistent rank-deficient system;
- an inconsistent system with the same coefficient rank as the preceding
  case;
- zero rows and row permutations;
- rectangular systems with fewer equations than variables;
- symbolic coefficients, including a documented generic-nonzero pivot;
- compatibility checks that legacy `linsolve` still returns `{}` for both
  non-unique statuses;
- dimension errors and nonlinear-equation errors.

Acceptance: callers can distinguish underdetermined from inconsistent input
without parsing exceptions or recomputing a rank, and existing source/API
behavior remains compatible.

## Change 2: safe external `FunctionWrapper` introspection

### Problem to solve

`FunctionWrapper` is an extension point, but subclasses cannot allocate a new
`TypeID`. Consequently `is_a<T>`'s exact-type documentation does not hold when
`T` is an external wrapper subclass: the inherited static type code is shared
by every wrapper.

Verified in the vendored tree: `basic-inl.h` implements `is_a` as
`T::type_code_id == b.get_type_code()`; symcsepp's `CustomFunc` subclasses
(`Expm1`, `Log1p`, ...) inherit `FunctionWrapper`'s type code, so
`is_a<Expm1>` returns true for ANY wrapper. `is_a_sub` is
`dynamic_cast`-based and (in the fork) guarded by `#if HAVE_SYMENGINE_RTTI`
(a CMake cache option defaulting to yes). Note that `dynamic_cast<const T*>`
alone gives SUBCLASS semantics — the exact-type helper below needs
`typeid(value) == typeid(T)`.

### API and documentation change

1. Document on `FunctionWrapper` and `is_a` that external subclasses must not
   be tested with `is_a<Derived>`.
2. Add a named helper in `functions.h`, implemented after `FunctionWrapper` is
   complete, for exact wrapper checks. A representative contract is:

   ```cpp
   template <class T>
   bool is_a_function_wrapper(const Basic &value);
   ```

   It must statically require `T` to derive from `FunctionWrapper` and use an
   exact RTTI check — concretely `typeid(value) == typeid(T)`, not a bare
   `dynamic_cast`, which would accept further subclasses of `T`. It should be
   available only when SymEngine is built with RTTI. Keep `is_a_sub<T>` for
   ordinary inheritance checks.
3. Do not silently change all `is_a<T>` calls to RTTI. Its type-code fast path
   is pervasive, and a global semantic/performance change is unnecessary.

### Tests

Define two wrapper subclasses in `test_functions.cpp` and verify that the
supported predicate accepts only the requested subclass, that
`is_a<FunctionWrapper>` accepts both, and that base-class `is_a_sub` behavior
is unchanged. Include the no-RTTI build behavior in compile-time guards.

Acceptance: extension authors have a documented, tested predicate that cannot
mistake one custom wrapper for another. symcsepp migration site once
available: the `dynamic_cast` chain in `symcsepp/symb_eqsys/scalar_solve.cpp`
(origin-invertible detection). The name-based dispatch in symcsepp's
LLVM/eval visitors is deliberately name-based and should NOT migrate.

## Change 3: make `srepr` identify function wrappers (do this first)

Teach `SreprPrinter` to include function identity rather than emitting only
the shared `FunctionWrapper` type code and arguments.

Do this change FIRST: it is small, and the companion symcsepp plan's work
package 4 (stable structural keys) is blocked until it is merged to `sep24`
and the submodule pinned in the super-project.

Verified defects in the current printer (`printers/srepr.cpp`): the generic
`bvisit` prints `type_code_name(get_type_code())` plus arguments only, so
`FunctionSymbol("f", x)` and `FunctionSymbol("g", x)` render identically and
every external wrapper renders as `FunctionWrapper(...)` — `Expm1(x)` and
`Log1p(x)` collide. Symbol names are quoted but not escaped
(`s << '"' << get_name() << '"'`), so a name containing `"` breaks key
uniqueness.

1. Add printer handling for `FunctionSymbol` and `FunctionWrapper` that
   includes `get_name()` and recursively prints arguments. Introduce a shared
   string-escaping helper and use it for both function names and symbol names.
2. Choose and document a stable format, for example
   `FunctionWrapper("Expm1", Symbol("x"))`. The format need not expose the
   C++ subclass name, which is neither portable nor available without RTTI;
   the symbolic function name is the relevant identity.
3. Ensure two wrappers with equal arguments and different names produce
   different representations.
4. Add cases to `symengine/tests/printing/test_srepr.cpp` for ordinary
   `FunctionSymbol`, two custom wrappers, multiple arguments, and names that
   require escaping. Also pin the representation of two distinct `Dummy`s
   with equal names (decide whether the dummy index participates, and
   document the choice — symcsepp's stable-key consumer feeds
   `cse_no_floats`-generated dummies through this printer).
5. Note the output change in the divergence ledger and the commit message:
   consumers may use the string as a deterministic key, even though it is
   not a serialization format.

Acceptance: `srepr` preserves the public identity of every function-symbol
node and remains deterministic across repeated calls.

## Explicit non-goals

- Do not submit anything to `symengine/symengine`; upstreaming is explicitly
  deferred (possibly forever). Writing in an upstreamable style is a
  rebase-cost measure, not a submission plan.
- Do not add the SymPy-shaped coefficient/decomposition helpers. SymEngine's
  `Add`/`Mul` accessors are sufficient for consumer utilities.
- Do not add `Expr.count`; traversal policy belongs to the consumer.
- Do not add CSE `ignore` or symbol-stream APIs for symcsepp. Its required
  semantics can be implemented by filtering, inlining, and renaming around
  the existing CSE result.
- Do not add a partial `powsimp` to the fork based only on symcsepp needs.
  It can be reconsidered separately if SymEngine defines assumptions and
  branch-safety requirements for a general API.
- Do not add `cosm1` as a core type. A `FunctionWrapper` is the intended
  extension mechanism.
- Do not design a general `doit()` evaluation protocol in this effort.
- Do not include the Eigen permutation issue; it is not a SymEngine defect.

## Sequence and verification

Land the three changes as independent merges to `sep24`; none depends on
another technically. Recommended order: change 3 first (smallest; unblocks
symcsepp work package 4), then change 2 (small), then change 1 (largest).
Each change carries its own tests inside the fork and stays free of
symcsepp specifics.

For every change:

```bash
cmake -B build -S . -DWITH_SYMENGINE_ASSERT=yes -DBUILD_TESTS=yes
cmake --build build --parallel
ctest --test-dir build --output-on-failure
./bin/test_format_local.sh
```

After updating the submodule in the super-project, rebuild symcsepp and change
its `symb_solve_fully` fallback to consume `linsolve_detailed`. Add integration
tests asserting that SymEngine's `underdetermined` and `inconsistent` statuses
remain distinct through the symcsepp error/reporting layer.

This plan is complete when all three changes are merged to `sep24` and
pinned by the super-project submodule, symcsepp consumes the detailed solve
result, the divergence ledger lists every fork-local feature including these
three, and the upstream log in `03-symb-eqsys-symbolic-core.md` records the
fork-fix disposition instead of describing these items as missing features.
