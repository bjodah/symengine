# Plan: upstream SymEngine core gaps found by symcsepp

## Objective

Upstream only the gaps from symcsepp's symbolic-core port that cannot be
implemented correctly as consumer-side utilities. Keep the proposed changes
generic to SymEngine and split them into independently reviewable pull
requests.

The
[companion symcsepp plan](../../../../symcsepp/.meta-docs/plans/symengine-gap-utilities.md)
owns expression decomposition, occurrence counting, CSE policy/naming, the
required `powsimp` subset, and symbolic `cosm1`. Those APIs should not be added
to SymEngine as part of this effort.

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

## Pull request 1: diagnostic linear-system results

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

Names may be adjusted during upstream review, but the status must be explicit
and machine-readable. Do not use an empty solution to encode three outcomes.
The existing `linsolve` functions delegate to the detailed implementation and
retain their current compatibility behavior: return the solution for
`unique`, and `{}` otherwise.

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

## Pull request 2: safe external `FunctionWrapper` introspection

### Problem to solve

`FunctionWrapper` is an extension point, but subclasses cannot allocate a new
`TypeID`. Consequently `is_a<T>`'s exact-type documentation does not hold when
`T` is an external wrapper subclass: the inherited static type code is shared
by every wrapper.

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
   exact RTTI check (or a design with equivalent correctness). It should be
   available only when SymEngine is built with RTTI. Keep `is_a_sub<T>` for
   ordinary inheritance checks.
3. Do not silently change all `is_a<T>` calls to RTTI. Its type-code fast path
   is pervasive, and a global semantic/performance change is unnecessary.
4. If upstream maintainers prefer no new helper, the minimum acceptable change
   is prominent API documentation plus tests demonstrating `is_a_sub<T>` as
   the supported external-wrapper predicate. In that outcome, explicitly
   document that exact discrimination requires `typeid`/`dynamic_cast`.

### Tests

Define two wrapper subclasses in `test_functions.cpp` and verify that the
supported predicate accepts only the requested subclass, that
`is_a<FunctionWrapper>` accepts both, and that base-class `is_a_sub` behavior
is unchanged. Include the no-RTTI build behavior in compile-time guards.

Acceptance: extension authors have a documented, tested predicate that cannot
mistake one custom wrapper for another.

## Pull request 3: make `srepr` identify function wrappers

Teach `SreprPrinter` to include function identity rather than emitting only
the shared `FunctionWrapper` type code and arguments.

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
   require escaping.
5. Note the output change in release notes because consumers may use the
   string as a deterministic key, even though it is not a serialization
   format.

Acceptance: `srepr` preserves the public identity of every function-symbol
node and remains deterministic across repeated calls.

## Explicit non-goals

- Do not add the SymPy-shaped coefficient/decomposition helpers. SymEngine's
  `Add`/`Mul` accessors are sufficient for consumer utilities.
- Do not add `Expr.count`; traversal policy belongs to the consumer.
- Do not add CSE `ignore` or symbol-stream APIs for symcsepp. Its required
  semantics can be implemented by filtering, inlining, and renaming around
  the existing CSE result.
- Do not upstream a partial `powsimp` based only on symcsepp needs. It can be
  reconsidered separately if SymEngine defines assumptions and branch-safety
  requirements for a general API.
- Do not add `cosm1` as a core type. A `FunctionWrapper` is the intended
  extension mechanism.
- Do not design a general `doit()` evaluation protocol in this effort.
- Do not include the Eigen permutation issue; it is not a SymEngine defect.

## Upstream sequence and verification

Submit the three pull requests independently; none should depend on another.
Each pull request should reference a minimal standalone SymEngine reproducer,
not the symcsepp implementation.

For every pull request:

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

This plan is complete when the accepted upstream APIs are released or pinned
by submodule commit, symcsepp consumes the detailed solve result, and the
upstream log records rejected or superseded proposals rather than continuing
to describe them as missing features.
