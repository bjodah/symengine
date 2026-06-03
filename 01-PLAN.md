# Gruntz Algorithm Implementation for SymEngine

## Overview

Implement the Gruntz algorithm for computing symbolic limits in SymEngine. The algorithm reduces any limit (`x → a`, `x → ±∞`) to the `x → ∞` case, classifies subexpressions by growth rate, introduces a small variable `w → 0`, and reads off the limit from the leading term of a series expansion.

Reference implementation: SymPy's `/opt-3/cpython-v3.13-apt-deb/lib/python3.13/site-packages/sympy/series/gruntz.py`

---

## Files to Create

| File | Purpose |
|------|---------|
| `/work/symengine/gruntz.h` | Header: `SubsSet` class, `gruntz()` and helper function declarations |
| `/work/symengine/gruntz.cpp` | ~800-1000 line implementation of the full algorithm |
| `/work/symengine/tests/basic/test_gruntz.cpp` | Test suite ported from SymPy |

## Files to Modify

| File | Change |
|------|--------|
| `/work/symengine/CMakeLists.txt` | Add `gruntz.cpp` to `SRC` list; add `gruntz.h` to `HEADERS` list |
| `/work/symengine/tests/basic/CMakeLists.txt` | Add `test_gruntz` executable target |

---

## Data Structures

### `SubsSet`

Maps MRV expressions to dummy variables, with rewrites tracking expression dependencies:

```cpp
struct SubsSet {
    // expr -> dummy variable
    std::map<RCP<const Basic>, RCP<const Symbol>, RCPBasicKeyLess> dict;
    // dummy -> rewritten expression (in terms of other dummies)
    std::map<RCP<const Symbol>, RCP<const Basic>, RCPBasicKeyLess> rewrites;
    
    RCP<const Symbol> insert(const RCP<const Basic> &expr);
    bool meets(const SubsSet &other) const;
    std::pair<SubsSet, RCP<const Basic>>
        union_set(const SubsSet &other, const RCP<const Basic> &exps) const;
    RCP<const Basic> do_subs(const RCP<const Basic> &e) const;
};
```

- `insert()`: Auto-create `Dummy()` if `expr` not yet in `dict`
- `meets()`: True if key sets intersect
- `union_set()`: Merge two sets, remapping duplicates to same dummy; adjust `exps` accordingly
- `do_subs()`: Reverse-substitute dummies → original expressions

---

## Functions to Implement (all in `namespace SymEngine`)

### 1. `int sign(const RCP<const Basic> &e, const RCP<const Symbol> &x)`

Returns `1`, `-1`, or `0` for sign of `e(x)` as `x → ∞`.

**Logic:**
- Constants: use `is_positive()/is_negative()/is_zero()`
- `x` itself → `+1`
- `Mul`: split factors, multiply signs
- `exp(...)` → `+1`
- `Pow(E, ...)` → `+1`
- `Log(pos_arg)` → `sign(arg - 1)`
- Fallback: `mrv_leadterm(e, x)` → `sign(c0, x)`

### 2. `char compare(const RCP<const Basic> &a, const RCP<const Basic> &b, const RCP<const Symbol> &x)`

Compares growth rates: returns `'<'`, `'='`, or `'>'`.

**Logic:**
- `la = log(a)`, `lb = log(b)`; special-case `exp` to use its exponent directly
- Compute `L = limitinf(la / lb, x)`
- `L == 0` → `'<'` (a grows slower)
- `L` is infinite → `'>'` (a grows faster)
- otherwise → `'='`

### 3. `std::pair<SubsSet, RCP<const Basic>> mrv(const RCP<const Basic> &e, const RCP<const Symbol> &x)`

**Most Rapidly Varying set.** Returns `(Omega, exps)` where `Omega` maps MRV subexpressions to dummies, and `exps` is `e` with those dummies substituted.

**Type dispatch:**
| Expression Type | Action |
|----------------|--------|
| No `x` dependency | `(empty, e)` |
| `e == x` | `({x: d1}, d1)` |
| `Mul`/`Add` | Split into independent/dependent parts; recurse on each; `mrv_max1` to combine |
| `Pow` (base != E) | If exponent has `x`: rewrite as `exp(exp*log(base))`, recurse. Else recurse on base, raise to exponent. |
| `Log` | Recurse on argument, wrap result in `log` |
| `Exp` / `Pow(E, ...)` | Check `limitinf(exp, x)`. If → ∞: create entry, recurse on exponent, merge with `mrv_max3`. Else recurse on exponent. |
| Other functions | Recurse on args (single-var only) |

### 4–5. `mrv_max3` / `mrv_max1`

Combine two MRV sets. `mrv_max3` uses `compare()` to pick the dominant set. `mrv_max1` computes the union first then calls `mrv_max3` with both orderings.

### 6–7. `moveup2` / `moveup`

When `x` itself is in the MRV set, replace `x → exp(x)` everywhere to "move up" to the exponential level.

### 8. `std::pair<RCP<const Basic>, RCP<const Basic>> rewrite(const RCP<const Basic> &e, const SubsSet &Omega, const RCP<const Symbol> &x, const RCP<const Symbol> &wsym)`

The most complex step. Rewrites `e` in terms of `w` (where `w → 0`):

1. Build dependency tree from `Omega.rewrites`, sort topologically
2. Pick representative `g` from Omega
3. If `g.exp → +∞`, set `wsym = 1/w` (so `w → 0`)
4. For each Omega element: `c = limitinf(f.exp / g.exp, x)`, rewrite as `exp(arg - c*g.exp) * w^c`
5. `xreplace` the rewritten forms into `e`
6. Handle fractional exponents: lcm of all denominators → raise `w` to that power
7. Return `(f, logw)` where `logw = log(w)`

### 9. `std::pair<RCP<const Basic>, RCP<const Basic>> mrv_leadterm(const RCP<const Basic> &e, const RCP<const Symbol> &x)`

Returns `(c0, e0)` such that `e ∼ c0 * w^e0` as `w → 0`:

1. `mrv(e, x)` → `(Omega, exps)`
2. `moveup2/moveup` if `x` in Omega
3. `rewrite(exps, Omega, x, w)` → `(f, logw)`
4. Expand `f` as series in `w` via `series(f, w, prec)`
5. Extract leading term from `series->as_dict()`: find minimum exponent, return `(coeff, exponent)`
6. Substitute `log(w) → logw` in coefficient for consistency

### 10. `RCP<const Basic> limitinf(const RCP<const Basic> &e, const RCP<const Symbol> &x)`

Computes `lim_{x→∞} e(x)`.

**Logic:**
- Constant → return `e`
- `mrv_leadterm(e, x)` → `(c0, e0)`
- `sig = sign(e0, x)`:
  - `sig == 1` (e0 > 0): return `0`
  - `sig == -1` (e0 < 0): return `sign(c0) * oo`
  - `sig == 0`: recurse: `return limitinf(c0, x)`

### 11. `RCP<const Basic> gruntz(const RCP<const Basic> &e, const RCP<const Symbol> &z, const RCP<const Basic> &z0, const std::string &dir = "+")`

Main entry point. Converts any limit to `x → ∞`:

| `z0` | Dir | Substitution |
|------|-----|-------------|
| `oo` | any | None |
| `-oo` | any | `z → -z` |
| finite | `"+"` | `z → z0 + 1/z` |
| finite | `"-"` | `z → z0 - 1/z` |

Then calls `limitinf(e0, z)` and returns the result.

---

## Helper Functions

### `has(const Basic &e, const Basic &x)`
Wrapper around `has_symbol` for convenience.

### `independent(const RCP<const Basic> &e, const RCP<const Symbol> &x, ...)`
Splits expression into x-free constant and x-dependent expression. For `Add`, groups terms by `has_symbol`. For `Mul`, splits factors similarly.

### `powdenest_simple(const RCP<const Basic> &e)`
Normalizes nested powers: `(a^b)^c → a^(b*c)`.

### `tractable_rewrite(const RCP<const Basic> &e)`
Rewrites `sin`/`cos`/etc. to exponentials for the limit computation.

---

## Memoization

Many functions (`limitinf`, `mrv_leadterm`, `sign`) are called repeatedly with the same arguments. We'll use `std::unordered_map` with expression hashing for memoization to avoid exponential blowup.

---

## Key Architectural Challenges

| Challenge | Solution |
|-----------|----------|
| Mutual recursion (compare ↔ limitinf ↔ mrv) | Memoization breaks cycles; depth limits as safeguard |
| Series with `log(w)` terms | Post-process series to substitute `log(w)` → `logw` |
| `as_independent` equivalent | Custom helper using `has_symbol` |
| `powdenest` equivalent | Simple recursive helper using `pow()` |
| Tractable rewrite (trig → exp) | Inline conversion in `limitinf` |
| No Python bindings | Only C++ needed in this repo |

---

## Tests

Port key tests from `sympy/series/tests/test_gruntz.py`:

**Simple limits:**
- `lim_{x→∞} x → ∞`
- `lim_{x→∞} 1/x → 0`
- `lim_{x→∞} 1/x^2 → 0`
- `lim_{x→∞} (x+1)/x → 1`

**Exponential limits:**
- `lim_{x→∞} exp(x)/x → ∞`
- `lim_{x→∞} exp(-x)*x → 0`
- `lim_{x→∞} (exp(x)+1)/exp(x) → 1`

**Logarithmic limits:**
- `lim_{x→∞} log(x)/x → 0`
- `lim_{x→∞} log(log(x))/log(x) → 0`

**Power limits:**
- `lim_{x→∞} x^x / exp(x) → ∞`
- `lim_{x→∞} x^(1/x) → 1`

**MRV interaction:**
- `lim_{x→∞} exp(x) / exp(x + exp(-x)) → 1`

**Left/right limits:**
- `lim_{x→0+} 1/x → ∞`
- `lim_{x→0-} 1/x → -∞`
- `lim_{x→0+} log(x) → -∞`

---

## Build Integration

### `/work/symengine/CMakeLists.txt`

Add to `SRC` list (around line 92, after `series.cpp`):
```cmake
    gruntz.cpp
```

Add to `HEADERS` list (around line 233, after `series_visitor.h`):
```cmake
    gruntz.h
```

### `/work/symengine/tests/basic/CMakeLists.txt`

Add (around line 18, after `test_series`):
```cmake
add_executable(test_gruntz test_gruntz.cpp)
target_link_libraries(test_gruntz symengine catch)
add_test(test_gruntz ${PROJECT_BINARY_DIR}/test_gruntz)
```

---

## Implementation Phases

| Phase | Files | Est. Lines | Description |
|-------|-------|------------|-------------|
| 1 | `gruntz.h`, `gruntz.cpp` | 400 | `SubsSet`, helpers, `sign`, `compare` |
| 2 | `gruntz.cpp` | 250 | `mrv`, `mrv_max1`, `mrv_max3`, `moveup` |
| 3 | `gruntz.cpp` | 200 | `rewrite`, `mrv_leadterm`, series integration |
| 4 | `gruntz.cpp` | 150 | `limitinf`, `gruntz` (main entry) |
| 5 | `test_gruntz.cpp`, CMake | 100 | Tests and build integration |
