# Gruntz Algorithm Implementation for SymEngine — Detailed Plan (v2)

> **Target audience:** Junior C++ developer familiar with SymEngine basics.
> **Reference:** SymPy's `gruntz.py` — `/opt-3/cpython-v3.13-apt-deb/lib/python3.13/site-packages/sympy/series/gruntz.py` (702 lines).
> **Thesis:** Dominik Gruntz, *"On computing limits in a symbolic manipulation system"* (ETH Zurich, 1996).

---

## Table of Contents

1. [Algorithm Overview](#1-algorithm-overview)
2. [Files to Create/Modify](#2-files-to-create--modify)
3. [Header File: `gruntz.h`](#3-header-file-gruntzh)
4. [Implementation: `gruntz.cpp` — Functions in Detail](#4-implementation-gruntzcpp--functions-in-detail)
   - [4.1 SubsSet Class](#41-subsset-class)
   - [4.2 `has()` Helper](#42-has-helper)
   - [4.3 `independent()` Helper](#43-independent-helper)
   - [4.4 `powdenest_simple()` Helper](#44-powdenest_simple-helper)
   - [4.5 `tractable_rewrite()` Helper](#45-tractable_rewrite-helper)
   - [4.6 `sign()`](#46-sign)
   - [4.7 `compare()`](#47-compare)
   - [4.8 `mrv()` — Most Rapidly Varying Set](#48-mrv--most-rapidly-varying-set)
   - [4.9 `mrv_max3()`](#49-mrv_max3)
   - [4.10 `mrv_max1()`](#410-mrv_max1)
   - [4.11 `moveup2()` and `moveup()`](#411-moveup2-and-moveup)
   - [4.12 `build_expression_tree()`](#412-build_expression_tree)
   - [4.13 `rewrite()`](#413-rewrite)
   - [4.14 `mrv_leadterm()`](#414-mrv_leadterm)
   - [4.15 `limitinf()`](#415-limitinf)
   - [4.16 `gruntz()` — Main Entry Point](#416-gruntz--main-entry-point)
5. [Build Integration](#5-build-integration)
6. [Test Suite: `test_gruntz.cpp`](#6-test-suite-test_gruntzcpp)
7. [Memoization Strategy](#7-memoization-strategy)
8. [Debugging and Verification](#8-debugging-and-verification)
9. [Implementation Phases](#9-implementation-phases)
10. [Common Pitfalls and Checklist](#10-common-pitfalls-and-checklist)

---

## 1. Algorithm Overview

The Gruntz algorithm computes `lim_{z→z₀} f(z)` by:

```
Step 1: Transform any limit to x → +∞
        (z→oo: identity,  z→-oo: z→-z,  z→a: z→a ± 1/x)

Step 2: Find the MRV (Most Rapidly Varying) subexpressions of f(x)
        → returns Omega = SubsSet mapping: {expr → dummy}, plus rewritten f

Step 3: Rewrite f in terms of a new variable w → 0
        → every element of Omega gets rewritten as w^c, then xreplace

Step 4: Compute series expansion of the rewritten f in w
        → extract leading term: c₀ * w^e₀

Step 5: Read off the limit:
        e₀ > 0  →  limit = 0
        e₀ < 0  →  limit = sign(c₀) * ∞
        e₀ = 0  →  recurse: limit = lim c₀
```

**Comparability classes** (slowest → fastest):
```
2  ~  3  ~  -5        [constants]
x  ~  x²  ~  1/x      [polynomial-like]
exp(x)  ~  exp(-x)    [exponential-like]
exp(x²)                [each new exponent gives new class]
exp(exp(x))            [doubly exponential]
```

---

## 2. Files to Create / Modify

| File | Action | Purpose |
|------|--------|---------|
| `/work/symengine/gruntz.h` | **Create** | Class and function declarations |
| `/work/symengine/gruntz.cpp` | **Create** | Full implementation (~800-1000 lines) |
| `/work/symengine/tests/basic/test_gruntz.cpp` | **Create** | Test suite (~200-300 lines) |
| `/work/symengine/CMakeLists.txt` | **Modify** | Add `gruntz.cpp` to `SRC`, `gruntz.h` to `HEADERS` |
| `/work/symengine/tests/basic/CMakeLists.txt` | **Modify** | Add `test_gruntz` executable target |

---

## 3. Header File: `gruntz.h`

**File:** `/work/symengine/gruntz.h`

```cpp
#ifndef SYMENGINE_GRUNTZ_H
#define SYMENGINE_GRUNTZ_H

#include <symengine/basic.h>
#include <symengine/symbol.h>
#include <symengine/dict.h>

namespace SymEngine {

// =========================================================================
// SubsSet — mappings from mrv expressions to dummy variables
// =========================================================================

struct SubsSet {
    // Mapping: original mrv expression → dummy Symbol
    map_basic_basic dict;

    // Mapping: dummy Symbol → rewritten expression (in terms of other dummies)
    map_basic_basic rewrites;

    // If expr is not in dict, auto-create a Dummy and insert it.
    // Returns the dummy bound to expr.
    RCP<const Symbol> insert(const RCP<const Basic> &expr);

    // Returns true if self and other have at least one expression in common.
    bool meets(const SubsSet &other) const;

    // Merge self with other.  If the same expression appears in both, remap
    // to self's dummy and adjust exps.  Returns (merged_set, adjusted_exps).
    std::pair<SubsSet, RCP<const Basic>>
    union_set(const SubsSet &other, const RCP<const Basic> &exps) const;

    // Reverse substitution: replace every dummy with its original expression.
    RCP<const Basic> do_subs(const RCP<const Basic> &e) const;

    // True when dict is empty.
    bool empty() const;

    // Number of entries in dict.
    size_t size() const;
};

// =========================================================================
// Helper functions
// =========================================================================

// True if b contains symbol x (wrapper around has_symbol).
bool has(const Basic &b, const Basic &x);

// Split e into (x-independent part, x-dependent part).
// For Add: group terms that don't contain x into const_part.
// For Mul: group factors that don't contain x into const_part.
// Returns (const_part, dependent_part).
RCP<const Basic> independent(const RCP<const Basic> &e,
                             const RCP<const Symbol> &x,
                             Ptr<RCP<const Basic>> &dependent);

// Normalize nested powers: (a^b)^c → a^(b*c).
RCP<const Basic> powdenest_simple(const RCP<const Basic> &e);

// Rewrite sin/cos/tanh/... into exponentials for limit computation.
RCP<const Basic> tractable_rewrite(const RCP<const Basic> &e);

// =========================================================================
// Core algorithm functions
// =========================================================================

// Sign of e(x) as x → +∞.  Returns integer(1), integer(-1), or integer(0).
RCP<const Integer> sign(const RCP<const Basic> &e,
                         const RCP<const Symbol> &x);

// Compare growth rates a vs b as x → +∞.
// Returns '<' (a slower), '=' (same class), or '>' (a faster).
char compare(const RCP<const Basic> &a,
             const RCP<const Basic> &b,
             const RCP<const Symbol> &x);

// Find the Most Rapidly Varying subexpressions set.
// Returns (Omega, exps) where Omega is the SubsSet and exps is e rewritten.
std::pair<SubsSet, RCP<const Basic>>
mrv(const RCP<const Basic> &e, const RCP<const Symbol> &x);

// Combine two MRV sets when they may be in different classes.
std::pair<SubsSet, RCP<const Basic>>
mrv_max1(const SubsSet &f, const SubsSet &g,
         const RCP<const Basic> &exps, const RCP<const Symbol> &x);

// Combine two MRV sets already known to be in the same class.
std::pair<SubsSet, RCP<const Basic>>
mrv_max3(const SubsSet &f, const RCP<const Basic> &expsf,
         const SubsSet &g, const RCP<const Basic> &expsg,
         const SubsSet &_union, const RCP<const Basic> &expsboth,
         const RCP<const Symbol> &x);

// Replace x → exp(x) in all Omega entries (used when x itself is in Omega).
SubsSet moveup2(const SubsSet &s, const RCP<const Symbol> &x);
RCP<const Basic> moveup(const RCP<const Basic> &e,
                         const RCP<const Symbol> &x);

// Build dependency tree for Omega entries (helper for rewrite).
// Returns nodes map: Dummy -> NodeInfo (ht, before_list).
std::map<RCP<const Basic>, int, RCPBasicKeyLess>
build_expression_tree(const SubsSet &Omega);

// Rewrite e in terms of w (→ 0).
// Returns (rewritten_expression, logw) where logw = log(w).
std::pair<RCP<const Basic>, RCP<const Basic>>
rewrite(const RCP<const Basic> &e,
        const SubsSet &Omega,
        const RCP<const Symbol> &x,
        const RCP<const Symbol> &wsym);

// Leading term of series expansion: returns (c0, e0) such that
// e ∼ c0 * w^e0 as w → 0.
std::pair<RCP<const Basic>, RCP<const Basic>>
mrv_leadterm(const RCP<const Basic> &e,
             const RCP<const Symbol> &x);

// Limit at infinity: lim_{x→∞} e(x).
RCP<const Basic> limitinf(const RCP<const Basic> &e,
                           const RCP<const Symbol> &x);

// Main entry: limit of e(z) as z → z0 from direction dir ("+" or "-").
RCP<const Basic> gruntz(const RCP<const Basic> &e,
                         const RCP<const Symbol> &z,
                         const RCP<const Basic> &z0,
                         const std::string &dir = "+");

} // namespace SymEngine

#endif // SYMENGINE_GRUNTZ_H
```

**Includes needed in `gruntz.cpp`:**
```cpp
#include <symengine/gruntz.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/log.h>
#include <symengine/constants.h>
#include <symengine/functions.h>
#include <symengine/integer.h>
#include <symengine/rational.h>
#include <symengine/symbol.h>
#include <symengine/number.h>
#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/visitor.h>
#include <symengine/series.h>
#include <symengine/series_generic.h>
#include <symengine/symengine_exception.h>

#include <unordered_map>
#include <vector>
#include <algorithm>
#include <string>
#include <cmath>
```

---

## 4. Implementation: `gruntz.cpp` — Functions in Detail

### 4.1 `SubsSet` Class

**SymPy ref:** `gruntz.py` lines 156-245, class `SubsSet(dict)`.

#### Data Members

```cpp
struct SubsSet {
    // map<expr, dummy>  —  the key set is Omega (mrv set)
    map_basic_basic dict;       // std::map, ordered by RCPBasicKeyLess

    // map<dummy, rewritten_expr>  —  how dummies rewrite in terms of others
    map_basic_basic rewrites;   // e.g. {d3: exp(x - d2)}
};
```

**`map_basic_basic`** is `std::map<RCP<const Basic>, RCP<const Basic>, RCPBasicKeyLess>` defined in `/work/symengine/dict.h`. It uses `Basic::__cmp__` for ordering.

#### `RCP<const Symbol> insert(const RCP<const Basic> &expr)`

```cpp
RCP<const Symbol> SubsSet::insert(const RCP<const Basic> &expr)
{
    auto it = dict.find(expr);
    if (it != dict.end()) {
        return rcp_static_cast<const Symbol>(it->second);
    }
    RCP<const Symbol> d = dummy();   // or Dummy("d")
    dict[expr] = d;
    return d;
}
```

**Important:** In C++, unlike Python, we cannot overload `operator[]` to auto-insert. So every lookup that may insert must go through `insert()`. This means callers of `mrv()` etc. must use `insert()` instead of `dict[expr]` when they intend to create a new entry.

**Dummy creation:** SymEngine's `dummy()` from `symbol.h` creates a `Dummy` with auto-incremented index. Use raw `Dummy()` constructor in a loop to get predictable naming (optional).

#### `bool meets(const SubsSet &other) const`

Check if key sets intersect:

```cpp
bool SubsSet::meets(const SubsSet &other) const
{
    for (const auto &kv : dict) {
        if (other.dict.find(kv.first) != other.dict.end())
            return true;
    }
    return false;
}
```

#### `std::pair<SubsSet, RCP<const Basic>> union_set(...) const`

Equivalent of Python `SubsSet.union(s2, exps)` (line 224).

```
Algorithm:
  1. res = copy of self
  2. tr = {}  // remapping table: old_dummy → res_dummy
  3. For each (expr, var) in other:
       a. If expr already in self:
            - Tr[var] = res[expr]
            - If exps != null: exps = exps.xreplace({var: res[expr]})
       b. Else:
            - res[expr] = var
  4. For each (var, rewr) in other.rewrites:
       - res.rewrites[var] = rewr.xreplace(tr)
  5. Return (res, exps)
```

**Implementation notes:**
- Use `Basic::xreplace()` with a `map_basic_basic` for substitution.
- The `tr` table maps `RCP<const Basic>` → `RCP<const Basic>`.
- `xreplace` unlike `subs` does simple structural replacement without any mathematical simplification — exactly what we want.

#### `RCP<const Basic> do_subs(const RCP<const Basic> &e) const`

```cpp
RCP<const Basic> SubsSet::do_subs(const RCP<const Basic> &e) const
{
    map_basic_basic rev;  // dummy → original expr
    for (const auto &kv : dict) {
        rev[kv.second] = kv.first;
    }
    return e->xreplace(rev);
}
```

#### `empty()` and `size()`

```cpp
bool SubsSet::empty() const { return dict.empty(); }
size_t SubsSet::size() const { return dict.size(); }
```

---

### 4.2 `has()` Helper

```cpp
bool has(const Basic &b, const Basic &x)
{
    return has_symbol(b, x);  // from visitor.h
}
```

Uses `has_symbol` from `/work/symengine/visitor.h`. This is a visitor that walks the expression tree checking for the presence of a `Symbol` or `Dummy`.

---

### 4.3 `independent()` Helper

**Purpose:** Split `Add`/`Mul` into x-free constant and x-dependent expression. Equivalent to Python's `e.as_independent(x)`.

```cpp
// Returns (x_free_part, x_dependent_part)
// Uses output parameter Ptr for the dependent part.
RCP<const Basic> independent(const RCP<const Basic> &e,
                             const RCP<const Symbol> &x,
                             Ptr<RCP<const Basic>> &dependent)
```

**Algorithm:**

```
Case Add(e):
  - coef = e.get_coef()         // the numeric coefficient
  - dict = e.get_dict()         // {term: coefficient}
  - const_dict = {}             // terms without x
  - var_dict = {}               // terms with x
  - For each (term, coeff) in dict:
      if has_symbol(term, x): var_dict[term] = coeff
      else:                    const_dict[term] = coeff
  - If const_dict is empty: return zero, dependent = e
  - If var_dict is empty:   return e,    dependent = zero
  - const_part = Add::from_dict(coef + coeffs from const_dict, const_dict)
  - var_part   = Add::from_dict(const_var_dict...)
  - return const_part, dependent = var_part

Case Mul(e):
  - coef = e.get_coef()
  - dict = e.get_dict()          // {factor: exponent}
  - const_factors = {one: integer(1)}, const_coef = integer(1)
  - var_dict = {}
  - var_coef = integer(1)
  - For each (factor, exp) in dict:
      f = pow(factor, exp)       // power the factor by its exponent
      if has_symbol(f, x):
          var_dict[factor] = exp
      else:
          const_coef = mulnum(const_coef, pow(f, exp))
          // or use Mul::dict_add_term
  - If var_dict empty: return e, dependent = one
  - If const factors empty (only coef): return coef or one, dependent = var_part
  - Return Mul::from_dict(const_coef, const_factors), dependent = Mul::from_dict(var_coef, var_dict)

Default:
  - If has_symbol(e, x): return one, dependent = e
  - Else: return e, dependent = one
```

**Simplification:** For most practical cases in Gruntz, we can use a simpler version that extracts the constant factor:

```cpp
RCP<const Basic> constant_factor(const RCP<const Basic> &e,
                                  const RCP<const Symbol> &x)
{
    if (is_a<Mul>(*e)) {
        const Mul &m = down_cast<const Mul &>(*e);
        RCP<const Number> const_coef = m.get_coef();
        map_basic_basic const_factors;
        map_basic_basic var_factors;
        for (const auto &kv : m.get_dict()) {
            if (has_symbol(*kv.first, *x))
                var_factors[kv.first] = kv.second;
            else
                const_factors[kv.first] = kv.second;
        }
        if (const_factors.empty() && eq(*const_coef, *one))
            return one;
        return Mul::from_dict(const_coef, const_factors);
    } else if (is_a<Add>(*e)) {
        // similar: group terms by has_symbol
        // ...
    }
    if (!has_symbol(*e, *x))
        return e;
    return one;
}
```

**Important:** The Gruntz algorithm uses `as_independent` in `mrv()` when handling `Mul`/`Add` nodes. The constant factor gets multiplied back after recursing on the dependent part. If implementing the full version is too complex, a simpler alternative is to just return `(one, e)` and let the algorithm work without this optimization (though it may be slower).

---

### 4.4 `powdenest_simple()` Helper

**Purpose:** Normalize `(a^b)^c` → `a^(b*c)`. Equivalent to Python's `powdenest(expr, force=True)`.

```
RCP<const Basic> powdenest_simple(const RCP<const Basic> &e)
{
    // Walk expr tree bottom-up
    if is_a<Pow>(*e):
        const Pow &p = down_cast<const Pow &>(*e);
        RCP<const Basic> base = powdenest_simple(p.get_base());
        RCP<const Basic> expn = powdenest_simple(p.get_exp());
        if is_a<Pow>(*base):
            // (a^b)^c → a^(b*c)
            const Pow &inner = down_cast<const Pow &>(*base);
            return pow(inner.get_base(), mul(inner.get_exp(), expn));
        return pow(base, expn);
    if is_a<Mul>(*e):
        // recurse on factors
        const Mul &m = down_cast<const Mul &>(*e);
        map_basic_basic new_dict;
        for (const auto &kv : m.get_dict()) {
            new_dict[powdenest_simple(kv.first)] = kv.second;
        }
        return Mul::from_dict(m.get_coef(), new_dict);
    if is_a<Add>(*e):
        // recurse on terms
        const Add &a = down_cast<const Add &>(*e);
        umap_basic_num new_dict;
        for (const auto &kv : a.get_dict()) {
            new_dict[powdenest_simple(kv.first)] = kv.second;
        }
        return Add::from_dict(a.get_coef(), new_dict);
    // Functions: recurse on args
    ...
    return e;
}
```

**Alternatively:** Use `expand()` — in many cases `expand()` already flattens nested powers.

---

### 4.5 `tractable_rewrite()` Helper

**Purpose:** Rewrite trig/hyperbolic functions to exponential form so the limit algorithm can handle them.

```cpp
RCP<const Basic> tractable_rewrite(const RCP<const Basic> &e)
```

**Transformations needed:**

| Function | exp Rewrite |
|----------|-------------|
| `sin(x)` | `(exp(I*x) - exp(-I*x)) / (2*I)` |
| `cos(x)` | `(exp(I*x) + exp(-I*x)) / 2` |
| `tan(x)` | `sin(x)/cos(x)` → exp form |
| `sinh(x)` | `(exp(x) - exp(-x)) / 2` |
| `cosh(x)` | `(exp(x) + exp(-x)) / 2` |
| `tanh(x)` | `(exp(x) - exp(-x)) / (exp(x) + exp(-x))` |
| `coth(x)` | `(exp(x) + exp(-x)) / (exp(x) - exp(-x))` |
| `log(x)` | keep as-is (already tractable) |
| `exp(x)` | keep as-is (already tractable) |

**Implementation:** Walk the expression tree bottom-up. When encountering a function with a known rewrite, replace it. Use `Basic::xreplace()` or recursion.

```cpp
RCP<const Basic> tractable_rewrite(const RCP<const Basic> &e)
{
    // Use rewrite_as_exp() from basic.h for trig functions
    // SymEngine already has rewrite_as_exp()!
    RCP<const Basic> result = e;
    return rewrite_as_exp(result);  // Handles sin, cos, tan, cot, etc.
}
```

**Key insight:** SymEngine already has `rewrite_as_exp()` in `/work/symengine/basic.h`! It handles `Sin`, `Cos`, `Tan`, `Cot`, `Csc`, `Sec`, `Sinh`, `Cosh`, `Tanh`, `Coth`, `Csch`, `Sech`. Just use it.

For hyperbolic functions not covered by `rewrite_as_exp()`, implement manually:
```cpp
if (is_a<Sinh>(*e)) {
    auto arg = down_cast<const Sinh &>(*e).get_arg();
    return div(sub(exp(arg), exp(neg(arg))), integer(2));
}
// similar for Cosh, Tanh, etc.
```

**When to apply:** In `limitinf()` before calling `mrv_leadterm()`. In `mrv()` for `Log` arguments that may contain trig.

---

### 4.6 `sign()` — Determine Sign as x → ∞

**SymPy ref:** `gruntz.py` lines 368-422.

```cpp
// Returns integer(1), integer(-1), or integer(0)
RCP<const Integer> sign(const RCP<const Basic> &e,
                         const RCP<const Symbol> &x)
```

**Algorithm:**

```
function sign(e, x):
    # --- Fast path: use numeric properties ---
    if is_a_Number(e):
        if e->is_positive(): return integer(1)
        if e->is_negative(): return integer(-1)
        if e->is_zero():     return integer(0)
        // fall through

    # --- No x dependency → constant → evaluate directly ---
    if not has_symbol(e, x):
        // For pure numbers, use is_positive/is_negative
        // For expressions without x, try to evaluate numerically
        // (if advanced numeric eval is not available, fall through)

    # --- x itself → +1 ---
    if eq(e, x):
        return integer(1)

    # --- exp(...) → always +1 ---
    if is_exp_like(e):
        return integer(1)

    # --- Mul: multiply signs of factors (split into two) ---
    if is_a<Mul>(*e):
        RCP<const Basic> a, b;
        down_cast<const Mul &>(*e).as_two_terms(outArg(a), outArg(b));
        int sa = sign(a, x);
        if sa == 0: return integer(0);
        return mul(sa, sign(b, x));  // sa * sb

    # --- Pow(E, ...) → +1 ---
    if is_a<Pow>(*e):
        const Pow &p = down_cast<const Pow &>(*e);
        if eq(*p.get_base(), *E):
            return integer(1)
        // For base^exp where base is positive → +1
        int sb = sign(p.get_base(), x);
        if sb == 1:
            return integer(1)
        // For even integer exponent: result is positive
        if is_a<Integer>(*p.get_exp()):
            int exp_val = down_cast<const Integer &>(*p.get_exp()).as_int();
            if exp_val % 2 == 0:
                return integer(1)  // even power → positive
    }
    # --- Log(pos_arg) → sign(arg - 1) ---
    if is_a<Log>(*e):
        RCP<const Basic> arg = down_cast<const Log &>(*e).get_arg();
        // Only valid if arg is known positive
        RCP<const Basic> diff = sub(arg, one);
        return sign(diff, x);
    }

    # --- Fallback: use mrv_leadterm ---
    RCP<const Basic> c0, e0;
    std::tie(c0, e0) = mrv_leadterm(e, x);
    return sign(c0, x);
```

**Key points:**
- Return type is `RCP<const Integer>` — use `integer(1)`, `integer(-1)`, `integer(0)`.
- For `Mul`, split into two factors with `as_two_terms()`, then multiply signs.
- `is_exp_like(e)` = `is_a<Pow>(*e) && eq(*p.get_base(), *E)`.
- The fallback calls `mrv_leadterm()` which itself calls `sign()` — this is mutual recursion. Memoization is essential to avoid infinite loops.

---

### 4.7 `compare()` — Compare Growth Rates

**SymPy ref:** `gruntz.py` lines 138-153.

```cpp
char compare(const RCP<const Basic> &a,
             const RCP<const Basic> &b,
             const RCP<const Symbol> &x)
```

**Algorithm:**

```
function compare(a, b, x):
    # Compute log(a) and log(b)
    la = log(a)
    lb = log(b)

    # Special case: if a = exp(p), use p directly (log(exp(p)) → p)
    if is_exp_like(a):
        la = exponent_of(a)    // for exp(p): la = p
    if is_exp_like(b):
        lb = exponent_of(b)    // for exp(p): lb = p

    # Compute L = lim_{x→∞} la / lb
    L = limitinf(div(la, lb), x)

    # 0 → slower, infinite → faster, otherwise → equal
    if eq(*L, *zero):
        return '<'
    if is_a<Infty>(*L):
        return '>'
    return '='
```

**Important:** `log(exp(p))` should simplify to `p`. Use `is_exp_like()` helper:
```cpp
static bool is_exp_like(const RCP<const Basic> &e)
{
    if (!is_a<Pow>(*e))
        return false;
    const Pow &p = down_cast<const Pow &>(*e);
    return eq(*p.get_base(), *E);
}

static RCP<const Basic> exponent_of(const RCP<const Basic> &e)
{
    // e is exp(something) = E^something
    return down_cast<const Pow &>(*e).get_exp();
}
```

**Edge cases:**
- `compare(a, b, x)` and `compare(b, a, x)` should give opposite results.
- If `la/lb` simplifies to a constant, `limitinf` returns it, and compare returns `'='`.
- If `a` or `b` does NOT contain `x`, it's in the constant class (slower than any polynomial+).

---

### 4.8 `mrv()` — Most Rapidly Varying Set

**SymPy ref:** `gruntz.py` lines 249-321.

This is the **core recursive decomposition** function.

```cpp
std::pair<SubsSet, RCP<const Basic>>
mrv(const RCP<const Basic> &e, const RCP<const Symbol> &x)
```

**Full algorithm (type-dispatch):**

```
function mrv(e, x):

    # --- BASELINE CASES ---

    # Case 1: e does not depend on x
    if not has_symbol(e, x):
        return (SubsSet(), e)

    # Case 2: e == x
    if eq(e, x):
        s = SubsSet()
        d = s.insert(x)         // x → dummy
        return (s, d)

    # --- COMPOSITE CASES ---

    # Case 3: e is Mul or Add
    if is_a<Mul>(e) or is_a<Add>(e):
        # Split into x-independent and x-dependent parts
        const_part = independent(e, x, outArg(dep_part))
        if const_part is not e:
            # Some parts are independent of x
            (s, expr) = mrv(dep_part, x)
            return (s, e.func(const_part, expr))  // reconstruct
        # All parts depend on x: split into two terms
        a, b = e.as_two_terms()
        (s1, e1) = mrv(a, x)
        (s2, e2) = mrv(b, x)
        return mrv_max1(s1, s2, reconstruct(e, a→e1, b→e2), x)

    # Case 4: e is Pow (base != E)
    if is_a<Pow>(e) and not is_exp_like(e):
        # Rewrite base^exp → exp(exp * log(base))
        e1 = integer(1)
        current = e
        while is_a<Pow>(current):
            base = current.get_base()
            e1 = mul(e1, current.get_exp())
            current = base
        if eq(*current, *one):
            return (SubsSet(), current)
        if has_symbol(e1, x):
            return mrv(exp(mul(e1, log(current))), x)
        else:
            (s, expr) = mrv(current, x)
            return (s, pow(expr, e1))

    # Case 5: e is Log
    if is_a<Log>(e):
        arg = e.get_arg()
        (s, expr) = mrv(arg, x)
        return (s, log(expr))

    # Case 6: e is exp (or Pow(E, ...))
    if is_exp_like(e):
        exp_arg = exponent_of(e)
        # Simplify: exp(log(p)) → p
        if is_a<Log>(exp_arg):
            return mrv(exp_arg.get_arg(), x)

        # Check if exp → ∞
        li = limitinf(exp_arg, x)
        if is_infinite(li):  # exp → ∞ or exp → 0
            s1 = SubsSet()
            e1 = s1.insert(e)         // add e to Omega
            (s2, e2) = mrv(exp_arg, x)
            su = s1.union_set(s2).first
            su.rewrites[e1] = exp(e2)  // e = exp(rewritten_arg)
            return mrv_max3(s1, e1, s2, exp(e2), su, e1, x)
        else:
            (s, expr) = mrv(exp_arg, x)
            return (s, exp(expr))

    # Case 7: Other single-argument functions
    if is single-arg function:
        (s, expr) = mrv(e.get_arg(), x)
        return (s, e.create(expr))  // or use factory

    # Case 8: Unhandled
    throw NotImplementedError("mrv: unhandled type")
```

**Implementation details:**

For **Case 3 (Mul/Add)**:
- `as_two_terms()` is available on both `Mul` and `Add`.
- `reconstruct()`: If `e = Mul(coef, {a: 2, b: 3})` and we changed `a → e1`, `b → e2`, rebuild as `Mul(coef, {e1: 2, e2: 3})`. Similarly for Add.
- For simplicity, the reconstruction can be done by creating a new expression from the two terms: for Add: `add(e1, e2)`, for Mul: `mul(e1, e2)` (ignoring the full factor structure), then multiplying back the constant part.

For **Case 4 (Pow)**:
- The loop `while is_a<Pow>(current)` collects nested powers into `e1 = exponent_product`.
- Example: `(x^2)^(x+1)` → `e1 = x+1`, `current = x^2` → then `e1 = (x+1)*2 = 2x+2`, `current = x`.
- After the loop, if base has x, rewrite as `exp(e1 * log(base))` and recurse.

For **Case 6 (exp)**:
- `is_infinite(li)`: check if `li` is an `Infty` object.
- The key logic: if the exponent goes to infinity, then the exp itself goes to either 0 or ∞, making it a candidate for the MRV set. We add it to Omega and recurse on the exponent.
- `su.rewrites[e1] = exp(e2)`: records that `e1` (in terms of x) rewrites to `exp(e2)` (where `e2` is the rewritten exponent).

---

### 4.9 `mrv_max3()` — Combine Two MRV Sets (Same Class)

**SymPy ref:** `gruntz.py` lines 324-350.

```cpp
std::pair<SubsSet, RCP<const Basic>>
mrv_max3(const SubsSet &f, const RCP<const Basic> &expsf,
         const SubsSet &g, const RCP<const Basic> &expsg,
         const SubsSet &_union, const RCP<const Basic> &expsboth,
         const RCP<const Symbol> &x)
```

**Algorithm:**

```
function mrv_max3(f, expsf, g, expsg, _union, expsboth, x):
    # Trivial cases: one set is empty
    if f.empty():    return (g, expsg)
    if g.empty():    return (f, expsf)

    # Sets intersect (share at least one expression) → use union
    if f.meets(g):   return (_union, expsboth)

    # Compare first elements
    c = compare(first_key(f), first_key(g), x)
    if c == '>':     return (f, expsf)     # f dominates
    if c == '<':     return (g, expsg)     # g dominates
    // c == '=':     return (_union, expsboth)  # same class → merge
```

**`first_key()`:** Extract the first expression from `f.dict`. Since `dict` is a `std::map` ordered by `RCPBasicKeyLess`, just use `f.dict.begin()->first`.

---

### 4.10 `mrv_max1()` — Combine Two MRV Sets

**SymPy ref:** `gruntz.py` lines 353-362.

```cpp
std::pair<SubsSet, RCP<const Basic>>
mrv_max1(const SubsSet &f, const SubsSet &g,
         const RCP<const Basic> &exps, const RCP<const Symbol> &x)
```

**Algorithm:**

```
function mrv_max1(f, g, exps, x):
    # Compute union of f and g, adjusting exps accordingly
    (u, b) = f.union_set(g, exps)

    # Compare both orderings to determine which dominates
    return mrv_max3(f, g.do_subs(exps),
                    g, f.do_subs(exps),
                    u, b, x)
```

**Key insight:** We pass `exps` twice — once with f's dummies replaced by g's, and once with g's replaced by f's. This ensures that the comparison uses expressions from the same "namespace" of dummies.

---

### 4.11 `moveup2()` and `moveup()`

**SymPy ref:** `gruntz.py` lines 473-483.

These are used when `x` itself appears in the MRV set Omega. We "move up" to the exponential level by replacing `x → exp(x)` everywhere.

```cpp
SubsSet moveup2(const SubsSet &s, const RCP<const Symbol> &x)
{
    SubsSet r;
    map_basic_basic xmap;
    xmap[x] = exp(x);

    for (const auto &kv : s.dict) {
        r.dict[kv.first->xreplace(xmap)] = kv.second;
    }
    for (const auto &kv : s.rewrites) {
        r.rewrites[kv.first] = kv.second->xreplace(xmap);
    }
    return r;
}
```

```cpp
RCP<const Basic> moveup(const RCP<const Basic> &e,
                         const RCP<const Symbol> &x)
{
    map_basic_basic xmap;
    xmap[x] = exp(x);
    return e->xreplace(xmap);
}
```

---

### 4.12 `build_expression_tree()` — Topological Sort for Rewrite

**SymPy ref:** `gruntz.py` lines 541-578.

**Purpose:** Determine the order in which Omega entries should be rewritten. Dependencies are tracked: if dummy `d3`'s rewrite expression contains `d2`, then `d2` must be rewritten before `d3`.

**Algorithm:**

```
function build_expression_tree(Omega):
    nodes = {}  // dummy_var → index (height)

    // Initialize all nodes with height 1
    for each (expr, var) in Omega:
        nodes[var] = 1

    // Compute heights: if var's rewrite depends on var2, increase var's height
    changed = true
    while changed:
        changed = false
        for each (expr, var) in Omega:
            if var in Omega.rewrites:
                r = Omega.rewrites[var]
                for each (_, var2) in Omega:
                    if has_symbol(r, var2):
                        new_h = nodes[var2] + 1
                        if new_h > nodes[var]:
                            nodes[var] = new_h
                            changed = true

    return nodes  // map: dummy → height
```

**Simplified approach:**

Instead of building explicit tree nodes, compute a height map. Height = 1 + max(height of all dependencies). This gives a topological ordering (higher height = depends on more things = should be processed first / later).

The Omega items are then sorted by height (descending) so that the most-dependent expressions are processed first.

```cpp
std::map<RCP<const Basic>, int, RCPBasicKeyLess>
build_expression_tree(const SubsSet &Omega)
{
    std::map<RCP<const Basic>, int, RCPBasicKeyLess> nodes;

    for (const auto &kv : Omega.dict) {
        nodes[kv.second] = 1;  // dummy → height
    }

    bool changed = true;
    while (changed) {
        changed = false;
        for (const auto &kv : Omega.dict) {
            auto rit = Omega.rewrites.find(kv.second);
            if (rit != Omega.rewrites.end()) {
                const auto &r = rit->second;
                for (const auto &kv2 : Omega.dict) {
                    if (has_symbol(*r, *kv2.second)) {
                        int new_h = nodes[kv2.second] + 1;
                        if (new_h > nodes[kv.second]) {
                            nodes[kv.second] = new_h;
                            changed = true;
                        }
                    }
                }
            }
        }
    }
    return nodes;
}
```

---

### 4.13 `rewrite()` — The Hardest Function

**SymPy ref:** `gruntz.py` lines 583-659.

**Purpose:** Given Omega (the MRV set) and expression `e`, rewrite everything in terms of a new variable `w` where `w → 0` as `x → ∞`.

```cpp
std::pair<RCP<const Basic>, RCP<const Basic>>
rewrite(const RCP<const Basic> &e,
        const SubsSet &Omega,
        const RCP<const Symbol> &x,
        const RCP<const Symbol> &wsym)
```

**Full algorithm:**

```
function rewrite(e, Omega, x, wsym):
    # --- Validation ---
    assert not Omega.empty()
    assert all items in Omega are exp-like

    rewrites_map = Omega.rewrites
    Omega_items = Omega.dict items as a vector of (expr, var) pairs

    # --- Step 1: Sort Omega by dependency height ---
    nodes = build_expression_tree(Omega)
    sort Omega_items by nodes[var] descending

    # --- Step 2: Determine sign of each g.exp and pick "w" ---
    # The last g (highest height, simplest) is the representative.
    g = last Omega_item
    sig = sign(g.exponent, x)    // sign of g's exponent as x → ∞

    if sig == 1:
        wsym = pow(wsym, integer(-1))   // wsym = 1/wsym (we need w → 0)

    # --- Step 3: Compute O2 (rewriting rules) ---
    O2 = []          // list of (var, rewritten_expr) pairs
    denominators = []

    for each (f, var) in Omega_items:
        f_exp = f.exponent    // f = exp(f_exp)

        # Compute c = lim f_exp / g_exp
        c = limitinf(div(f_exp, g_exp), x)

        # Track denominators for LCM later
        if is_a<Rational>(c):
            denominators.push(c.denominator)

        # Determine arg to use in rewriting
        arg = f_exp
        if var in rewrites_map:
            rw = rewrites_map[var]
            if is_exp_like(rw):
                arg = rw.exponent
            else:
                raise ValueError("rewrites value should be exp")

        # Rewritten form: exp(arg - c * g_exp) * wsym^c
        rewritten = mul(exp(sub(arg, mul(c, g_exp))),
                        pow(wsym, c))
        O2.push((var, rewritten))

    # --- Step 4: Substitute O2 into e ---
    f = powdenest_simple(e)   // or powsimp(e, combine='exp')
    for each (var, rewritten) in O2:
        f = f.xreplace({var: rewritten})

    # Verify no Omega variables remain in f
    for each (_, var) in Omega_items:
        assert not has_symbol(f, var)
        // NOTE: In C++, using SYMENGINE_ASSERT

    # --- Step 5: Compute logw ---
    logw = g_exp           // = exponent of g
    if sig == 1:
        logw = neg(logw)   // log(w) → log(1/w) = -log(w)

    # --- Step 6: Handle fractional exponents (LCM normalization) ---
    exponent = lcm(all denominators, default=1)
    f = f.subs(wsym, pow(wsym, exponent))
    logw = div(logw, integer(exponent))

    # --- Step 7: Simplify ---
    # In Python: f = bottom_up(f, lambda w: getattr(w, 'normal', lambda: w)())
    # In C++: apply basic simplification (expand, cancel, etc.)
    f = expand(f)

    return (f, logw)
```

**Key implementation details:**

1. **`is_exp_like(f)`:** Check `is_a<Pow>(*f) && eq(*p.get_base(), *E)`.

2. **`exponent_of(expr)`:** For `exp(arg)` = `pow(E, arg)`, return `arg`. For a non-exp expression, we need to decide — in `rewrite()`, all Omega items MUST be exp-like (checked in validation).

3. **LCM of denominators:** Need a helper to compute LCM of integers:
   ```cpp
   integer_class lcm_of(const std::vector<integer_class> &nums) {
       if (nums.empty()) return integer_class(1);
       integer_class result = nums[0];
       for (size_t i = 1; i < nums.size(); i++)
           result = mp_lcm(result, nums[i]);  // from mp_gcd
       return result;
   }
   ```

4. **`xreplace` vs `subs`:** Always use `xreplace` for structural replacement — it does NOT try to mathematically simplify. This is crucial because we want exact replacement, not simplification.

5. **The `bottom_up` final pass:** In SymPy, after substitution, a `bottom_up` traversal applies `.normal()` to each subexpression. In C++, we can apply `expand()` for a similar effect, or skip it for the initial implementation and handle edge cases as they arise.

6. **Sorting Omega_items by height:** After calling `build_expression_tree`, create a vector of `(expr, var)` pairs from `Omega.dict`, then sort by `nodes[var]` descending.

---

### 4.14 `mrv_leadterm()` — Leading Term of Series

**SymPy ref:** `gruntz.py` lines 489-538.

```cpp
std::pair<RCP<const Basic>, RCP<const Basic>>
mrv_leadterm(const RCP<const Basic> &e,
             const RCP<const Symbol> &x)
```

**Algorithm:**

```
function mrv_leadterm(e, x):
    # --- Step 1: Compute Omega via mrv ---
    if not has_symbol(e, x):
        return (e, integer(0))

    Omega, exps = mrv(e, x)

    if Omega.empty():
        return (exps, integer(0))

    # --- Step 2: If x is in Omega, move up ---
    if Omega.dict contains x:
        Omega = moveup2(Omega, x)
        exps = moveup(exps, x)

    # --- Step 3: Rewrite in terms of w (positive dummy) ---
    w = dummy("w", positive=true)   // or just dummy("w")
    f, logw = rewrite(exps, Omega, x, w)

    # --- Step 4: Compute series expansion in w ---
    prec = 6   // starting precision; increase if needed

    series_result = series(f, w, prec)
    series_dict = series_result->as_dict()

    # Find the lowest exponent (leading term)
    # series_dict is map<int, RCP<const Basic>> sorted by key
    lead_exp = series_dict.begin()->first     // lowest exponent
    lead_coeff = series_dict.begin()->second  // coefficient

    # --- Step 5: Substitute log(w) → logw ---
    # In SymPy: lt[0].subs(log(w), logw)
    map_basic_basic logsub;
    logsub[log(w)] = logw;
    lead_coeff = lead_coeff->xreplace(logsub);

    return (lead_coeff, integer(lead_exp));
```

**Important details:**

1. **Series precision (`prec`):** Start with 6. If `lead_coeff` is 0 for the lowest exponent (due to cancellation), increase `prec` and retry. Actually, the series `as_dict()` returns all non-zero coefficients, so the first entry is guaranteed to be non-zero. However, if cancellation occurs, the true leading term might need higher precision to discover — increase `prec` if the lead coefficient still depends on `w` (indicating incomplete cancellation).

2. **Dealing with log(w) in coefficients:** The series expansion produces terms like `log(w) * w^k`. After extracting the leading coefficient, we substitute `log(w) → logw` so the result is in terms of `logw` (which equals `log(g_exp)`).

3. **`series()` function:** Uses `SymEngine::series(ex, var, prec)` from `series.h`. This returns `RCP<const SeriesCoeffInterface>`. Call `as_dict()` to get the `umap_int_basic` (which is `std::map<int, RCP<const Basic>>` sorted by exponent).

4. **Positive dummy `w`:** In SymPy, `w = Dummy("w", positive=True)` — this is important because `log(w*2)` should expand to `log(w) + log(2)`. In SymEngine, `Dummy` doesn't carry assumptions, so we rely on the series code handling it correctly.

5. **`moveup` condition:** Check if `Omega.dict` contains `x` as a key. In SymEngine: iterate through Omega.dict and check if any key equals `x`.

**Fallback for series failure:**

If `series()` throws `NotImplementedError`, fall back to trying higher precision:
```cpp
for (int prec = 6; prec <= 24; prec *= 2) {
    try {
        auto ser = series(f, w, prec);
        // extract leading term
        if (lead_coeff_found) return ...;
    } catch (const NotImplementedError &) {
        // try higher precision
    }
}
throw NotImplementedError("Cannot compute series for leadterm");
```

---

### 4.15 `limitinf()` — Limit at Infinity

**SymPy ref:** `gruntz.py` lines 428-470.

```cpp
RCP<const Basic> limitinf(const RCP<const Basic> &e,
                           const RCP<const Symbol> &x)
```

**Algorithm:**

```
function limitinf(e, x):
    # --- Skip if constant ---
    if not has_symbol(e, x):
        return e

    # --- Apply tractable rewrite ---
    e = tractable_rewrite(e)

    # --- Normalize powers ---
    e = powdenest_simple(e)

    # --- Compute leading term ---
    c0, e0 = mrv_leadterm(e, x)

    # --- Determine limit from leading term ---
    sig = sign(e0, x)

    if sig == integer(1):
        return integer(0)       // e0 > 0 → limit = 0

    if sig == integer(-1):
        // e0 < 0 → limit = ±∞
        sc = sign(c0, x)
        if eq(*sc, *zero):
            throw std::runtime_error("Leading term sign should not be zero")
        if eq(*sc, *integer(1)):
            return Inf           // +∞
        else:
            return NegInf        // -∞

    if sig == integer(0):
        // e0 == 0 → recurse on c0
        return limitinf(c0, x)

    throw std::runtime_error("Unexpected sign result");
```

**Edge cases:**
- If `c0` is complex (contains `I`), the sign function needs special handling.
- The recursion on `c0` when `e0 == 0` is guaranteed to terminate because each successive `limitinf` call operates on an expression from a **lower comparability class**.

---

### 4.16 `gruntz()` — Main Entry Point

**SymPy ref:** `gruntz.py` lines 662-702.

```cpp
RCP<const Basic> gruntz(const RCP<const Basic> &e,
                         const RCP<const Symbol> &z,
                         const RCP<const Basic> &z0,
                         const std::string &dir)
```

**Algorithm:**

```
function gruntz(e, z, z0, dir="+"):
    # --- Validate: z must be a Symbol ---
    if not is_a<Symbol>(z):
        throw NotImplementedError("Second argument must be a Symbol")

    # --- Transform to x → +∞ ---
    if is_a<Infty>(z0):
        inf = down_cast<const Infty>(z0)
        if inf->is_positive_infinity():
            e0 = e
        elif inf->is_negative_infinity():
            e0 = e.subs(z, neg(z))
        else:
            throw NotImplementedError("Unsigned infinity not supported")
    else:
        # Finite z0: z → z0 ± 1/z
        if dir == "+":
            e0 = e.subs(z, add(z0, div(one, z)))   // z0 + 1/z
        elif dir == "-":
            e0 = e.subs(z, sub(z0, div(one, z)))   // z0 - 1/z
        else:
            throw NotImplementedError("dir must be '+' or '-'")

    # --- Compute limit ---
    return limitinf(e0, z)
```

**Key: `subs()` API**

SymEngine's `Basic::subs()` takes a `map_basic_basic`:
```cpp
map_basic_basic d;
d[z] = replacement;
RCP<const Basic> result = e->subs(d);
```

---

## 5. Build Integration

### 5.1 `/work/symengine/CMakeLists.txt`

**SRC list** — insert `gruntz.cpp` between `functions.cpp` and `infinity.cpp`:

Find the block (approximately lines 54-56):
```cmake
    functions.cpp
    infinity.cpp
    integer.cpp
```

Change to:
```cmake
    functions.cpp
    gruntz.cpp
    infinity.cpp
    integer.cpp
```

**HEADERS list** — insert `gruntz.h` between `functions.h` and `infinity.h`:

Find the block (approximately lines 184-186):
```cmake
    functions.h
    infinity.h
    integer.h
```

Change to:
```cmake
    functions.h
    gruntz.h
    infinity.h
    integer.h
```

### 5.2 `/work/symengine/tests/basic/CMakeLists.txt`

Append at the end of the file (after the last `add_test` block, approximately line 142):

```cmake
add_executable(test_gruntz test_gruntz.cpp)
target_link_libraries(test_gruntz symengine catch)
add_test(test_gruntz ${PROJECT_BINARY_DIR}/test_gruntz)
```

---

## 6. Test Suite: `test_gruntz.cpp`

**File:** `/work/symengine/tests/basic/test_gruntz.cpp`

Uses the **Catch2** test framework (included via `#include "catch.hpp"`).

### Test Structure

```cpp
#include "catch.hpp"
#include <symengine/basic.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/constants.h>
#include <symengine/functions.h>
#include <symengine/integer.h>
#include <symengine/symbol.h>
#include <symengine/series.h>
#include <symengine/gruntz.h>

using SymEngine::RCP;
using SymEngine::Basic;
using SymEngine::Symbol;
using SymEngine::Integer;
using SymEngine::Number;
using SymEngine::Add;
using SymEngine::Mul;
using SymEngine::Pow;
using SymEngine::Log;
using SymEngine::Infty;
using SymEngine::Dummy;
using SymEngine::SubsSet;
using SymEngine::add;
using SymEngine::mul;
using SymEngine::pow;
using SymEngine::div;
using SymEngine::sub;
using SymEngine::exp;
using SymEngine::log;
using SymEngine::neg;
using SymEngine::symbol;
using SymEngine::integer;
using SymEngine::dummy;
using SymEngine::one;
using SymEngine::zero;
using SymEngine::minus_one;
using SymEngine::E;
using SymEngine::Inf;
using SymEngine::NegInf;
using SymEngine::is_a;
using SymEngine::eq;
using SymEngine::down_cast;
using SymEngine::make_rcp;
```

### Test Categories

#### 6.1 Basic Limits (x → ∞)

```cpp
TEST_CASE("Basic limits x->oo", "[gruntz]")
{
    auto x = symbol("x");

    // x → ∞
    auto r = gruntz(x, x, Inf);
    REQUIRE(eq(*r, *Inf));

    // 1/x → 0
    r = gruntz(div(integer(1), x), x, Inf);
    REQUIRE(eq(*r, *integer(0)));

    // 1/x^2 → 0
    r = gruntz(div(integer(1), pow(x, integer(2))), x, Inf);
    REQUIRE(eq(*r, *integer(0)));

    // (x+1)/x → 1
    r = gruntz(div(add(x, integer(1)), x), x, Inf);
    REQUIRE(eq(*r, *integer(1)));
}
```

#### 6.2 Exponential Limits

```cpp
TEST_CASE("Exponential limits", "[gruntz]")
{
    auto x = symbol("x");

    // exp(x)/x → ∞
    auto r = gruntz(div(exp(x), x), x, Inf);
    REQUIRE(eq(*r, *Inf));

    // exp(-x)*x → 0
    r = gruntz(mul(exp(neg(x)), x), x, Inf);
    REQUIRE(eq(*r, *integer(0)));

    // (exp(x)+1)/exp(x) → 1
    r = gruntz(div(add(exp(x), integer(1)), exp(x)), x, Inf);
    REQUIRE(eq(*r, *integer(1)));
}
```

#### 6.3 Logarithmic Limits

```cpp
TEST_CASE("Logarithmic limits", "[gruntz]")
{
    auto x = symbol("x");

    // log(x)/x → 0
    auto r = gruntz(div(log(x), x), x, Inf);
    REQUIRE(eq(*r, *integer(0)));

    // log(log(x))/log(x) → 0
    r = gruntz(div(log(log(x)), log(x)), x, Inf);
    REQUIRE(eq(*r, *integer(0)));
}
```

#### 6.4 Power Limits

```cpp
TEST_CASE("Power limits", "[gruntz]")
{
    auto x = symbol("x");

    // x^x / exp(x) → ∞
    auto r = gruntz(div(pow(x, x), exp(x)), x, Inf);
    REQUIRE(eq(*r, *Inf));

    // x^(1/x) → 1
    r = gruntz(pow(x, div(integer(1), x)), x, Inf);
    REQUIRE(eq(*r, *integer(1)));
}
```

#### 6.5 MRV Interaction

```cpp
TEST_CASE("MRV interaction", "[gruntz]")
{
    auto x = symbol("x");

    // exp(x) / exp(x + exp(-x)) → 1
    auto r = gruntz(div(exp(x), exp(add(x, exp(neg(x))))), x, Inf);
    REQUIRE(eq(*r, *integer(1)));
}
```

#### 6.6 Left/Right (Finite) Limits

```cpp
TEST_CASE("Finite directional limits", "[gruntz]")
{
    auto x = symbol("x");

    // lim_{x→0+} 1/x → ∞
    auto r = gruntz(div(integer(1), x), x, integer(0), "+");
    REQUIRE(eq(*r, *Inf));

    // lim_{x→0-} 1/x → -∞
    r = gruntz(div(integer(1), x), x, integer(0), "-");
    REQUIRE(eq(*r, *NegInf));

    // lim_{x→0+} log(x) → -∞
    r = gruntz(log(x), x, integer(0), "+");
    REQUIRE(eq(*r, *NegInf));
}
```

#### 6.7 Negative Infinity

```cpp
TEST_CASE("Negative infinity limits", "[gruntz]")
{
    auto x = symbol("x");

    // lim_{x→-∞} exp(x) → 0
    auto r = gruntz(exp(x), x, NegInf);
    REQUIRE(eq(*r, *integer(0)));
}
```

#### 6.8 Unit Tests for Internal Functions

```cpp
// --- Sign ---
TEST_CASE("sign() function", "[gruntz_internal]")
{
    auto x = symbol("x");

    REQUIRE(eq(*sign(integer(5), x), *integer(1)));
    REQUIRE(eq(*sign(integer(-3), x), *integer(-1)));
    REQUIRE(eq(*sign(integer(0), x), *integer(0)));
    REQUIRE(eq(*sign(x, x), *integer(1)));
    REQUIRE(eq(*sign(exp(x), x), *integer(1)));
    REQUIRE(eq(*sign(exp(neg(x)), x), *integer(1)));
    REQUIRE(eq(*sign(log(x), x), *integer(1)));  // log(x) > 0 for x > 1
}

// --- Compare ---
TEST_CASE("compare() function", "[gruntz_internal]")
{
    auto x = symbol("x");

    // 2 < x < exp(x) < exp(x^2) < exp(exp(x))
    REQUIRE(compare(integer(2), x, x) == '<');
    REQUIRE(compare(x, exp(x), x) == '<');
    REQUIRE(compare(exp(x), exp(pow(x, integer(2))), x) == '<');
    REQUIRE(compare(exp(pow(x, integer(2))), exp(exp(x)), x) == '<');

    // Equality within class: x ~ x^2
    REQUIRE(compare(x, pow(x, integer(2)), x) == '=');
    REQUIRE(compare(pow(x, integer(2)), x, x) == '=');
}

// --- SubsSet ---
TEST_CASE("SubsSet operations", "[gruntz_internal]")
{
    auto x = symbol("x");
    SubsSet s;

    // Insert and retrieve
    auto d1 = s.insert(x);
    REQUIRE(s.dict.size() == 1);
    auto d1_again = s.insert(x);
    REQUIRE(eq(*d1, *d1_again));  // same dummy returned

    // meets
    SubsSet t;
    t.insert(exp(x));
    REQUIRE(!s.meets(t));
    t.insert(x);  // now both have x
    REQUIRE(s.meets(t));

    // do_subs
    SubsSet sub;
    auto d = sub.insert(exp(x));
    auto expr = pow(d, integer(2));  // d^2
    auto result = sub.do_subs(expr);
    REQUIRE(eq(*result, *pow(exp(x), integer(2))));
}
```

#### 6.9 Gruntz Thesis Examples (Pages 122-123)

```cpp
TEST_CASE("Gruntz thesis examples", "[gruntz_thesis]")
{
    auto x = symbol("x");

    // Example from thesis: limit of (exp(x*exp(-x)) - 1)/x
    // Add more examples from thesis pp.122-123 as appropriate
}
```

---

## 7. Memoization Strategy

Many functions (`limitinf`, `mrv_leadterm`, `sign`) are called repeatedly with the same arguments. Without memoization, the algorithm can blow up exponentially.

### Approach

Use `std::unordered_map` with a combined hash of the expression and the symbol.

```cpp
namespace {
    using CacheKey = std::pair<size_t, size_t>;  // hash(e) + hash(x)

    CacheKey make_key(const RCP<const Basic> &e, const RCP<const Symbol> &x) {
        return {e->hash(), x->hash()};
    }
}
```

### Functions to Memoize

| Function | Cache Type | Notes |
|----------|-----------|-------|
| `sign(e, x)` | `map<CacheKey, RCP<const Integer>>` | Must clear between top-level calls |
| `limitinf(e, x)` | `map<CacheKey, RCP<const Basic>>` | Must clear between top-level calls |
| `mrv_leadterm(e, x)` | `map<CacheKey, pair<RCP, RCP>>` | Must clear between top-level calls |

### Thread Safety

For single-threaded use, static local maps are fine. For multi-threaded, use `thread_local` or pass a context object.

### Cache Clearing

All caches must be cleared at the beginning of `gruntz()` (the top-level entry point) to avoid stale entries from previous calls.

```cpp
RCP<const Basic> gruntz(const RCP<const Basic> &e,
                         const RCP<const Symbol> &z,
                         const RCP<const Basic> &z0,
                         const std::string &dir)
{
    // Clear memoization caches
    _gruntz_clear_caches();

    // ... rest of gruntz
}
```

### Implementation Pattern

```cpp
// In anonymous namespace at top of gruntz.cpp
namespace {

struct GruntzCache {
    std::unordered_map<size_t, RCP<const Basic>> limitinf_cache;
    std::unordered_map<size_t, RCP<const Integer>> sign_cache;
    std::unordered_map<size_t,
        std::pair<RCP<const Basic>, RCP<const Basic>>> mrv_leadterm_cache;

    void clear() {
        limitinf_cache.clear();
        sign_cache.clear();
        mrv_leadterm_cache.clear();
    }
};

// Thread-local cache (for multi-threaded safety)
thread_local GruntzCache tls_cache;

} // anonymous namespace
```

### Hash Key Construction

For memoization, the key should cover both `e` and `x`:

```cpp
static size_t make_key2(const RCP<const Basic> &e, const RCP<const Symbol> &x) {
    // Combine two hashes
    size_t h1 = e->hash();
    size_t h2 = x->hash();
    return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
}
```

---

## 8. Debugging and Verification

### Debug Print Macro

Define a conditional debug print (like SymPy's `SYMPY_DEBUG`):

```cpp
#ifdef GRUNTZ_DEBUG
#define GRUNTZ_LOG(msg) std::cerr << "[gruntz] " << msg << std::endl
#else
#define GRUNTZ_LOG(msg) ((void)0)
#endif
```

Compile with `-DGRUNTZ_DEBUG` for verbose output.

### Verification Strategy

1. **Compare with SymPy:** For each test case, compute the limit in both SymEngine and SymPy and compare results.

2. **Unit test each function independently:**
   - `sign()` — known signs for simple expressions
   - `compare()` — known ordering of classes
   - `mrv()` — known MRV sets for simple expressions
   - `rewrite()` — known rewrites
   - `limitinf()` — known limits

3. **Run SymPy reference test file:** Use `/opt-3/cpython-v3.13-apt-deb/lib/python3.13/site-packages/sympy/series/tests/test_gruntz.py` as a reference for expected results.

### Debugging Recursive Functions

Since the algorithm is deeply recursive, add a `depth` counter:

```cpp
thread_local int gruntz_depth = 0;

#define GRUNTZ_ENTER(fn) \
    gruntz_depth++; \
    GRUNTZ_LOG(std::string(gruntz_depth*2, ' ') + "+-" + fn)

#define GRUNTZ_EXIT(fn, result) \
    GRUNTZ_LOG(std::string(gruntz_depth*2, ' ') + "   => " + result); \
    gruntz_depth--
```

---

## 9. Implementation Phases

### Phase 1: Foundation (~400 lines)

**Files:** `gruntz.h`, `gruntz.cpp` (first 400 lines)

**Deliverables:**
- [ ] `SubsSet` class with all methods (`insert`, `meets`, `union_set`, `do_subs`, `empty`, `size`)
- [ ] `has()` helper
- [ ] `independent()` helper (simplified version: return `one` for const part if complex)
- [ ] `powdenest_simple()` helper
- [ ] `tractable_rewrite()` helper (using existing `rewrite_as_exp`)
- [ ] `sign()` function
- [ ] `compare()` function

**Build:** `gruntz.cpp` compiles. CMakeLists updated.

**Tests:** `sign()` and `compare()` unit tests. `SubsSet` unit tests.

### Phase 2: MRV Computation (~250 lines)

**Files:** `gruntz.cpp` (next 250 lines)

**Deliverables:**
- [ ] `mrv()` — the core recursive decomposition
- [ ] `mrv_max3()` — combine same-class MRV sets
- [ ] `mrv_max1()` — combine any two MRV sets
- [ ] `moveup2()` / `moveup()` — move up to exponential level

**Tests:** `mrv()` unit tests with known MRV sets.

### Phase 3: Rewrite and Lead Term (~200 lines)

**Files:** `gruntz.cpp` (next 200 lines)

**Deliverables:**
- [ ] `build_expression_tree()` — dependency ordering
- [ ] `rewrite()` — the complex rewriting function
- [ ] `mrv_leadterm()` — leading term extraction

**Tests:** `rewrite()` unit tests, `mrv_leadterm()` unit tests.

### Phase 4: Limit Computation and Main Entry (~150 lines)

**Files:** `gruntz.cpp` (next 150 lines), `test_gruntz.cpp`

**Deliverables:**
- [ ] `limitinf()` — limit at infinity
- [ ] `gruntz()` — main entry point with limit transformation
- [ ] Memoization (optional at this phase; add if needed for performance)
- [ ] Full integration test suite

**Tests:** All test categories: basic, exponential, logarithmic, power, MRV, directional, negative infinity.

### Phase 5: Polish and Edge Cases

- [ ] Handle complex expressions (with `I`)
- [ ] Handle trig/hyperbolic functions via `tractable_rewrite`
- [ ] Handle `exp(log(...))` simplification
- [ ] Performance optimization (memoization if not yet added)
- [ ] Edge case: expressions that are identically 0/∞
- [ ] Error handling: throw `NotImplementedError` for unsupported types

---

## 10. Common Pitfalls and Checklist

### Pitfalls to Avoid

| # | Pitfall | Solution |
|---|---------|----------|
| 1 | Using `subs` instead of `xreplace` | Use `xreplace` for structural replacement; `subs` may mathematically simplify, which breaks the algorithm |
| 2 | Not handling `exp(log(x))` → `x` | In `mrv()` case for exp, check if exponent is `Log` and simplify immediately |
| 3 | Forgetting to memoize `limitinf` | Mutual recursion between `compare` → `limitinf` → `mrv_leadterm` → `sign` → `compare` will loop infinitely |
| 4 | Confusing `E^arg` vs `exp` class | There is no `Exp` class. `exp(x)` = `pow(E, x)`. Check `is_a<Pow>` and `base == E` |
| 5 | Not handling empty MRV set | Check `Omega.empty()` before calling `rewrite` |
| 6 | Integer overflow in LCM | Use `integer_class` (GMP mpz) for arbitrary precision |
| 7 | Series precision too low | Start with 6, increase if lead coefficient depends on `w` or is zero |
| 8 | Not clearing caches between top-level calls | Memoization caches must be thread-local and cleared in `gruntz()` |
| 9 | Casting `RCP<const Basic>` to specific types | Use `down_cast<const T &>(*ptr)` for checked casts |
| 10 | `as_two_terms()` not available | Both `Add` and `Mul` have `as_two_terms(Ptr<RCP<const Basic>> a, Ptr<RCP<const Basic>> b)` |

### Pre-Implementation Checklist

- [ ] Read SymPy `gruntz.py` fully (702 lines)
- [ ] Understand the exponential comparability class hierarchy
- [ ] Understand the role of `w` and why `w → 0`
- [ ] Understand why `xreplace` not `subs`
- [ ] Trace through one complete example on paper: `lim_{x→∞} exp(x) / exp(x + exp(-x))`
- [ ] Understand `mrv_max3` vs `mrv_max1` difference
- [ ] Understand `moveup2` — when and why x must be exponentiated

### Post-Implementation Checklist

- [ ] All test cases from Section 6 pass
- [ ] `make` or `cmake --build build` succeeds with no errors
- [ ] `ctest -R test_gruntz` passes all tests
- [ ] Compare a few results against SymPy for verification
- [ ] No memory leaks (valgrind if available)

---

## Appendix A: SymEngine API Quick Reference

| Operation | C++ Code |
|-----------|----------|
| Create symbol | `auto x = symbol("x");` |
| Create dummy | `auto d = dummy("d");` |
| Integer constant | `integer(42)` |
| Constants | `zero`, `one`, `minus_one`, `E`, `Inf`, `NegInf` |
| Addition | `add(a, b)` |
| Subtraction | `sub(a, b)` |
| Multiplication | `mul(a, b)` |
| Division | `div(a, b)` |
| Power | `pow(a, b)` |
| Negation | `neg(a)` |
| Exponential | `exp(a)` = `pow(E, a)` |
| Logarithm | `log(a)` |
| Sine/Cosine | `sin(a)`, `cos(a)` |
| Hyperbolic | `sinh(a)`, `cosh(a)`, `tanh(a)` |
| Rewrite trig → exp | `rewrite_as_exp(expr)` |
| Check contains symbol | `has_symbol(expr, x)` |
| Tree substitution | `expr->xreplace({{old, new}})` |
| Math substitution | `expr->subs({{old, new}})` |
| Expand expression | `expand(expr)` |
| Type check | `is_a<Pow>(*expr)` |
| Down-cast | `down_cast<const Pow &>(*expr)` |
| Compare equality | `eq(*a, *b)` |
| String print | `std::cout << *expr` |
| Series expansion | `series(expr, var, prec)` → `.as_dict()` |
| Throw not-implemented | `throw NotImplementedError("msg");` |
| Debug assert | `SYMENGINE_ASSERT(condition);` |
| Create RCP | `make_rcp<const T>(args...)` |

---

## Appendix B: SymPy → SymEngine Type Mapping

| SymPy | SymEngine | Notes |
|-------|-----------|-------|
| `S.Zero` | `integer(0)` / `zero` | Predefined constant |
| `S.One` | `integer(1)` / `one` | Predefined constant |
| `S.NegativeOne` | `integer(-1)` / `minus_one` | Predefined constant |
| `S.Exp1` | `E` | Euler's number |
| `S.Infinity` / `oo` | `Inf` | Unsigned infinity |
| `S.NegativeInfinity` / `-oo` | `NegInf` | Negative infinity |
| `Symbol('x')` | `symbol("x")` | |
| `Dummy('d')` | `dummy("d")` | |
| `exp(x)` | `exp(x)` = `pow(E, x)` | No separate Exp class |
| `log(x)` | `log(x)` | |
| `x + y` | `add(x, y)` | |
| `x * y` | `mul(x, y)` | |
| `x / y` | `div(x, y)` | |
| `x ** y` | `pow(x, y)` | |
| `-x` | `neg(x)` | |
| `x.as_two_terms()` | `a.as_two_terms(outArg(t1), outArg(t2))` | Uses `Ptr<>` |
| `e.has(x)` | `has_symbol(e, x)` | From visitor.h |
| `e.subs({a: b})` | `e->subs({{a, b}})` | |
| `e.xreplace({a: b})` | `e->xreplace({{a, b}})` | |
| `e.is_Pow` | `is_a<Pow>(*e)` | |
| `e.is_Mul` | `is_a<Mul>(*e)` | |
| `e.is_Add` | `is_a<Add>(*e)` | |
| `isinstance(e, log)` | `is_a<Log>(*e)` | |
| `e.base` (Pow) | `down_cast<const Pow &>(*e).get_base()` | |
| `e.exp` (Pow) | `down_cast<const Pow &>(*e).get_exp()` | |
| `e.args[0]` (Log) | `down_cast<const Log &>(*e).get_arg()` | |
| `Series.leadterm(w)` | `series(e, w, prec)->as_dict()` | Extract first entry |
| `ilcm(a, b, c)` | Manual LCM using `mp_lcm` | |

---

## Appendix C: Example Trace

For `lim_{x→∞} exp(x) / exp(x + exp(-x)) → 1`:

```
gruntz(e = exp(x)/exp(x + exp(-x)), z = x, z0 = oo)
 |
 +-limitinf(e = exp(x)/exp(x + exp(-x)), x)
    |
    +-mrv_leadterm(e, x)
       |
       +-mrv(e, x)
       |  // e = pow(E,x) / pow(E, x + pow(E, neg(x)))
       |  // is Mul -> split into two terms
       |  // t1 = exp(x), t2 = 1/exp(x + exp(-x))
       |  //
       |  // mrv(exp(x), x) -> Omega1 = {exp(x)}, exps1 = d1
       |  // mrv(1/exp(x+exp(-x)), x) -> Omega2 = {exp(x+exp(-x))}, exps2 = 1/d2
       |  //
       |  // mrv_max1(Omega1, Omega2, d1 * 1/d2, x)
       |  //   -> union computes {exp(x): d1, exp(x+exp(-x)): d2}
       |  //   -> mrv_max3 compares exp(x) vs exp(x+exp(-x))
       |  //      -> compare: L = limitinf(x/(x+exp(-x)), x) = 1
       |  //      -> result: = (equal class)
       |  //      -> return (union, d1/d2)
       |
       +-rewrite(exps = d1/d2, Omega = {exp(x):d1, exp(x+exp(-x)):d2}, x, w)
          // Omega sorted: [exp(x+exp(-x)), exp(x)]  (higher height first)
          // g = exp(x), sig = sign(x, x) = 1 -> wsym = 1/w
          //
          // For f=exp(x+exp(-x)):
          //   c = limitinf((x+exp(-x))/x, x) = 1
          //   arg = x + exp(-x)
          //   rewritten = exp((x+exp(-x)) - 1*x) * w^1 = exp(exp(-x)) * w
          //
          // For f=exp(x):
          //   c = limitinf(x/x, x) = 1
          //   arg = x
          //   rewritten = exp(x - 1*x) * w = exp(0) * w = w
          //
          // Substitute: d1 → w, d2 → exp(exp(-x))*w
          // But wait — this isn't right yet. The rewrite needs logw.
          //
          // After full rewrite: f = 1/exp(exp(-x))
          // Hmm, this is complex. See SymPy trace for details.
          //
          // Final rewritten: f = exp(-w), logw = -x
          |
          +-series(f = exp(-w), w, 6)
             -> series = 1 - w + w^2/2 - w^3/6 + ...
             -> leadterm = (1, 0)
             -> return (c0=1, e0=0)
       |
       +-sign(e0=0, x) = 0
       +-limitinf(1, x) = 1
       |
       return 1
```

---

*End of Plan*
