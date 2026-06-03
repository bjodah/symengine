# Gruntz Algorithm Implementation for SymEngine — Detailed Plan (v3)

> **Target audience:** Junior C++ developer familiar with SymEngine basics.
> **Reference:** SymPy's `gruntz.py` — `/opt-3/cpython-v3.13-apt-deb/lib/python3.13/site-packages/sympy/series/gruntz.py` (703 lines).
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
   - [4.8 `mrv()`](#48-mrv--most-rapidly-varying-set)
   - [4.9 `mrv_max3()`](#49-mrv_max3)
   - [4.10 `mrv_max1()`](#410-mrv_max1)
   - [4.11 `moveup2()` and `moveup()`](#411-moveup2-and-moveup)
   - [4.12 `build_expression_tree()`](#412-build_expression_tree)
   - [4.13 `rewrite()`](#413-rewrite)
   - [4.14 `mrv_leadterm()`](#414-mrv_leadterm)
   - [4.15 `limitinf()`](#415-limitinf)
   - [4.16 `gruntz()`](#416-gruntz--main-entry-point)
5. [Build Integration](#5-build-integration)
6. [Test Suite: `test_gruntz.cpp`](#6-test-suite-test_gruntzcpp)
7. [Memoization Strategy](#7-memoization-strategy)
8. [Debugging and Verification](#8-debugging-and-verification)
9. [Implementation Phases](#9-implementation-phases)
10. [Common Pitfalls and Checklist](#10-common-pitfalls-and-checklist)

---

## 1. Algorithm Overview

The Gruntz algorithm computes $\lim_{z \to z_0} f(z)$ by:

```
Step 1: Transform any limit to x -> +oo
        (z->oo: identity,  z->-oo: z->-z,  z->a: z->a ± 1/x)

Step 2: Find the MRV (Most Rapidly Varying) subexpressions of f(x)
        -> returns Omega = SubsSet mapping: {expr -> dummy}, plus rewritten f

Step 3: Rewrite f in terms of a new variable w -> 0
        -> every element of Omega gets rewritten in terms of w and logw, then xreplace

Step 4: Compute series expansion of the rewritten f in w
        -> extract leading term: c0 * w^e0 (requires searching unordered map for minimum exponent)

Step 5: Read off the limit:
        e0 > 0  ->  limit = 0
        e0 < 0  ->  limit = sign(c0) * oo
        e0 = 0  ->  recurse: limit = lim c0
```

**Comparability classes** (slowest $\to$ fastest):
```
2  ~  3  ~  -5        [constants]
x  ~  x²  ~  1/x      [polynomial-like]
exp(x)  ~  exp(-x)    [exponential-like]
exp(x²)               [each new exponent gives new class]
exp(exp(x))           [doubly exponential]
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
// Returns const_part, dependent is populated with the dependent part.
// Matches standard SymEngine style for output pointer arguments.
RCP<const Basic> independent(const RCP<const Basic> &e,
                             const RCP<const Symbol> &x,
                             const Ptr<RCP<const Basic>> &dependent);

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
// Returns nodes map: Dummy -> height
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

---

## 4. Implementation: `gruntz.cpp` — Functions in Detail

### 4.1 `SubsSet` Class

#### `RCP<const Symbol> insert(const RCP<const Basic> &expr)`

Creates a dummy symbol safely:

```cpp
RCP<const Symbol> SubsSet::insert(const RCP<const Basic> &expr)
{
    auto it = dict.find(expr);
    if (it != dict.end()) {
        return rcp_static_cast<const Symbol>(it->second);
    }
    // Predictable dummy naming
    RCP<const Symbol> d = dummy(); 
    dict[expr] = d;
    return d;
}
```

#### `bool meets(const SubsSet &other) const`

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

```cpp
std::pair<SubsSet, RCP<const Basic>>
SubsSet::union_set(const SubsSet &other, const RCP<const Basic> &exps) const
{
    SubsSet res = *this;
    map_basic_basic tr;
    RCP<const Basic> res_exps = exps;

    for (const auto &kv : other.dict) {
        const auto &expr = kv.first;
        const auto &var = kv.second;
        auto it = res.dict.find(expr);
        if (it != res.dict.end()) {
            tr[var] = it->second;
            if (res_exps != nullptr) {
                res_exps = res_exps->xreplace(tr);
            }
        } else {
            res.dict[expr] = var;
        }
    }
    for (const auto &kv : other.rewrites) {
        res.rewrites[kv.first] = kv.second->xreplace(tr);
    }
    return {res, res_exps};
}
```

#### `RCP<const Basic> do_subs(const RCP<const Basic> &e) const`

```cpp
RCP<const Basic> SubsSet::do_subs(const RCP<const Basic> &e) const
{
    map_basic_basic rev;
    for (const auto &kv : dict) {
        rev[kv.second] = kv.first;
    }
    return e->xreplace(rev);
}
```

---

### 4.2 `has()` Helper

```cpp
bool has(const Basic &b, const Basic &x)
{
    return has_symbol(b, x);
}
```

---

### 4.3 `independent()` Helper

Splits terms canonicalizing `umap_basic_num` and `map_basic_basic` using `std::move`.

```cpp
RCP<const Basic> independent(const RCP<const Basic> &e,
                             const RCP<const Symbol> &x,
                             const Ptr<RCP<const Basic>> &dependent)
{
    if (is_a<Add>(*e)) {
        const Add &a = down_cast<const Add &>(*e);
        RCP<const Number> coef = a.get_coef();
        umap_basic_num const_dict;
        umap_basic_num var_dict;
        for (const auto &kv : a.get_dict()) {
            if (has_symbol(*kv.first, *x)) {
                var_dict[kv.first] = kv.second;
            } else {
                const_dict[kv.first] = kv.second;
            }
        }
        if (const_dict.empty()) {
            *dependent = e;
            return zero;
        }
        if (var_dict.empty()) {
            *dependent = zero;
            return e;
        }
        *dependent = Add::from_dict(zero, std::move(var_dict));
        return Add::from_dict(coef, std::move(const_dict));
    }
    
    if (is_a<Mul>(*e)) {
        const Mul &m = down_cast<const Mul &>(*e);
        RCP<const Number> coef = m.get_coef();
        map_basic_basic const_dict;
        map_basic_basic var_dict;
        for (const auto &kv : m.get_dict()) {
            if (has_symbol(*kv.first, *x)) {
                // Mul::dict_add_term parameter order is (dict, exp, base)
                Mul::dict_add_term(var_dict, kv.second, kv.first);
            } else {
                Mul::dict_add_term(const_dict, kv.second, kv.first);
            }
        }
        if (const_dict.empty()) {
            *dependent = e;
            return one;
        }
        if (var_dict.empty()) {
            *dependent = one;
            return e;
        }
        *dependent = Mul::from_dict(one, std::move(var_dict));
        return Mul::from_dict(coef, std::move(const_dict));
    }

    if (has_symbol(*e, *x)) {
        *dependent = e;
        return one;
    }
    *dependent = one;
    return e;
}
```

---

### 4.4 `powdenest_simple()` Helper

Implemented cleanly via `TransformVisitor` subclassing.

```cpp
class PowDenestVisitor : public TransformVisitor {
public:
    using TransformVisitor::bvisit;

    void bvisit(const Pow &x) override {
        auto base = apply(x.get_base());
        auto expn = apply(x.get_exp());
        if (is_a<Pow>(*base)) {
            const Pow &inner = down_cast<const Pow &>(*base);
            result_ = pow(inner.get_base(), mul(inner.get_exp(), expn));
        } else {
            result_ = pow(base, expn);
        }
    }
};

RCP<const Basic> powdenest_simple(const RCP<const Basic> &e)
{
    PowDenestVisitor visitor;
    return visitor.apply(e);
}
```

---

### 4.5 `tractable_rewrite()` Helper

```cpp
RCP<const Basic> tractable_rewrite(const RCP<const Basic> &e)
{
    return rewrite_as_exp(e);
}
```

---

### 4.6 `sign()` — Determine Sign as x → ∞

```cpp
RCP<const Integer> sign(const RCP<const Basic> &e,
                         const RCP<const Symbol> &x)
{
    if (is_a_Number(*e)) {
        const Number &num = down_cast<const Number &>(*e);
        if (num.is_positive()) return integer(1);
        if (num.is_negative()) return integer(-1);
        return integer(0);
    }

    if (!has_symbol(*e, *x)) {
        // Evaluate numerical value if expression has no x dependencies
        // (fallback to leadterm if not possible)
    }

    if (eq(*e, *x)) {
        return integer(1);
    }

    if (is_a<Pow>(*e)) {
        const Pow &p = down_cast<const Pow &>(*e);
        if (eq(*p.get_base(), *E)) {
            return integer(1); // exp(...) is positive
        }
        auto sb = sign(p.get_base(), x);
        if (eq(*sb, *one)) {
            return integer(1);
        }
        if (is_a<Integer>(*p.get_exp())) {
            int exp_val = down_cast<const Integer &>(*p.get_exp()).as_int();
            if (exp_val % 2 == 0) return integer(1);
        }
    }

    if (is_a<Mul>(*e)) {
        RCP<const Basic> a, b;
        down_cast<const Mul &>(*e).as_two_terms(outArg(a), outArg(b));
        auto sa = sign(a, x);
        if (eq(*sa, *zero)) return integer(0);
        return rcp_static_cast<const Integer>(mul(sa, sign(b, x)));
    }

    if (is_a<Log>(*e)) {
        RCP<const Basic> arg = down_cast<const Log &>(*e).get_arg();
        return sign(sub(arg, one), x);
    }

    // Leadterm fallback
    RCP<const Basic> c0, e0;
    std::tie(c0, e0) = mrv_leadterm(e, x);
    return sign(c0, x);
}
```

---

### 4.7 `compare()` — Compare Growth Rates

```cpp
char compare(const RCP<const Basic> &a,
             const RCP<const Basic> &b,
             const RCP<const Symbol> &x)
{
    RCP<const Basic> la = log(a);
    RCP<const Basic> lb = log(b);

    if (is_a<Pow>(*a) && eq(*down_cast<const Pow &>(*a).get_base(), *E)) {
        la = down_cast<const Pow &>(*a).get_exp();
    }
    if (is_a<Pow>(*b) && eq(*down_cast<const Pow &>(*b).get_base(), *E)) {
        lb = down_cast<const Pow &>(*b).get_exp();
    }

    RCP<const Basic> L = limitinf(div(la, lb), x);

    if (eq(*L, *zero)) return '<';
    if (is_a<Infty>(*L)) return '>';
    return '=';
}
```

---

### 4.8 `mrv()` — Most Rapidly Varying Set

```cpp
std::pair<SubsSet, RCP<const Basic>>
mrv(const RCP<const Basic> &e, const RCP<const Symbol> &x)
{
    if (!has_symbol(*e, *x)) {
        return {SubsSet(), e};
    }

    if (eq(*e, *x)) {
        SubsSet s;
        auto d = s.insert(x);
        return {s, d};
    }

    if (is_a<Mul>(*e) || is_a<Add>(*e)) {
        RCP<const Basic> dep_part;
        RCP<const Basic> const_part = independent(e, x, outArg(dep_part));
        if (const_part != one && const_part != zero) {
            SubsSet s;
            RCP<const Basic> expr;
            std::tie(s, expr) = mrv(dep_part, x);
            if (is_a<Mul>(*e)) {
                return {s, mul(const_part, expr)};
            } else {
                return {s, add(const_part, expr)};
            }
        }
        
        RCP<const Basic> a, b;
        if (is_a<Mul>(*e)) {
            down_cast<const Mul &>(*e).as_two_terms(outArg(a), outArg(b));
        } else {
            down_cast<const Add &>(*e).as_two_terms(outArg(a), outArg(b));
        }
        
        SubsSet s1, s2;
        RCP<const Basic> e1, e2;
        std::tie(s1, e1) = mrv(a, x);
        std::tie(s2, e2) = mrv(b, x);
        
        RCP<const Basic> reconstructed = is_a<Mul>(*e) ? mul(e1, e2) : add(e1, e2);
        return mrv_max1(s1, s2, reconstructed, x);
    }

    if (is_a<Pow>(*e) && !eq(*down_cast<const Pow &>(*e).get_base(), *E)) {
        RCP<const Basic> e1 = one;
        RCP<const Basic> current = e;
        while (is_a<Pow>(*current)) {
            const Pow &p = down_cast<const Pow &>(*current);
            e1 = mul(e1, p.get_exp());
            current = p.get_base();
        }
        if (eq(*current, *one)) {
            return {SubsSet(), current};
        }
        if (has_symbol(*e1, *x)) {
            return mrv(exp(mul(e1, log(current))), x);
        } else {
            SubsSet s;
            RCP<const Basic> expr;
            std::tie(s, expr) = mrv(current, x);
            return {s, pow(expr, e1)};
        }
    }

    if (is_a<Log>(*e)) {
        RCP<const Basic> arg = down_cast<const Log &>(*e).get_arg();
        SubsSet s;
        RCP<const Basic> expr;
        std::tie(s, expr) = mrv(arg, x);
        return {s, log(expr)};
    }

    // Exponential Case
    if (is_a<Pow>(*e) && eq(*down_cast<const Pow &>(*e).get_base(), *E)) {
        RCP<const Basic> exp_arg = down_cast<const Pow &>(*e).get_exp();
        if (is_a<Log>(*exp_arg)) {
            return mrv(down_cast<const Log &>(*exp_arg).get_arg(), x);
        }

        RCP<const Basic> li = limitinf(exp_arg, x);
        if (is_a<Infty>(*li)) {
            SubsSet s1;
            auto e1 = s1.insert(e);
            SubsSet s2;
            RCP<const Basic> e2;
            std::tie(s2, e2) = mrv(exp_arg, x);
            
            auto su = s1.union_set(s2, nullptr).first;
            su.rewrites[e1] = exp(e2);
            return mrv_max3(s1, e1, s2, exp(e2), su, e1, x);
        } else {
            SubsSet s;
            RCP<const Basic> expr;
            std::tie(s, expr) = mrv(exp_arg, x);
            return {s, exp(expr)};
        }
    }

    throw NotImplementedError("mrv: unhandled function type");
}
```

---

### 4.9 `mrv_max3()`

```cpp
std::pair<SubsSet, RCP<const Basic>>
mrv_max3(const SubsSet &f, const RCP<const Basic> &expsf,
         const SubsSet &g, const RCP<const Basic> &expsg,
         const SubsSet &_union, const RCP<const Basic> &expsboth,
         const RCP<const Symbol> &x)
{
    if (f.empty()) return {g, expsg};
    if (g.empty()) return {f, expsf};
    if (f.meets(g)) return {_union, expsboth};

    char c = compare(f.dict.begin()->first, g.dict.begin()->first, x);
    if (c == '>') return {f, expsf};
    if (c == '<') return {g, expsg};
    return {_union, expsboth};
}
```

---

### 4.10 `mrv_max1()`

```cpp
std::pair<SubsSet, RCP<const Basic>>
mrv_max1(const SubsSet &f, const SubsSet &g,
         const RCP<const Basic> &exps, const RCP<const Symbol> &x)
{
    SubsSet u;
    RCP<const Basic> b;
    std::tie(u, b) = f.union_set(g, exps);
    return mrv_max3(f, g.do_subs(exps),
                    g, f.do_subs(exps),
                    u, b, x);
}
```

---

### 4.11 `moveup2()` and `moveup()`

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

RCP<const Basic> moveup(const RCP<const Basic> &e,
                         const RCP<const Symbol> &x)
{
    map_basic_basic xmap;
    xmap[x] = exp(x);
    return e->xreplace(xmap);
}
```

---

### 4.12 `build_expression_tree()`

```cpp
std::map<RCP<const Basic>, int, RCPBasicKeyLess>
build_expression_tree(const SubsSet &Omega)
{
    std::map<RCP<const Basic>, int, RCPBasicKeyLess> nodes;

    for (const auto &kv : Omega.dict) {
        nodes[kv.second] = 1;
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

### 4.13 `rewrite()`

LCM denominator helper fold over integer denominators:

```cpp
static integer_class compute_denominators_lcm(const std::vector<integer_class> &denoms) {
    integer_class result(1);
    for (const auto &d : denoms) {
        result = mp_lcm(result, d);
    }
    return result;
}
```

```cpp
std::pair<RCP<const Basic>, RCP<const Basic>>
rewrite(const RCP<const Basic> &e,
        const SubsSet &Omega,
        const RCP<const Symbol> &x,
        const RCP<const Symbol> &wsym)
{
    SYMENGINE_ASSERT(!Omega.empty());
    
    auto nodes = build_expression_tree(Omega);
    std::vector<std::pair<RCP<const Basic>, RCP<const Basic>>> Omega_items;
    for (const auto &kv : Omega.dict) {
        Omega_items.push_back({kv.first, kv.second});
    }
    
    // Sort Omega descending by tree dependency heights
    std::sort(Omega_items.begin(), Omega_items.end(),
              [&nodes](const auto &a, const auto &b) {
                  return nodes[a.second] > nodes[b.second];
              });

    // The simplest exp-like term in MRV class is the representative
    RCP<const Basic> g = Omega_items.back().first;
    RCP<const Pow> g_pow = rcp_static_cast<const Pow>(g);
    RCP<const Basic> g_exp = g_pow->get_exp();

    auto sig = sign(g_exp, x);
    
    RCP<const Basic> w_val = wsym;
    if (eq(*sig, *one)) {
        w_val = pow(wsym, minus_one); // If representative goes to oo, replace with 1/w
    }

    std::vector<std::pair<RCP<const Basic>, RCP<const Basic>>> O2;
    std::vector<integer_class> denominators;

    for (const auto &item : Omega_items) {
        RCP<const Basic> f = item.first;
        RCP<const Basic> var = item.second;
        RCP<const Basic> f_exp = rcp_static_cast<const Pow>(f)->get_exp();

        RCP<const Basic> c = limitinf(div(f_exp, g_exp), x);

        if (is_a<Rational>(*c)) {
            denominators.push_back(down_cast<const Rational &>(*c).as_rational_class().get_den());
        }

        RCP<const Basic> arg = f_exp;
        auto rit = Omega.rewrites.find(var);
        if (rit != Omega.rewrites.end()) {
            arg = rcp_static_cast<const Pow>(rit->second)->get_exp();
        }

        RCP<const Basic> rewritten = mul(exp(sub(arg, mul(c, g_exp))), pow(w_val, c));
        O2.push_back({var, rewritten});
    }

    RCP<const Basic> f_expr = powdenest_simple(e);
    for (const auto &sub_rule : O2) {
        map_basic_basic repl;
        repl[sub_rule.first] = sub_rule.second;
        f_expr = f_expr->xreplace(repl);
    }

    // Verify all mrv dummies are completely removed
    for (const auto &item : Omega_items) {
        SYMENGINE_ASSERT(!has_symbol(*f_expr, *item.second));
    }

    RCP<const Basic> logw = g_exp;
    if (eq(*sig, *one)) {
        logw = neg(logw);
    }

    integer_class lcm_val = compute_denominators_lcm(denominators);
    map_basic_basic w_norm;
    w_norm[wsym] = pow(wsym, integer(lcm_val));
    f_expr = f_expr->subs(w_norm);
    logw = div(logw, integer(lcm_val));

    f_expr = expand(f_expr);
    return {f_expr, logw};
}
```

---

### 4.14 `mrv_leadterm()`

Implements the **Singular Log substitution workaround** and correct **Unordered Map iteration** to find the minimum exponent.

```cpp
std::pair<RCP<const Basic>, RCP<const Basic>>
mrv_leadterm(const RCP<const Basic> &e,
             const RCP<const Symbol> &x)
{
    if (!has_symbol(*e, *x)) {
        return {e, zero};
    }

    SubsSet Omega;
    RCP<const Basic> exps;
    std::tie(Omega, exps) = mrv(e, x);

    if (Omega.empty()) {
        return {exps, zero};
    }

    if (Omega.dict.find(x) != Omega.dict.end()) {
        Omega = moveup2(Omega, x);
        exps = moveup(exps, x);
    }

    auto w = dummy("w"); // No arguments assumptions allowed in C++
    RCP<const Basic> f, logw;
    std::tie(f, logw) = rewrite(exps, Omega, x, w);

    // Singular Log Workaround: Substitutes log(w) with _logw_dummy
    auto logw_dummy = symbol("_logw_dummy");
    map_basic_basic log_repl;
    log_repl[log(w)] = logw_dummy;
    
    RCP<const Basic> f_for_series = f->xreplace(log_repl);

    unsigned int prec = 6;
    auto series_result = series(f_for_series, w, prec);
    auto series_dict = series_result->as_dict(); // Returns umap_int_basic (unordered)

    // Robust search for the leading term (minimum exponent with non-zero coefficient)
    int lead_exp = -1;
    RCP<const Basic> lead_coeff = zero;
    bool found = false;
    
    for (const auto &kv : series_dict) {
        if (!eq(*kv.second, *zero)) {
            if (!found || kv.first < lead_exp) {
                lead_exp = kv.first;
                lead_coeff = kv.second;
                found = true;
            }
        }
    }
    
    if (!found) {
        lead_coeff = zero;
        lead_exp = 0;
    }

    // Restore actual logw expression into lead coefficient
    map_basic_basic log_restore;
    log_restore[logw_dummy] = logw;
    lead_coeff = lead_coeff->xreplace(log_restore);

    return {lead_coeff, integer(lead_exp)};
}
```

---

### 4.15 `limitinf()`

Uses `is_a_Symbol` instead of strict `is_a<Symbol>`.

```cpp
RCP<const Basic> limitinf(const RCP<const Basic> &e,
                           const RCP<const Symbol> &x)
{
    if (!has_symbol(*e, *x)) {
        return e;
    }

    RCP<const Basic> e_rewritten = tractable_rewrite(e);
    e_rewritten = powdenest_simple(e_rewritten);

    RCP<const Basic> c0, e0;
    std::tie(c0, e0) = mrv_leadterm(e_rewritten, x);

    auto sig = sign(e0, x);

    if (eq(*sig, *one)) {
        return zero; // e0 > 0: lim = 0
    }
    
    if (eq(*sig, *minus_one)) {
        auto sc = sign(c0, x);
        if (eq(*sc, *zero)) {
            throw std::runtime_error("Leading term sign should not be zero");
        }
        if (eq(*sc, *one)) return Inf;
        return NegInf;
    }

    if (eq(*sig, *zero)) {
        return limitinf(c0, x); // e0 = 0: recurse on c0
    }

    throw std::runtime_error("Unexpected sign result");
}
```

---

### 4.16 `gruntz()` — Main Entry Point

```cpp
RCP<const Basic> gruntz(const RCP<const Basic> &e,
                         const RCP<const Symbol> &z,
                         const RCP<const Basic> &z0,
                         const std::string &dir)
{
    // Validate: handles Symbol and Dummy correctly
    if (!is_a_Symbol(*z)) {
        throw NotImplementedError("Second argument must be a Symbol");
    }

    // Clear thread-local memoization caches
    _gruntz_clear_caches();

    RCP<const Basic> e0;
    if (is_a<Infty>(*z0)) {
        const Infty &inf = down_cast<const Infty &>(*z0);
        if (inf.is_positive_infinity()) {
            e0 = e;
        } else if (inf.is_negative_infinity()) {
            map_basic_basic subst_map;
            subst_map[z] = neg(z);
            e0 = e->subs(subst_map);
        } else {
            throw NotImplementedError("Unsigned infinity not supported");
        }
    } else {
        map_basic_basic subst_map;
        if (dir == "+") {
            subst_map[z] = add(z0, div(one, z));
        } else if (dir == "-") {
            subst_map[z] = sub(z0, div(one, z));
        } else {
            throw NotImplementedError("dir must be '+' or '-'");
        }
        e0 = e->subs(subst_map);
    }

    return limitinf(e0, z);
}
```

---

## 5. Build Integration

Identical to Plan (v2). See Section 5 of [02-PLAN-GRUNTZ-v2.md](file:///work/02-PLAN-GRUNTZ-v2.md).

---

## 6. Test Suite: `test_gruntz.cpp`

Identical to Plan (v2). See Section 6 of [02-PLAN-GRUNTZ-v2.md](file:///work/02-PLAN-GRUNTZ-v2.md).

---

## 7. Memoization Strategy

Using `RCPBasicHash` and `RCPBasicKeyEq` for thread-local, collision-safe maps keyed by `RCP<const Basic>` expression key:

```cpp
namespace {

struct GruntzCache {
    std::unordered_map<RCP<const Basic>, RCP<const Basic>, RCPBasicHash, RCPBasicKeyEq> limitinf_cache;
    std::unordered_map<RCP<const Basic>, RCP<const Integer>, RCPBasicHash, RCPBasicKeyEq> sign_cache;
    std::unordered_map<RCP<const Basic>, std::pair<RCP<const Basic>, RCP<const Basic>>, RCPBasicHash, RCPBasicKeyEq> mrv_leadterm_cache;

    void clear() {
        limitinf_cache.clear();
        sign_cache.clear();
        mrv_leadterm_cache.clear();
    }
};

thread_local GruntzCache tls_cache;

void _gruntz_clear_caches() {
    tls_cache.clear();
}

} // anonymous namespace
```

---

## 8. Debugging and Verification

Identical to Plan (v2). See Section 8 of [02-PLAN-GRUNTZ-v2.md](file:///work/02-PLAN-GRUNTZ-v2.md).

---

## 9. Implementation Phases

Identical to Plan (v2). See Section 9 of [02-PLAN-GRUNTZ-v2.md](file:///work/02-PLAN-GRUNTZ-v2.md).

---

## 10. Common Pitfalls and Checklist

### Pitfalls to Avoid

| # | Pitfall | Solution |
|---|---------|----------|
| 1 | Iterating unordered series dictionary sequentially | Search `umap_int_basic` to locate the minimum exponent. |
| 2 | Evaluating `log(w)` at `w=0` | Use the Singular Log dummy replacement workaround prior to running `series()`. |
| 3 | Assigning `Pow` to `RCP<const Symbol>` | Use `RCP<const Basic> w_val` intermediate variables for inverted dummy values. |
| 4 | Passing out temporary `Ptr` as non-const ref | Declare `independent` parameter as `const Ptr<RCP<const Basic>> &dependent`. |
| 5 | Swapping `Mul::dict_add_term` arguments | Exponent goes *second*, term goes *third*. |
| 6 | Unsafe combined hash cache keys | Key maps simply on `RCP<const Basic>` since $x$ is invariant during recursive limit computations. |
| 7 | `is_a<Symbol>` rejecting `Dummy` | Use `is_a_Symbol` to validate the variable. |
| 8 | Complex manual tree walkers | Subclass `TransformVisitor` to implement `powdenest_simple` cleanly. |

---
*End of Plan*
