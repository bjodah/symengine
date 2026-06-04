# Review of Plan v3: Gruntz Algorithm Implementation for SymEngine

> Review of [04-PLAN-GRUNTZ-v3.md](./04-PLAN-GRUNTZ-v3.md) against the existing SymEngine codebase at `/work/symengine/`.

---

## 1. Executive Summary

The plan is **substantially correct and implementable**. The vast majority of the SymEngine infrastructure it relies on already exists. The primary work is **new code** (the `SubsSet` class and all algorithm functions), not modification of existing SymEngine internals. Below, every dependency from the plan is validated against the actual codebase.

---

## 2. What Already Exists (Verified)

| Dependency | Location | Status |
|---|---|---|
| `has_symbol(b, x)` | `visitor.h:154` | ✅ Exists |
| `dummy()` / `dummy(name)` | `symbol.h:93-106` | ✅ Exists, returns `RCP<const Dummy>` (subclass of `Symbol`) |
| `is_a_Symbol(b)` | `symbol.h:108` | ✅ Exists, accepts `Symbol` and `Dummy` |
| `is_a_Number(b)` | `number.h:130` | ✅ Exists |
| `is_a<T>(b)` | `basic.h:252` | ✅ Exists |
| `down_cast<T>(b)` | `symengine_casts.h` | ✅ Exists |
| `rcp_static_cast<T>(p)` | `symengine_rcp.h` | ✅ Exists |
| `RCPBasicHash` | `basic.h:210` | ✅ Exists |
| `RCPBasicKeyEq` | `basic.h:219` | ✅ Exists |
| `RCPBasicKeyLess` | `basic.h:228` / `dict.h:28` | ✅ Exists |
| `map_basic_basic` | `dict.h:57` | ✅ Exists (`std::map<... , RCPBasicKeyLess>`) |
| `umap_basic_num` | `dict.h:35` | ✅ Exists |
| `umap_int_basic` | `dict.h:37` | ✅ Exists |
| `set_basic` | `dict.h:50` | ✅ Exists |
| `Subs` class (symbolic) | `functions.h:738` | ✅ Exists but **unrelated** to the Gruntz `SubsSet` |
| `series()` function | `series.h:716` | ✅ Exists: `series(expr, var, prec)` → `SeriesCoeffInterface` |
| `SeriesCoeffInterface::as_dict()` | `series.h:21` | ✅ Exists, returns `umap_int_basic` |
| `rewrite_as_exp()` | `basic.h:278`, `rewrite.cpp:125` | ✅ Exists, converts trig/hyperbolic → exponentials |
| `xreplace(map)` | `basic.h:193` | ✅ Exists, direct replacement without evaluation |
| `subs(map)` | `basic.h:191` | ✅ Exists |
| `expand(expr)` | `basic.h:270` | ✅ Exists |
| `Inf`, `NegInf` | `constants.h:73-74` | ✅ Exists as `RCP<const Infty>` |
| `Infty` class | `infinity.h:20` | ✅ Exists with `is_positive_infinity()`, `is_negative_infinity()` |
| `zero`, `one`, `minus_one` | `constants.h:59-61` | ✅ Exists as `RCP<const Integer>` |
| `E` (Euler's number) | `constants.h:67` | ✅ Exists as `RCP<const Constant>` |
| `Add`, `Mul`, `Pow` classes | `add.h`, `mul.h`, `pow.h` | ✅ Exists |
| `Log` class | `functions.h:529` | ✅ Exists |
| `Add::from_dict(coef, dict)` | `add.h:83` | ✅ Exists, takes `umap_basic_num&&` |
| `Add::dict_add_term(d, coef, t)` | `add.h:93` | ✅ Exists |
| `Add::get_coef()` | `add.h:142` | ✅ Exists, returns `RCP<const Number>` |
| `Add::get_dict()` | `add.h:148` | ✅ Exists, returns `const umap_basic_num&` |
| `Add::as_two_terms(a, b)` | `add.h:114` | ✅ Exists |
| `Mul::from_dict(coef, dict)` | `mul.h:96` | ✅ Exists, takes `map_basic_basic&&` |
| `Mul::dict_add_term(d, exp, base)` | `mul.h:99` | ✅ Exists (order: dict, exponent, base) |
| `Mul::get_coef()` | `mul.h:125` | ✅ Exists, returns `RCP<const Number>` |
| `Mul::get_dict()` | `mul.h:129` | ✅ Exists, returns `const map_basic_basic&` |
| `Mul::as_two_terms(a, b)` | `mul.h:113` | ✅ Exists |
| `Pow::get_base()` / `get_exp()` | `pow.h:37,42` | ✅ Exists |
| `Log::get_arg()` | via `OneArgFunction::get_arg()` at `functions.h:36` | ✅ Exists |
| `Integer::as_int()` | `integer.h:41` | ✅ Exists |
| `Rational::as_rational_class()` | `rational.h:50` | ✅ Exists, returns `const rational_class&` |
| `integer(n)` | `integer.h:203` | ✅ Exists |
| `eq(a, b)` | `basic.h:243` | ✅ Exists |
| `neg()`, `add()`, `sub()`, `mul()`, `div()` | `add.h`, `mul.h` | ✅ All exist |
| `pow()`, `exp()` | `pow.h:51,54` | ✅ Exists |
| `log()` | `functions.h:545` | ✅ Exists |
| `TransformVisitor` | `visitor.h:261` + `visitor.cpp:167` | ✅ Exists with `apply()`, `result_`, bvisits for Add/Mul/Pow/etc. |
| `outArg(ptr)` | `symengine_rcp.h:88` | ✅ Exists |
| `Ptr<T>` | `symengine_rcp.h` | ✅ Exists |
| `NotImplementedError` | `symengine_exception.h:63` | ✅ Exists |
| `SymEngineException` | `symengine_exception.h:30` | ✅ Exists |
| `SYMENGINE_ASSERT` | `symengine_assert.h` | ✅ Exists |
| `mp_lcm(result, a, b)` | `mp_class.h:247` | ✅ Exists |
| `integer_class` | `mp_class.h` | ✅ Exists (platform-dependent: `mpz_class`, `fmpz_wrapper`, or `piranha::integer`) |
| `rational_class` | `mp_class.h` | ✅ Exists |
| `get_den(rational)` / `get_num(rational)` | `mp_class.h:266,271` | ✅ Exists |

---

## 3. What Needs to Be Created (Entirely New)

| Item | Plan Section | Notes |
|---|---|---|
| `SubsSet` class | §4.1 | No equivalent in SymEngine. Must implement `insert()`, `meets()`, `union_set()`, `do_subs()`, `empty()`, `size()`. |
| `has()` helper | §4.2 | Thin wrapper around `has_symbol()`. Trivial. |
| `independent()` helper | §4.3 | Splits Add/Mul into x-dependent and x-independent parts. Not in SymEngine. |
| `powdenest_simple()` | §4.4 | Uses `TransformVisitor` subclass. Not in SymEngine. |
| `tractable_rewrite()` | §4.5 | One-liner calling `rewrite_as_exp()`. Trivial. |
| `sign()` | §4.6 | New function. Determines sign of expression as x → ∞. |
| `compare()` | §4.7 | New function. Compares growth rates. |
| `mrv()` | §4.8 | Core MRV algorithm. New function. |
| `mrv_max3()` / `mrv_max1()` | §4.9-4.10 | New functions. Combine MRV sets. |
| `moveup2()` / `moveup()` | §4.11 | New functions. Transform x → exp(x). |
| `build_expression_tree()` | §4.12 | New function. Dependency analysis for Omega. |
| `rewrite()` | §4.13 | Core rewrite step. New function. |
| `mrv_leadterm()` | §4.14 | Computes leading term via series. New function. |
| `limitinf()` | §4.15 | Limit at infinity. New function. |
| `gruntz()` | §4.16 | Main entry point. New function. |
| Memoization infrastructure | §7 | New thread-local cache. |
| Test file | §6 | New file `test_gruntz.cpp` (~200–300 lines). |
| Build integration | §5 | Modify two `CMakeLists.txt` files. |

---

## 4. Specific Issues Found in the Plan

### 4.1 Plan's `sign()` — Missing Cache and Potential Infinite Recursion

**Plan §4.6, lines 478–482:**
```cpp
// Leadterm fallback
RCP<const Basic> c0, e0;
std::tie(c0, e0) = mrv_leadterm(e, x);
return sign(c0, x);
```

`sign()` calls `mrv_leadterm()` as a fallback, which in turn can call `sign()` (through `compare()` → `limitinf()` → `sign()`). Without memoization being actually wired into the functions (the plan declares caches but never reads/writes them in the function bodies), this risks deep or infinite mutual recursion. **The implementation must ensure cache lookups/stores in `sign()`, `limitinf()`, and `mrv_leadterm()`.**

### 4.2 Plan's `sign()` — Not All Cases Handled

The plan's `sign()` handles `Mul`, `Pow` with `E` base, `Log`, and the leadterm fallback. There are several gaps:

- **`Sin` / `Cos` / `Tan` / `Sinh` / `Cosh` / `Tanh`**: These should be rare (since `tractable_rewrite` / `rewrite_as_exp` converts them) but `sign()` could be called on them before the rewrite. Need explicit handling or a fallback through `tractable_rewrite` → `mrv_leadterm`.
- **`FunctionSymbol`**: Custom functions may appear. The leadterm fallback should catch these.
- **`is_a<Pow>(*e) && !eq(*base, *E)`**: e.g., `x^2`, `(-x)^3`. The plan's code for this case (lines 450–462) only handles `E` base and `Integer` exponent. Other cases fall through to `mrv_leadterm`. This is acceptable but should be documented.

**Recommendation**: Add a `tractable_rewrite` call in the `sign()` fallback, or ensure the leadterm fallback always catches everything. Also add explicit cache read/write.

### 4.3 Plan's `compare()` — SymEngine `Infty` vs `Basic` Comparison

**Plan §4.7, line 507:**
```cpp
if (is_a<Infty>(*L)) return '>';
```

`L = limitinf(div(la, lb), x)` returns `RCP<const Basic>`. `limitinf()` returns `Inf` (which is `RCP<const Infty>`) or `NegInf` or `zero` or a constant. Checking `is_a<Infty>(*L)` is correct but should also handle `NegInf`:

```cpp
if (is_a<Infty>(*L) && down_cast<const Infty &>(*L).is_positive_infinity()) return '>';
if (is_a<Infty>(*L) && down_cast<const Infty &>(*L).is_negative_infinity()) return '<';
```

As written, `NegInf` would fall through and return `'='`, which is wrong.

### 4.4 Plan's `rewrite()` — Unsafe `rcp_static_cast<const Pow>`

**Plan §4.13, lines 764–765:**
```cpp
RCP<const Pow> g_pow = rcp_static_cast<const Pow>(g);
RCP<const Basic> g_exp = g_pow->get_exp();
```

This assumes `g` is always a `Pow`. After `moveup2` in `mrv_leadterm()`, `x` has been replaced by `exp(x)` in the Omega dict, so all keys should contain `exp(...)` subexpressions. However:

- If `moveup2` was NOT called (because `x` was not in Omega), the representative `g` should still be a `Pow` since all MRV terms are exponentials. But this depends on correct behavior of `mrv()`.
- A safer approach would be to use `as_base_exp` from SymEngine or a check before casting.

**Recommendation**: Add `SYMENGINE_ASSERT(is_a<Pow>(*g))` before the cast, or use `Mul::as_base_exp(g, outArg(exp), outArg(base))` to safely decompose.

### 4.5 Plan's `rewrite()` — `sign()` on `g_exp`

**Plan §4.13, line 767:**
```cpp
auto sig = sign(g_exp, x);
```

`g_exp` is the exponent of the representative exponential (e.g., for `exp(x)`, `g_exp = x`). `sign(x, x) = 1`. For `exp(-x)`, `g_exp = -x`, `sign(-x, x) = -1`. This is correct.

However, after `moveup2`, if the original term was `log(x)`, the initial Omega entry was `x`, which becomes `exp(x)` after `moveup2`. So `g = exp(x)`, `g_exp = x`. This is fine.

If the original Omega entry was already an exponential like `exp(x^2)`, `g = exp(x^2)`, `g_exp = x^2`. Then `sign(x^2, x) = 1`. This is correct.

But there's a subtle case: what if `g_exp` is a product like `x * 2` (from `exp(2*x)`)? `sign(2*x, x) = 2 * sign(x, x) = 2 * 1 = 1`. Wait, `sign(2*x, x)` would call `sign(2, x) = 1`, `sign(x, x) = 1`, then `mul(1, 1) = 1`. So `sig = 1`. This is correct.

But what about `sign(-x, x)`? `sign(-x, x)` → `sign(-1, x) = -1`, `sign(x, x) = 1`, `mul(-1, 1) = -1`. So `sig = -1`. Correct: `exp(-x) -> 0` as `x -> oo`, so we want `w` not `1/w`.

### 4.6 Plan's `rewrite()` — Empty Denominator Vector

**Plan §4.13, lines 775, 815:**
```cpp
std::vector<integer_class> denominators;
...
integer_class lcm_val = compute_denominators_lcm(denominators);
```

If `denominators` is empty (all `c` are integers, not rational), `compute_denominators_lcm` starts with `result = 1` and the loop never executes, returning `1`. This is correct and safe.

### 4.7 Plan's `limitinf()` — `sign(e0, x)` with Non-Integer e0

**Plan §4.15, line 918:**
```cpp
auto sig = sign(e0, x);
```

`e0` is the second element of the `mrv_leadterm` return pair, which is `integer(lead_exp)`. So `e0` is always an `Integer`, which is a `Number`. The `sign()` function handles numbers at the top. This is correct.

However, the `sign()` of an integer `0` returns `integer(0)`, so the `eq(*sig, *zero)` branch on line 933 correctly recurses on `c0`.

### 4.8 Plan's `gruntz()` — Handling `dir` Parameter

**Plan §4.16, lines 974–978:**
```cpp
if (dir == "+") {
    subst_map[z] = add(z0, div(one, z));
} else if (dir == "-") {
    subst_map[z] = sub(z0, div(one, z));
}
```

When `z0` is an `Infty`, the `dir` parameter is ignored (the code returns before reaching this block). For finite `z0`, this transforms `z -> z0 ± 1/x`. This is the standard Gruntz approach. Correct.

### 4.9 Plan's `mrv()` — Missing `Cot` / `Csc` / `Sec` / `Coth` / `Csch` / `Sech`

The plan's `mrv()` handles `Add`, `Mul`, `Pow`, `Log`, and `exp`. Other functions are handled by `tractable_rewrite()` in `limitinf()`. This is correct — `mrv()` doesn't need to handle them directly since `limitinf()` calls `tractable_rewrite()` first.

### 4.10 Plan's `independent()` — `Mul::dict_add_term` Parameter Order

**Plan §4.3, line 359:**
```cpp
Mul::dict_add_term(var_dict, kv.second, kv.first);
```

Verified against `mul.h:99`:
```cpp
static void dict_add_term(map_basic_basic &d, const RCP<const Basic> &exp,
                          const RCP<const Basic> &t);
```

Parameter order is `(dict, exp, base/t)`. The plan's call is correct: `kv.first` is the base (from `Mul::get_dict()` which maps `base → exp`), `kv.second` is the exponent. ✅

---

## 5. Additional Shortcomings

### 5.1 Series Backend Availability

SymEngine supports three series backends: flint, piranha, and generic. The `series()` function template only works if one of these is compiled. The plan should mention checking `HAVE_SYMENGINE_FLINT`, `HAVE_SYMENGINE_PIRANHA`, or falling back to the generic backend. The `series_generic.h` header provides a pure-C++ implementation that requires no external dependencies.

### 5.2 `PowDenestVisitor` — TransformVisitor Recursion

**Plan §4.4**: The `PowDenestVisitor` calls `apply(x.get_base())` which recursively transforms the base. If the base itself contains a Pow like `(a^b)^c`, the recursion correctly handles it. But the plan's visitor only overrides `bvisit(const Pow &x)`. For all other types (Add, Mul, Symbol, etc.), it falls back to `TransformVisitor::bvisit(...)`. The `TransformVisitor`'s `bvisit(const Pow &x)` (at `visitor.cpp:196`) recursively applies to base and exponent and reconstructs via `pow()` if changed. The `PowDenestVisitor`'s override replaces this behavior. This is correct.

### 5.3 `TransformVisitor` Base Class Patterns

`RewriteAsExp` in `rewrite.cpp` uses `BaseVisitor<RewriteAsExp, TransformVisitor>`. The plan's `PowDenestVisitor` uses the simpler `public TransformVisitor`. Both styles work, but the plan's style (`using TransformVisitor::bvisit`) requires bringing base class bvisits into scope. This is fine but should be tested.

**Recommendation**: Follow the same pattern as `RewriteAsExp` for consistency:
```cpp
class PowDenestVisitor : public BaseVisitor<PowDenestVisitor, TransformVisitor>
```

### 5.4 Memoization Not Wired Into Functions

The plan declares cache maps in §7 but the function implementations in §4 never check or populate them. The implementation must ensure:

- `limitinf()` checks `tls_cache.limitinf_cache` before computing.
- `sign()` checks `tls_cache.sign_cache`.
- `mrv_leadterm()` checks `tls_cache.mrv_leadterm_cache`.
- All three store their result before returning.

Without this, the algorithm will be exponentially slow on many real expressions.

### 5.5 `SubsSet::union_set()` — Potential Shared-Pointer Leaks

**Plan §4.1, lines 269–288**: The `union_set()` method copies `*this` at the start (`SubsSet res = *this;`). The `dict` and `rewrites` members are `map_basic_basic` (std::map), and the copy will increment RCP reference counts. This is fine for correctness but may be slow for large sets. Not a real concern for the expected set sizes.

### 5.6 `build_expression_tree()` — O(n²) Complexity with No Cycle Detection

**Plan §4.12**: The nested `for` loop over `Omega.dict` twice makes this O(n²) where n = |Omega|. For typical limits, n is small (2–5). This is acceptable.

### 5.7 No `LambertW` or Special Function Handling

SymEngine has `LambertW`, `Zeta`, `Gamma`, `Erf`, etc. The Gruntz algorithm doesn't natively handle these. `tractable_rewrite()` (via `rewrite_as_exp()`) only handles trig/hyperbolic functions. The plan should mention that limits with special functions may fail with `NotImplementedError` unless additional `mrv()` cases or rewrite rules are added.

### 5.8 `SubsSet::do_subs()` Uses `xreplace`, Not `subs`

**Plan §4.1, lines 295–302**: This is the correct choice. `xreplace` performs a direct syntactic replacement of dummies with their original expressions, without evaluation. Using `subs()` could incorrectly simplify expressions (e.g., `log(exp(x))` → `x`).

### 5.9 `series()` with Dummy Symbol

**Plan §4.14, line 866:**
```cpp
auto series_result = series(f_for_series, w, prec);
```

`w` is a `Dummy` created via `dummy("w")`. The `series()` function (series.h:716) takes `const RCP<const Symbol> &var`. Since `Dummy` inherits from `Symbol`, `RCP<const Dummy>` implicitly converts to `RCP<const Symbol>`. The `SeriesVisitor::bvisit(const Symbol &x)` (series_visitor.h:269) checks `x.get_name() == varname` — the Dummy's name is `"w"`, so this works.

### 5.10 `Eq` Comparison with `E`

**Plan §4.6, line 452, and elsewhere:**
```cpp
if (eq(*p.get_base(), *E))
```

`E` is `RCP<const Constant>`. `p.get_base()` returns `RCP<const Basic>`. Dereferencing both with `*` gives `const Constant &` vs `const Basic &`. `eq()` takes `const Basic &, const Basic &` (basic.h:243). `Constant` inherits from `Basic`, so this works correctly. ✅

---

## 6. Recommendations

### 6.1 Wire Memoization Into All Functions

Add cache lookups and stores to `sign()`, `limitinf()`, and `mrv_leadterm()`. Example pattern for `limitinf`:

```cpp
RCP<const Basic> limitinf(const RCP<const Basic> &e, const RCP<const Symbol> &x) {
    // check cache
    auto it = tls_cache.limitinf_cache.find(e);
    if (it != tls_cache.limitinf_cache.end()) return it->second;

    // ... existing computation ...

    // store in cache before returning
    tls_cache.limitinf_cache[e] = result;
    return result;
}
```

### 6.2 Use Consistent Visitor Pattern

Follow the existing pattern in `rewrite.cpp` for `PowDenestVisitor`:

```cpp
class PowDenestVisitor : public BaseVisitor<PowDenestVisitor, TransformVisitor>
```

### 6.3 Add Safety Checks

- `SYMENGINE_ASSERT(is_a<Pow>(*g))` in `rewrite()` before the `rcp_static_cast`.
- Handle `NegInf` correctly in `compare()`.
- Test edge cases: `exp(-x)`, `exp(x) / x`, `x * exp(-x)`, `log(x) / x`.

### 6.4 Series Backend Selection

The plan should mention selecting the series backend. The `series_generic.h` backend requires no external dependencies and is always available. Add a note to check `HAVE_SYMENGINE_FLINT` / `HAVE_SYMENGINE_PIRANHA` and fall back to the generic backend.

### 6.5 Mark Series Expansion Prec as Tunable

The plan hard-codes `prec = 6` in `mrv_leadterm()`. This is reasonable for correctness but may need to be increased for more complex expressions. Consider making it a parameter or increasing to ~8.

### 6.6 Test Cases to Include

Based on the review, the test suite should include at minimum:

1. Polynomial limits: `x → oo`, `x^2 → oo`, `1/x → 0`
2. Exponential limits: `exp(x) → oo`, `exp(-x) → 0`
3. Log limits: `log(x) → oo`, `log(x)/x → 0`
4. Mixed: `x*exp(-x) → 0`, `exp(x)/x → oo`, `(exp(x) - exp(-x)) / (exp(x) + exp(-x)) → 1`
5. Two-sided limits at finite points: `1/x` at 0 from `+` and `-`
6. Limits at `-oo`: `exp(x)` as `x → -oo`
7. Tractable rewrite: `sin(x)/x → 0` (as `x → oo`)
8. Nested exponentials: `exp(exp(x)) → oo`
9. Limits with rational exponents: `x^(1/2) → oo`
10. `gruntz()` entry point with various `z0` and `dir` values

### 6.7 Header Include Hygiene

The plan's `gruntz.h` (§3) includes only `<symengine/basic.h>`, `<symengine/symbol.h>`, `<symengine/dict.h>`. Based on the actual dependencies, the implementation file will need:
- `<symengine/add.h>` — `Add::from_dict`, `Add::dict_add_term`
- `<symengine/mul.h>` — `Mul::from_dict`, `Mul::dict_add_term`, `mul`, `div`, `neg`
- `<symengine/pow.h>` — `Pow`, `pow`, `exp`
- `<symengine/functions.h>` — `Log`, `log`
- `<symengine/constants.h>` — `zero`, `one`, `minus_one`, `E`, `Inf`, `NegInf`
- `<symengine/series.h>` — `series()`, `SeriesCoeffInterface`
- `<symengine/visitor.h>` — `has_symbol`, `TransformVisitor`
- `<symengine/symengine_exception.h>` — `NotImplementedError`
- `<symengine/symengine_assert.h>` — `SYMENGINE_ASSERT`

### 6.8 Potential Integer Overflow in Series Exponent Search

**Plan §4.14, line 870:** `int lead_exp = -1;`
The `umap_int_basic` key is `int`. For series with very high precision, the exponent could overflow `int`. This is unlikely at `prec = 6` but should be noted.

---

## 7. Conclusion

**The plan is sound and implementable.** All SymEngine dependencies it relies on exist. The approximately 800–1000 lines of new code (SubsSet + ~15 functions) plus tests (~200–300 lines) is a realistic estimate. The main risks are:

1. **Memoization not actually wired** into the function bodies — must be added at implementation time.
2. **`compare()` mishandles `NegInf`** — needs a fix.
3. **Potentially slow for complex expressions** — series expansion is the bottleneck.
4. **Series backend availability** — the generic backend should be used as a fallback.

With the issues in §4 and §5 addressed during implementation, the Gruntz algorithm can be successfully integrated into SymEngine.
