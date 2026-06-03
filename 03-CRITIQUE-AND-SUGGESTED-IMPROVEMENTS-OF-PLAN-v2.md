# Critique and Suggested Improvements of Gruntz Implementation Plan (v2)

This document provides an in-depth critical review of the proposed plan [02-PLAN-GRUNTZ-v2.md](file:///work/02-PLAN-GRUNTZ-v2.md) for implementing the Gruntz limit algorithm in SymEngine. The critique identifies major algorithmic bugs, type-mismatch compilation errors, API inconsistencies with the existing SymEngine codebase, and proposes robust design enhancements.

---

## 1. Algorithmic Errors

### 1.1 Unordered Series Coefficient Map (`umap_int_basic`)
*   **The Issue in the Plan:**  
    Section 4.14, Step 4 assumes the series dictionary returned by `as_dict()` is sorted:
    ```cpp
    lead_exp = series_dict.begin()->first     // lowest exponent
    lead_coeff = series_dict.begin()->second  // coefficient
    ```
    And states: *"Call `as_dict()` to get the `umap_int_basic` (which is `std::map<int, RCP<const Basic>>` sorted by exponent)."*
*   **The Reality:**  
    In SymEngine, `umap_int_basic` is defined as:
    ```cpp
    typedef std::unordered_map<int, RCP<const Basic>> umap_int_basic;
    ```
    Because it is an **unordered** map, `series_dict.begin()` points to an arbitrary element, not the leading term. This will cause non-deterministic behavior and incorrect limit calculations.
*   **Proposed Fix:**  
    Iterate through the dictionary to find the minimum key (exponent) with a non-zero coefficient:
    ```cpp
    int lead_exp = -1;
    RCP<const Basic> lead_coeff;
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
        // Handle identically zero expression at current precision
        lead_exp = 0;
        lead_coeff = zero;
    }
    ```

### 1.2 Missing Support for Logarithmic Expansion Variables
*   **The Issue in the Plan:**  
    In Gruntz, when we rewrite an expression in terms of a dummy variable $w \to 0$, log terms (e.g., $\log(x)$ where $x = 1/w$) are transformed into expressions containing $\log(w)$. The series expansion is then computed as $w \to 0$.
*   **The Reality:**  
    SymEngine's `series()` expansion does not support expanding expressions containing $\log(w)$ in terms of $w$, because $\log(0)$ is singular. Unlike SymPy, which has a `logx` parameter to treat $\log(w)$ as a separate dummy variable during expansion, SymEngine's `series` will throw an error or crash.
*   **Proposed Fix:**  
    Before calling `series()`, temporarily substitute all instances of $\log(w)$ with a fresh dummy symbol (e.g., `_logw_dummy`). Since `_logw_dummy` does not depend on $w$, the series engine will treat it as a constant coefficient:
    ```cpp
    auto logw_dummy = symbol("_logw_dummy");
    map_basic_basic log_repl;
    log_repl[log(w)] = logw_dummy;
    
    // Replace log(w) with dummy before expansion
    RCP<const Basic> f_for_series = f->xreplace(log_repl);
    auto series_result = series(f_for_series, w, prec);
    
    // Extract lead coefficient and restore actual logw expression
    map_basic_basic log_restore;
    log_restore[logw_dummy] = logw; // logw computed from rewrite
    lead_coeff = lead_coeff->xreplace(log_restore);
    ```

---

## 2. Compilation and Type Errors

### 2.1 Type Mismatch in `rewrite()` for `wsym`
*   **The Issue in the Plan:**  
    Section 4.13, Step 2 states:
    ```cpp
    if sig == 1:
        wsym = pow(wsym, integer(-1))   // wsym = 1/wsym (we need w → 0)
    ```
    However, the signature of `rewrite` declares `wsym` as `const RCP<const Symbol> &wsym`. A `Pow` expression (`1/wsym`) cannot be assigned to a variable of type `RCP<const Symbol>`. This will fail to compile.
*   **Proposed Fix:**  
    Keep `wsym` as a `Symbol` for final substitutions, but introduce an intermediate `RCP<const Basic> w_val` representing its value in the rewritten terms:
    ```cpp
    RCP<const Basic> w_val = wsym;
    if (eq(*sig, *one)) {
        w_val = pow(wsym, minus_one);
    }
    // ...
    // In loop:
    rewritten = mul(exp(sub(arg, mul(c, g_exp))), pow(w_val, c));
    ```

### 2.2 Invalid `Ptr` Signature in `independent()`
*   **The Issue in the Plan:**  
    The signature proposed in `gruntz.h` is:
    ```cpp
    RCP<const Basic> independent(const RCP<const Basic> &e,
                                 const RCP<const Symbol> &x,
                                 Ptr<RCP<const Basic>> &dependent);
    ```
    In `gruntz.cpp`, it is called as:
    ```cpp
    const_part = independent(e, x, outArg(dep_part));
    ```
    `outArg(dep_part)` returns a temporary `Ptr` (an rvalue). In C++, a temporary rvalue cannot bind to a non-const lvalue reference (`Ptr<...>&`). This will result in a compiler error.
*   **Proposed Fix:**  
    Change the signature to match standard SymEngine idioms:
    ```cpp
    RCP<const Basic> independent(const RCP<const Basic> &e,
                                 const RCP<const Symbol> &x,
                                 const Ptr<RCP<const Basic>> &dependent);
    ```

### 2.3 Wrong Parameter Order for `Mul::dict_add_term`
*   **The Issue in the Plan:**  
    Section 4.3 uses `Mul::dict_add_term` as:
    ```cpp
    Mul::dict_add_term(const_factors, factor, exp);
    ```
*   **The Reality:**  
    The signature of `Mul::dict_add_term` in `/work/symengine/mul.h` is:
    ```cpp
    static void dict_add_term(map_basic_basic &d, const RCP<const Basic> &exp,
                              const RCP<const Basic> &t);
    ```
    Note that the exponent (`exp`) is the *second* parameter, and the term/base (`t`) is the *third* parameter. The plan has them reversed, which would construct mathematically incorrect terms.
*   **Proposed Fix:**  
    Ensure the correct argument ordering:
    ```cpp
    Mul::dict_add_term(const_factors, exp, factor);
    ```

### 2.4 Lack of assumptions in SymEngine's `Dummy`
*   **The Issue in the Plan:**  
    Section 4.14, Step 3 creates a dummy variable using:
    ```cpp
    w = dummy("w", positive=true)
    ```
*   **The Reality:**  
    SymEngine's `dummy()` function does not accept an assumptions parameter like SymPy. The only signatures are:
    ```cpp
    inline RCP<const Dummy> dummy();
    inline RCP<const Dummy> dummy(const std::string &name);
    inline RCP<const Dummy> dummy(const std::string &name, size_t dummy_index);
    ```
*   **Proposed Fix:**  
    Change to:
    ```cpp
    auto w = dummy("w");
    ```

---

## 3. Design and Performance Enhancements

### 3.1 Using `TransformVisitor` for `powdenest_simple`
*   **The Issue in the Plan:**  
    Section 4.4 proposes a manual tree-walking algorithm for `powdenest_simple()`. Reconstructing expressions manually for all basic types (like MultiArgFunctions, etc.) is complex, verbose, and fragile.
*   **Proposed Fix:**  
    Inherit from SymEngine's built-in `TransformVisitor`. We only need to override `bvisit(const Pow &x)` to collapse nested powers, and let `TransformVisitor` handle the generic recursion automatically:
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
    
    RCP<const Basic> powdenest_simple(const RCP<const Basic> &e) {
        PowDenestVisitor visitor;
        return visitor.apply(e);
    }
    ```

### 3.2 Memoization Key Simplification and Safety
*   **The Issue in the Plan:**  
    Section 7 recommends using a manual `size_t` combined hash as the key in cache maps to combine `e` and `x`. This is prone to hash collisions.
*   **Proposed Fix:**  
    Since the symbol $x$ is invariant during all recursive steps of a single top-level `gruntz` limit invocation, the cache keys do not need to contain $x$. Using just the expression `e` (`RCP<const Basic>`) is sufficient and simplifies the cache to:
    ```cpp
    std::unordered_map<RCP<const Basic>, RCP<const Basic>, RCPBasicHash, RCPBasicKeyEq> limitinf_cache;
    std::unordered_map<RCP<const Basic>, RCP<const Integer>, RCPBasicHash, RCPBasicKeyEq> sign_cache;
    std::unordered_map<RCP<const Basic>, std::pair<RCP<const Basic>, RCP<const Basic>>, RCPBasicHash, RCPBasicKeyEq> mrv_leadterm_cache;
    ```
    This uses `RCPBasicKeyEq` which checks actual equality (`eq(*a, *b)`) instead of hashes, guaranteeing 100% safety against collisions.

### 3.3 Strict Type Check Warning (`is_a<Symbol>`)
*   **The Issue in the Plan:**  
    Section 4.16, Step 1 uses `is_a<Symbol>(z)` to validate the limit variable. However, `is_a<T>` checks strict type equality. Since `Dummy` has type ID `SYMENGINE_DUMMY` instead of `SYMENGINE_SYMBOL`, `is_a<Symbol>(z)` will return `false` if `z` is a `Dummy`, rejecting valid dummy variables.
*   **Proposed Fix:**  
    Use the helper `is_a_Symbol` from `symbol.h`, which correctly handles both `Symbol` and `Dummy` classes:
    ```cpp
    if (!is_a_Symbol(*z)) {
        throw NotImplementedError("Second argument must be a Symbol");
    }
    ```

### 3.4 LCM Helper Type
*   **The Issue in the Plan:**  
    Section 4.13, Step 6 uses `exponent = lcm(all denominators, default=1)`. In SymEngine, `ilcm` is not defined for vectors of numbers.
*   **Proposed Fix:**  
    Implement a helper to fold `mp_lcm` (from `mp_wrapper.h` / `mp_class.h`) over the denominators:
    ```cpp
    integer_class compute_denominators_lcm(const std::vector<integer_class> &denoms) {
        integer_class result(1);
        for (const auto &d : denoms) {
            result = mp_lcm(result, d);
        }
        return result;
    }
    ```

---

## 4. Summary of Suggested File Modifications

### `gruntz.h`
*   Correct `independent` signature: `Ptr<RCP<const Basic>> &dependent` $\to$ `const Ptr<RCP<const Basic>> &dependent`.
*   Correct `powdenest_simple` helper declaration.

### `gruntz.cpp`
*   Implement `PowDenestVisitor` using `TransformVisitor`.
*   Update `independent` logic for correct `Add` and `Mul` canonical construction and parameters ordering.
*   Update `rewrite` variable types (splitting `wsym` into `wsym` and `w_val`).
*   Fix `mrv_leadterm`'s `series_dict` sorting bug by iterating to find the minimum exponent, and add the `_logw_dummy` substitution workaround.
*   Fix `limitinf` variable validation check.
*   Update the thread-local cache map types to key off `RCP<const Basic>` with `RCPBasicHash` and `RCPBasicKeyEq`.
