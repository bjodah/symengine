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
