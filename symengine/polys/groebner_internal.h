#ifndef SYMENGINE_POLYS_GROEBNER_INTERNAL_H
#define SYMENGINE_POLYS_GROEBNER_INTERNAL_H

#include <symengine/polys/groebner.h>
#include <symengine/rational.h>
#include <symengine/expression.h>
#include <symengine/symbol.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/integer.h>
#include <symengine/symengine_exception.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>

namespace SymEngine {

using exponent_t = unsigned int;
using Monomial = std::vector<exponent_t>;

inline exponent_t degree(const Monomial &m) {
    exponent_t deg = 0;
    for (auto e : m) {
        deg += e;
    }
    return deg;
}

inline bool divides(const Monomial &a, const Monomial &b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] > b[i]) return false;
    }
    return true;
}

inline Monomial quotient(const Monomial &b, const Monomial &a) {
    Monomial res(b.size());
    for (size_t i = 0; i < b.size(); ++i) {
        res[i] = b[i] - a[i];
    }
    return res;
}

inline Monomial multiply(const Monomial &a, const Monomial &b) {
    Monomial res(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        res[i] = a[i] + b[i];
    }
    return res;
}

inline Monomial lcm(const Monomial &a, const Monomial &b) {
    Monomial res(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        res[i] = std::max(a[i], b[i]);
    }
    return res;
}

inline bool compare_monomial_less(const Monomial &a, const Monomial &b, MonomialOrder order) {
    if (order == MonomialOrder::Lex) {
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i] != b[i]) {
                return a[i] < b[i];
            }
        }
        return false;
    } else if (order == MonomialOrder::GrLex) {
        exponent_t deg_a = degree(a);
        exponent_t deg_b = degree(b);
        if (deg_a != deg_b) {
            return deg_a < deg_b;
        }
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i] != b[i]) {
                return a[i] < b[i];
            }
        }
        return false;
    } else { // MonomialOrder::DegRevLex
        exponent_t deg_a = degree(a);
        exponent_t deg_b = degree(b);
        if (deg_a != deg_b) {
            return deg_a < deg_b;
        }
        for (size_t i = a.size(); i > 0; --i) {
            size_t idx = i - 1;
            if (a[idx] != b[idx]) {
                return a[idx] > b[idx];
            }
        }
        return false;
    }
}

template <typename Coeff>
struct Term {
    Monomial monomial;
    Coeff coeff;
};

struct RationalCoeffDomain {
    using Coeff = rational_class;
    
    Coeff zero() const { return 0; }
    Coeff one() const { return 1; }
    bool is_zero(const Coeff &c) const { return c == 0; }
    bool is_one(const Coeff &c) const { return c == 1; }
    Coeff neg(const Coeff &c) const { return -c; }
    Coeff add(const Coeff &a, const Coeff &b) const { return a + b; }
    Coeff sub(const Coeff &a, const Coeff &b) const { return a - b; }
    Coeff mul(const Coeff &a, const Coeff &b) const { return a * b; }
    Coeff div(const Coeff &a, const Coeff &b) const { return a / b; }
    Coeff normalize(const Coeff &c) const { return c; }
};

struct ExpressionCoeffDomain {
    using Coeff = Expression;

    Coeff zero() const { return Expression(SymEngine::zero); }
    Coeff one() const { return Expression(SymEngine::one); }
    bool is_zero(const Coeff &c) const { return eq(*c.get_basic(), *SymEngine::zero); }
    bool is_one(const Coeff &c) const { return eq(*c.get_basic(), *SymEngine::one); }
    Coeff neg(const Coeff &c) const { return Expression(SymEngine::mul(integer(-1), c.get_basic())); }
    Coeff add(const Coeff &a, const Coeff &b) const {
        return expand(Expression(SymEngine::add(a.get_basic(), b.get_basic())));
    }
    Coeff sub(const Coeff &a, const Coeff &b) const {
        return expand(Expression(SymEngine::sub(a.get_basic(), b.get_basic())));
    }
    Coeff mul(const Coeff &a, const Coeff &b) const {
        return expand(Expression(SymEngine::mul(a.get_basic(), b.get_basic())));
    }
    Coeff div(const Coeff &a, const Coeff &b) const {
        return expand(Expression(SymEngine::div(a.get_basic(), b.get_basic())));
    }
    Coeff normalize(const Coeff &c) const { return c; }
};

template <typename Coeff, typename Domain>
class GPoly {
public:
    std::vector<Term<Coeff>> terms;
    Domain dom;
    MonomialOrder order;

    GPoly() : order(MonomialOrder::DegRevLex) {}
    GPoly(MonomialOrder ord) : order(ord) {}
    GPoly(const std::vector<Term<Coeff>> &t, MonomialOrder ord) : terms(t), order(ord) {
        canonicalize();
    }

    void canonicalize() {
        std::vector<Term<Coeff>> cleaned;
        for (auto &t : terms) {
            if (!dom.is_zero(t.coeff)) {
                cleaned.push_back(t);
            }
        }
        
        MonomialOrder ord = order;
        std::sort(cleaned.begin(), cleaned.end(), [ord](const Term<Coeff> &a, const Term<Coeff> &b) {
            return compare_monomial_less(b.monomial, a.monomial, ord);
        });

        terms.clear();
        for (auto &t : cleaned) {
            if (!terms.empty() && terms.back().monomial == t.monomial) {
                terms.back().coeff = dom.add(terms.back().coeff, t.coeff);
                if (dom.is_zero(terms.back().coeff)) {
                    terms.pop_back();
                }
            } else {
                terms.push_back(t);
            }
        }
    }

    bool is_zero() const {
        return terms.empty();
    }

    const Monomial &leading_monomial() const {
        SYMENGINE_ASSERT(!is_zero());
        return terms[0].monomial;
    }

    const Coeff &leading_coeff() const {
        SYMENGINE_ASSERT(!is_zero());
        return terms[0].coeff;
    }

    void make_monic() {
        if (is_zero()) return;
        Coeff lc = leading_coeff();
        if (dom.is_one(lc)) return;
        for (auto &t : terms) {
            t.coeff = dom.div(t.coeff, lc);
        }
        canonicalize();
    }

    void scale(const Coeff &c) {
        if (dom.is_zero(c)) {
            terms.clear();
            return;
        }
        for (auto &t : terms) {
            t.coeff = dom.mul(t.coeff, c);
        }
        canonicalize();
    }

    void mul_monomial(const Monomial &m) {
        for (auto &t : terms) {
            t.monomial = multiply(t.monomial, m);
        }
    }

    void add_poly(const GPoly &other) {
        SYMENGINE_ASSERT(order == other.order);
        terms.insert(terms.end(), other.terms.begin(), other.terms.end());
        canonicalize();
    }

    void sub_poly(const GPoly &other) {
        SYMENGINE_ASSERT(order == other.order);
        for (const auto &t : other.terms) {
            terms.push_back({t.monomial, dom.neg(t.coeff)});
        }
        canonicalize();
    }
};

} // namespace SymEngine

#endif // SYMENGINE_POLYS_GROEBNER_INTERNAL_H
