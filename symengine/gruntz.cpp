#include <symengine/gruntz.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/functions.h>
#include <symengine/constants.h>
#include <symengine/series.h>
#include <symengine/visitor.h>
#include <symengine/symengine_exception.h>
#include <symengine/symengine_assert.h>
#include <symengine/rational.h>
#include <symengine/integer.h>
#include <symengine/infinity.h>
#include <symengine/symbol.h>
#include <symengine/eval_double.h>
#include <symengine/tuple.h>
#include <symengine/mp_class.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include <utility>
#include <unordered_map>

namespace SymEngine {

// =========================================================================
// Thread-local Memoization Caches
// =========================================================================

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

bool is_a_Function(const Basic &b) {
    return dynamic_cast<const Function *>(&b) != nullptr;
}

RCP<const Basic> recreate_function(const RCP<const Basic> &e, const vec_basic &args) {
    const Function *func_ptr = dynamic_cast<const Function *>(e.get());
    if (func_ptr != nullptr) {
        if (auto one_arg = dynamic_cast<const OneArgFunction *>(func_ptr)) {
            return one_arg->create(args);
        }
        else if (auto two_arg = dynamic_cast<const TwoArgFunction *>(func_ptr)) {
            return two_arg->create(args);
        }
        else if (auto multi_arg = dynamic_cast<const MultiArgFunction *>(func_ptr)) {
            return multi_arg->create(args);
        }
    }
    throw NotImplementedError("recreate_function: not a recognized Function type");
}

} // anonymous namespace

// =========================================================================
// SubsSet Member Functions
// =========================================================================

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

bool SubsSet::meets(const SubsSet &other) const
{
    for (const auto &kv : dict) {
        if (other.dict.find(kv.first) != other.dict.end())
            return true;
    }
    return false;
}

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
            if (not res_exps.is_null()) {
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

RCP<const Basic> SubsSet::do_subs(const RCP<const Basic> &e) const
{
    map_basic_basic rev;
    for (const auto &kv : dict) {
        rev[kv.second] = kv.first;
    }
    return e->xreplace(rev);
}

bool SubsSet::empty() const
{
    return dict.empty();
}

size_t SubsSet::size() const
{
    return dict.size();
}

// =========================================================================
// Helper functions
// =========================================================================

bool has(const Basic &b, const Basic &x)
{
    return has_symbol(b, x);
}

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
            if (has_symbol(*kv.first, *x) || has_symbol(*kv.second, *x)) {
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

class PowDenestVisitor : public BaseVisitor<PowDenestVisitor, TransformVisitor> {
public:
    using TransformVisitor::bvisit;

    void bvisit(const Pow &x) {
        auto base = apply(x.get_base());
        auto expn = apply(x.get_exp());
        if (is_a<Pow>(*base)) {
            const Pow &inner = down_cast<const Pow &>(*base);
            result_ = pow(inner.get_base(), mul(inner.get_exp(), expn));
        } else {
            if (base == x.get_base() && expn == x.get_exp()) {
                result_ = x.rcp_from_this();
            } else {
                result_ = pow(base, expn);
            }
        }
    }
};

RCP<const Basic> powdenest_simple(const RCP<const Basic> &e)
{
    PowDenestVisitor visitor;
    return visitor.apply(e);
}

class RewriteHyperbolic : public BaseVisitor<RewriteHyperbolic, TransformVisitor> {
public:
    using TransformVisitor::bvisit;

    RewriteHyperbolic() : BaseVisitor<RewriteHyperbolic, TransformVisitor>() {}

    void bvisit(const Sinh &x) {
        auto arg = apply(x.get_arg());
        result_ = div(sub(exp(arg), exp(neg(arg))), integer(2));
    }

    void bvisit(const Cosh &x) {
        auto arg = apply(x.get_arg());
        result_ = div(add(exp(arg), exp(neg(arg))), integer(2));
    }

    void bvisit(const Tanh &x) {
        auto arg = apply(x.get_arg());
        result_ = div(sub(exp(arg), exp(neg(arg))), add(exp(arg), exp(neg(arg))));
    }

    void bvisit(const Coth &x) {
        auto arg = apply(x.get_arg());
        result_ = div(add(exp(arg), exp(neg(arg))), sub(exp(arg), exp(neg(arg))));
    }
};

RCP<const Basic> tractable_rewrite(const RCP<const Basic> &e)
{
    RewriteHyperbolic visitor;
    return visitor.apply(e);
}

int get_min_exponent(const RCP<const Basic> &e, const RCP<const Symbol> &w) {
    if (!has_symbol(*e, *w)) {
        return 0;
    }
    if (eq(*e, *w)) {
        return 1;
    }
    if (is_a<Pow>(*e)) {
        const Pow &p = down_cast<const Pow &>(*e);
        if (eq(*p.get_base(), *w)) {
            if (is_a<Integer>(*p.get_exp())) {
                return static_cast<int>(down_cast<const Integer &>(*p.get_exp()).as_int());
            }
            if (is_a<Rational>(*p.get_exp())) {
                double d = eval_double(*p.get_exp());
                return static_cast<int>(std::floor(d));
            }
        }
        return 0;
    }
    if (is_a<Mul>(*e)) {
        const Mul &m = down_cast<const Mul &>(*e);
        int sum_exp = 0;
        for (const auto &kv : m.get_dict()) {
            if (eq(*kv.first, *w)) {
                if (is_a<Integer>(*kv.second)) {
                    sum_exp += static_cast<int>(down_cast<const Integer &>(*kv.second).as_int());
                } else if (is_a<Rational>(*kv.second)) {
                    double d = eval_double(*kv.second);
                    sum_exp += static_cast<int>(std::floor(d));
                }
            } else {
                sum_exp += get_min_exponent(pow(kv.first, kv.second), w);
            }
        }
        return sum_exp;
    }
    if (is_a<Add>(*e)) {
        const Add &a = down_cast<const Add &>(*e);
        int min_exp = 999999;
        for (const auto &kv : a.get_dict()) {
            int term_exp = get_min_exponent(kv.first, w);
            if (term_exp < min_exp) {
                min_exp = term_exp;
            }
        }
        if (a.get_coef()->is_zero() == false) {
            if (0 < min_exp) {
                min_exp = 0;
            }
        }
        return min_exp;
    }
    return 0;
}

class LogSimplifier : public BaseVisitor<LogSimplifier, TransformVisitor> {
private:
    RCP<const Symbol> w_;
public:
    using TransformVisitor::bvisit;
    LogSimplifier(const RCP<const Symbol> &w) : w_(w) {}

    void bvisit(const Log &x) {
        auto arg = apply(x.get_arg());
        if (eq(*arg, *w_)) {
            result_ = x.rcp_from_this();
        }
        else if (is_a<Pow>(*arg)) {
            const Pow &p = down_cast<const Pow &>(*arg);
            if (eq(*p.get_base(), *w_)) {
                result_ = mul(p.get_exp(), log(w_));
            } else if (is_a<Pow>(*p.get_base()) && eq(*down_cast<const Pow&>(*p.get_base()).get_base(), *w_)) {
                const Pow &inner = down_cast<const Pow &>(*p.get_base());
                result_ = mul(mul(inner.get_exp(), p.get_exp()), log(w_));
            } else {
                result_ = log(arg);
            }
        }
        else {
            result_ = log(arg);
        }
    }
};

RCP<const Basic> simplify_logs(const RCP<const Basic> &e, const RCP<const Symbol> &w) {
    LogSimplifier visitor(w);
    return visitor.apply(e);
}

// =========================================================================
// Core algorithm functions
// =========================================================================

RCP<const Integer> sign(const RCP<const Basic> &e,
                         const RCP<const Symbol> &x)
{
    // Cache lookup
    auto it = tls_cache.sign_cache.find(e);
    if (it != tls_cache.sign_cache.end()) {
        return it->second;
    }

    RCP<const Integer> res;

    if (is_a_Number(*e)) {
        const Number &num = down_cast<const Number &>(*e);
        if (num.is_positive()) res = integer(1);
        else if (num.is_negative()) res = integer(-1);
        else res = integer(0);
    }
    else if (!has_symbol(*e, *x)) {
        try {
            double d = eval_double(*e);
            if (d > 0) res = integer(1);
            else if (d < 0) res = integer(-1);
            else res = integer(0);
        } catch (...) {
            // Fallback if eval_double fails: default to positive or throw
            res = integer(1);
        }
    }
    else if (eq(*e, *x)) {
        res = integer(1);
    }
    else if (is_a<Pow>(*e)) {
        const Pow &p = down_cast<const Pow &>(*e);
        if (eq(*p.get_base(), *E)) {
            res = integer(1); // exp(...) is positive
        } else {
            auto sb = sign(p.get_base(), x);
            if (eq(*sb, *one)) {
                res = integer(1);
            } else if (eq(*sb, *minus_one)) {
                if (is_a<Integer>(*p.get_exp())) {
                    auto exp_val = down_cast<const Integer &>(*p.get_exp()).as_int();
                    if (exp_val % 2 == 0) res = integer(1);
                    else res = integer(-1);
                } else {
                    // Non-integer exponent of a negative base is not handled
                    RCP<const Basic> c0, e0;
                    std::tie(c0, e0) = mrv_leadterm(e, x);
                    if (eq(*e, *c0)) {
                        throw std::runtime_error("Infinite recursion detected in sign()");
                    }
                    res = sign(c0, x);
                }
            } else {
                RCP<const Basic> c0, e0;
                std::tie(c0, e0) = mrv_leadterm(e, x);
                if (eq(*e, *c0)) {
                    throw std::runtime_error("Infinite recursion detected in sign()");
                }
                res = sign(c0, x);
            }
        }
    }
    else if (is_a<Mul>(*e)) {
        RCP<const Basic> a, b;
        down_cast<const Mul &>(*e).as_two_terms(outArg(a), outArg(b));
        auto sa = sign(a, x);
        if (eq(*sa, *zero)) {
            res = integer(0);
        } else {
            res = rcp_static_cast<const Integer>(mul(sa, sign(b, x)));
        }
    }
    else if (is_a<Log>(*e)) {
        RCP<const Basic> arg = down_cast<const Log &>(*e).get_arg();
        res = sign(sub(arg, one), x);
    }
    else {
        // Leadterm fallback
        RCP<const Basic> c0, e0;
        std::tie(c0, e0) = mrv_leadterm(e, x);
        if (eq(*e, *c0)) {
            throw std::runtime_error("Infinite recursion detected in sign()");
        }
        res = sign(c0, x);
    }

    tls_cache.sign_cache[e] = res;
    return res;
}

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

    char res = '=';
    if (eq(*L, *zero)) res = '<';
    else if (is_a<Infty>(*L)) {
        res = '>';
    }
    return res;
}

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
        if (const_part.get() != one.get() && const_part.get() != zero.get()) {
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
        while (is_a<Pow>(*current) && !eq(*down_cast<const Pow &>(*current).get_base(), *E)) {
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
            
            auto su = s1.union_set(s2, null).first;
            su.rewrites[e1] = exp(e2);
            return mrv_max3(s1, e1, s2, exp(e2), su, e1, x);
        } else {
            SubsSet s;
            RCP<const Basic> expr;
            std::tie(s, expr) = mrv(exp_arg, x);
            return {s, exp(expr)};
        }
    }

    // General fallback for all functions
    if (is_a_Function(*e)) {
        vec_basic args = e->get_args();
        if (!args.empty()) {
            SubsSet s_acc;
            RCP<const Basic> e_acc;
            std::tie(s_acc, e_acc) = mrv(args[0], x);
            
            vec_basic acc_args = {e_acc};
            for (size_t i = 1; i < args.size(); ++i) {
                SubsSet s_i;
                RCP<const Basic> e_i;
                std::tie(s_i, e_i) = mrv(args[i], x);
                
                vec_basic temp_vec = acc_args;
                temp_vec.push_back(e_i);
                RCP<const Basic> temp_tuple = tuple(temp_vec);
                
                std::pair<SubsSet, RCP<const Basic>> combined = mrv_max1(s_acc, s_i, temp_tuple, x);
                s_acc = combined.first;
                acc_args = combined.second->get_args();
            }
            
            RCP<const Basic> reconstructed = recreate_function(e, acc_args);
            return {s_acc, reconstructed};
        }
    }

    throw NotImplementedError("mrv: unhandled function type");
}

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

static integer_class compute_denominators_lcm(const std::vector<integer_class> &denoms) {
    integer_class result(1);
    for (const auto &d : denoms) {
        mp_lcm(result, result, d);
    }
    return result;
}

static void collect_puiseux_denominators(const RCP<const Basic> &e, const RCP<const Symbol> &w, std::vector<integer_class> &dens) {
    if (is_a<Pow>(*e)) {
        const Pow &p = down_cast<const Pow &>(*e);
        if (eq(*p.get_base(), *w)) {
            if (is_a<Rational>(*p.get_exp())) {
                const Rational &r = down_cast<const Rational &>(*p.get_exp());
                dens.push_back(get_den(r.as_rational_class()));
            }
        }
        collect_puiseux_denominators(p.get_base(), w, dens);
        collect_puiseux_denominators(p.get_exp(), w, dens);
    } else {
        for (const auto &arg : e->get_args()) {
            collect_puiseux_denominators(arg, w, dens);
        }
    }
}

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
    SYMENGINE_ASSERT(is_a<Pow>(*g));
    RCP<const Pow> g_pow = rcp_static_cast<const Pow>(g);
    RCP<const Basic> g_exp = g_pow->get_exp();

    auto sig = sign(g_exp, x);
    SYMENGINE_ASSERT(eq(*sig, *one) || eq(*sig, *minus_one));
    
    RCP<const Basic> w_val = wsym;
    if (eq(*sig, *one)) {
        w_val = pow(wsym, minus_one); // If representative goes to oo, replace with 1/w
    }

    std::vector<std::pair<RCP<const Basic>, RCP<const Basic>>> O2;
    std::vector<integer_class> denominators;

    for (const auto &item : Omega_items) {
        RCP<const Basic> f = item.first;
        RCP<const Basic> var = item.second;
        SYMENGINE_ASSERT(is_a<Pow>(*f));
        RCP<const Basic> f_exp = rcp_static_cast<const Pow>(f)->get_exp();

        RCP<const Basic> c = limitinf(div(f_exp, g_exp), x);

        if (is_a<Rational>(*c)) {
            const Rational &c_rat = down_cast<const Rational &>(*c);
            denominators.push_back(get_den(c_rat.as_rational_class()));
        }

        RCP<const Basic> arg = f_exp;
        auto rit = Omega.rewrites.find(var);
        if (rit != Omega.rewrites.end()) {
            SYMENGINE_ASSERT(is_a<Pow>(*(rit->second)));
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

std::pair<RCP<const Basic>, RCP<const Basic>>
mrv_leadterm(const RCP<const Basic> &e,
             const RCP<const Symbol> &x)
{
    // Cache lookup
    auto it = tls_cache.mrv_leadterm_cache.find(e);
    if (it != tls_cache.mrv_leadterm_cache.end()) {
        return it->second;
    }

    std::pair<RCP<const Basic>, RCP<const Basic>> res;

    if (!has_symbol(*e, *x)) {
        res = {e, zero};
        tls_cache.mrv_leadterm_cache[e] = res;
        return res;
    }

    SubsSet Omega;
    RCP<const Basic> exps;
    std::tie(Omega, exps) = mrv(e, x);

    if (Omega.empty()) {
        res = {exps, zero};
        tls_cache.mrv_leadterm_cache[e] = res;
        return res;
    }

    if (Omega.dict.find(x) != Omega.dict.end()) {
        Omega = moveup2(Omega, x);
        exps = moveup(exps, x);
    }

    auto w = dummy("w");
    RCP<const Basic> f, logw;
    std::tie(f, logw) = rewrite(exps, Omega, x, w);

    // Replace the limit variable x with its representation in terms of w, i.e., -log(w)
    map_basic_basic x_repl;
    x_repl[x] = neg(log(w));
    f = f->xreplace(x_repl);

    // Simplify log(w^c) expressions to avoid poles/log(1/w) in series expansion
    f = simplify_logs(f, w);

    // Singular Log Workaround: Substitutes log(w) with _logw_dummy
    auto logw_dummy = symbol("_logw_dummy");
    map_basic_basic log_repl;
    log_repl[log(w)] = logw_dummy;
    
    RCP<const Basic> f_for_series = f->xreplace(log_repl);

    // Find any fractional powers of w in f_for_series
    std::vector<integer_class> dens;
    collect_puiseux_denominators(f_for_series, w, dens);
    integer_class lcm_val = compute_denominators_lcm(dens);

    auto t = dummy("t");
    bool is_puiseux = (lcm_val > 1);
    if (is_puiseux) {
        map_basic_basic w_repl;
        w_repl[w] = pow(t, integer(lcm_val));
        f_for_series = f_for_series->xreplace(w_repl);
        f_for_series = powdenest_simple(f_for_series);
        f_for_series = expand(f_for_series);
    }

    auto var_for_series = is_puiseux ? t : w;

    // Find the minimum exponent of var_for_series in f_for_series
    int min_exp = get_min_exponent(f_for_series, var_for_series);
    
    if (min_exp < 0) {
        // Multiply by var_for_series^(-min_exp) to clear the pole
        f_for_series = mul(f_for_series, pow(var_for_series, integer(-min_exp)));
        f_for_series = expand(f_for_series);
    }

    unsigned int prec = 6;
    auto series_result = series(f_for_series, var_for_series, prec);
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

    // Adjust the lead exponent back if we cleared a pole
    if (min_exp < 0) {
        lead_exp += min_exp;
    }

    RCP<const Basic> lead_exp_final;
    if (is_puiseux) {
        lead_exp_final = Rational::from_two_ints(*integer(lead_exp), *integer(lcm_val));
    } else {
        lead_exp_final = integer(lead_exp);
    }

    // Restore actual logw expression into lead coefficient
    map_basic_basic log_restore;
    log_restore[logw_dummy] = logw;
    lead_coeff = lead_coeff->xreplace(log_restore);

    res = {lead_coeff, lead_exp_final};
    tls_cache.mrv_leadterm_cache[e] = res;
    return res;
}

RCP<const Basic> limitinf(const RCP<const Basic> &e,
                           const RCP<const Symbol> &x)
{
    // Cache lookup
    auto it = tls_cache.limitinf_cache.find(e);
    if (it != tls_cache.limitinf_cache.end()) {
        return it->second;
    }

    RCP<const Basic> res;

    if (!has_symbol(*e, *x)) {
        res = e;
        tls_cache.limitinf_cache[e] = res;
        return res;
    }

    RCP<const Basic> e_rewritten = tractable_rewrite(e);
    e_rewritten = powdenest_simple(e_rewritten);

    RCP<const Basic> c0, e0;
    std::tie(c0, e0) = mrv_leadterm(e_rewritten, x);

    auto sig = sign(e0, x);

    if (eq(*sig, *one)) {
        res = zero; // e0 > 0: lim = 0
    }
    else if (eq(*sig, *minus_one)) {
        auto sc = sign(c0, x);
        if (eq(*sc, *zero)) {
            throw std::runtime_error("Leading term sign should not be zero");
        }
        if (eq(*sc, *one)) res = Inf;
        else res = NegInf;
    }
    else if (eq(*sig, *zero)) {
        res = limitinf(c0, x); // e0 = 0: recurse on c0
    }
    else {
        throw std::runtime_error("Unexpected sign result");
    }

    tls_cache.limitinf_cache[e] = res;
    return res;
}

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

} // namespace SymEngine
