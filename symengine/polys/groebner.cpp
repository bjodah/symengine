#include <symengine/polys/groebner.h>
#include <symengine/polys/groebner_internal.h>
#include <symengine/rational.h>
#include <symengine/expression.h>
#include <symengine/symbol.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/integer.h>
#include <symengine/visitor.h>
#include <symengine/solve.h>
#include <symengine/tuple.h>
#include <symengine/symengine_exception.h>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <set>

namespace SymEngine {

// Overloaded coefficient converters to handle compile-time type dispatch in templates cleanly in C++11.
inline rational_class coeff_from_basic(const RCP<const Basic> &final_coeff_expr, const RationalCoeffDomain &dom) {
    if (is_a<Integer>(*final_coeff_expr)) {
        return down_cast<const Integer &>(*final_coeff_expr).as_integer_class();
    } else if (is_a<Rational>(*final_coeff_expr)) {
        return down_cast<const Rational &>(*final_coeff_expr).as_rational_class();
    } else {
        throw SymEngineException("Non-rational coefficient in rational mode: " + final_coeff_expr->__str__());
    }
}

inline Expression coeff_from_basic(const RCP<const Basic> &final_coeff_expr, const ExpressionCoeffDomain &dom) {
    return Expression(final_coeff_expr);
}

inline RCP<const Basic> coeff_to_basic(const rational_class &c, const RationalCoeffDomain &dom) {
    return Rational::from_mpq(c);
}

inline RCP<const Basic> coeff_to_basic(const Expression &expr, const ExpressionCoeffDomain &dom) {
    return expr.get_basic();
}

bool has_symbolic_parameters(const vec_basic &polys, const vec_sym &variables) {
    set_basic vars;
    for (const auto &v : variables) {
        vars.insert(v);
    }
    for (const auto &poly : polys) {
        set_basic syms = free_symbols(*poly);
        for (const auto &s : syms) {
            if (vars.find(s) == vars.end()) {
                return true;
            }
        }
    }
    return false;
}

void parse_term(const RCP<const Basic> &term_expr,
                const vec_sym &variables,
                Monomial &monomial,
                RCP<const Basic> &coeff_expr) {
    size_t n = variables.size();
    monomial.assign(n, 0);
    coeff_expr = one;

    std::vector<std::pair<RCP<const Basic>, RCP<const Basic>>> factors;
    
    if (is_a<Mul>(*term_expr)) {
        const Mul &m = down_cast<const Mul &>(*term_expr);
        coeff_expr = m.get_coef();
        for (const auto &it : m.get_dict()) {
            factors.push_back({it.first, it.second});
        }
    } else {
        if (is_a<Pow>(*term_expr)) {
            const Pow &p = down_cast<const Pow &>(*term_expr);
            factors.push_back({p.get_base(), p.get_exp()});
        } else {
            factors.push_back({term_expr, one});
        }
    }

    std::unordered_map<RCP<const Basic>, size_t, RCPBasicHash, RCPBasicKeyEq> var_indices;
    for (size_t i = 0; i < n; ++i) {
        var_indices[variables[i]] = i;
    }

    for (const auto &f : factors) {
        RCP<const Basic> base = f.first;
        RCP<const Basic> exponent = f.second;
        
        auto it = var_indices.find(base);
        if (it != var_indices.end()) {
            if (is_a<Integer>(*exponent)) {
                long ex_val = down_cast<const Integer &>(*exponent).as_int();
                if (ex_val < 0) {
                    throw SymEngineException("Non-polynomial variable exponent (negative)");
                }
                monomial[it->second] += (unsigned int)ex_val;
            } else {
                throw SymEngineException("Non-polynomial variable exponent (non-integer)");
            }
        } else {
            if (eq(*exponent, *one)) {
                coeff_expr = mul(coeff_expr, base);
            } else {
                coeff_expr = mul(coeff_expr, pow(base, exponent));
            }
        }
    }
}

template <typename Coeff, typename Domain>
GPoly<Coeff, Domain> basic_to_gpoly(const RCP<const Basic> &poly, const vec_sym &variables, MonomialOrder order) {
    Domain dom;
    GPoly<Coeff, Domain> res(order);

    RCP<const Basic> expanded = expand(poly);

    auto parse_and_add_term = [&](const RCP<const Basic> &term_expr, const RCP<const Basic> &sum_coef) {
        Monomial mon;
        RCP<const Basic> coeff_expr;
        parse_term(term_expr, variables, mon, coeff_expr);
        RCP<const Basic> final_coeff_expr = mul(sum_coef, coeff_expr);
        
        Coeff c = coeff_from_basic(final_coeff_expr, dom);
        res.terms.push_back({mon, c});
    };

    if (is_a<Add>(*expanded)) {
        const Add &add_expr = down_cast<const Add &>(*expanded);
        if (!add_expr.get_coef()->is_zero()) {
            parse_and_add_term(one, add_expr.get_coef());
        }
        for (const auto &it : add_expr.get_dict()) {
            parse_and_add_term(it.first, it.second);
        }
    } else {
        parse_and_add_term(expanded, one);
    }

    res.canonicalize();
    return res;
}

template <typename Coeff, typename Domain>
RCP<const Basic> gpoly_to_basic(const GPoly<Coeff, Domain> &gpoly, const vec_sym &variables) {
    vec_basic sum_terms;
    for (const auto &t : gpoly.terms) {
        RCP<const Basic> coeff_expr = coeff_to_basic(t.coeff, gpoly.dom);

        RCP<const Basic> mon_expr = one;
        for (size_t i = 0; i < t.monomial.size(); ++i) {
            unsigned int ex = t.monomial[i];
            if (ex == 1) {
                mon_expr = mul(mon_expr, variables[i]);
            } else if (ex > 1) {
                mon_expr = mul(mon_expr, pow(variables[i], integer(ex)));
            }
        }
        sum_terms.push_back(mul(coeff_expr, mon_expr));
    }
    if (sum_terms.empty()) {
        return integer(0);
    }
    return add(sum_terms);
}

template <typename Coeff, typename Domain>
GPoly<Coeff, Domain> normal_form(GPoly<Coeff, Domain> p, const std::vector<GPoly<Coeff, Domain>> &G, const GroebnerOptions &options, unsigned &reduction_steps) {
    GPoly<Coeff, Domain> remainder(p.order);
    Domain dom = p.dom;
    while (!p.is_zero()) {
        if (options.cancellation_token && options.cancellation_token->is_cancelled()) {
            throw SymEngineException("cancelled");
        }
        if (options.max_reduction_steps > 0 && reduction_steps >= options.max_reduction_steps) {
            throw SymEngineException("limit exceeded: max reduction steps");
        }
        
        bool reduced = false;
        const Monomial &lm_p = p.leading_monomial();
        for (const auto &g : G) {
            if (g.is_zero()) continue;
            const Monomial &lm_g = g.leading_monomial();
            if (divides(lm_g, lm_p)) {
                Monomial q_mon = quotient(lm_p, lm_g);
                Coeff q_coeff = dom.div(p.leading_coeff(), g.leading_coeff());
                
                GPoly<Coeff, Domain> term_g = g;
                term_g.mul_monomial(q_mon);
                term_g.scale(q_coeff);
                p.sub_poly(term_g);
                
                reduced = true;
                reduction_steps++;
                break;
            }
        }
        if (!reduced) {
            remainder.terms.push_back(p.terms[0]);
            p.terms.erase(p.terms.begin());
        }
    }
    remainder.canonicalize();
    return remainder;
}

template <typename Coeff, typename Domain>
std::vector<GPoly<Coeff, Domain>> interreduce(std::vector<GPoly<Coeff, Domain>> G, const GroebnerOptions &options, unsigned &reduction_steps) {
    bool changed = true;
    while (changed) {
        changed = false;
        for (size_t i = 0; i < G.size(); ++i) {
            GPoly<Coeff, Domain> g = G[i];
            G.erase(G.begin() + i);
            GPoly<Coeff, Domain> reduced_g = normal_form(g, G, options, reduction_steps);
            if (!reduced_g.is_zero()) {
                reduced_g.make_monic();
                G.insert(G.begin() + i, reduced_g);
                if (reduced_g.terms.size() != g.terms.size() || reduced_g.leading_monomial() != g.leading_monomial()) {
                    changed = true;
                }
            } else {
                changed = true;
                --i;
            }
        }
    }
    return G;
}

template <typename Coeff, typename Domain>
std::vector<GPoly<Coeff, Domain>> buchberger(std::vector<GPoly<Coeff, Domain>> F, const GroebnerOptions &options, GroebnerStats &stats) {
    Domain dom;
    unsigned s_pairs_processed = 0;
    unsigned reductions_to_zero = 0;
    unsigned max_basis_size = 0;
    unsigned reduction_steps = 0;
    
    std::vector<GPoly<Coeff, Domain>> G = F;
    for (auto &g : G) {
        g.make_monic();
    }
    if (options.interreduce_input) {
        G = interreduce(G, options, reduction_steps);
    }
    
    struct Pair {
        size_t i, j;
        Monomial lcm_mon;
        exponent_t lcm_deg;
    };
    std::vector<Pair> pairs;
    auto add_pairs_for = [&](size_t i) {
        for (size_t j = 0; j < i; ++j) {
            Monomial lm_i = G[i].leading_monomial();
            Monomial lm_j = G[j].leading_monomial();
            Monomial l = lcm(lm_i, lm_j);
            pairs.push_back({j, i, l, degree(l)});
        }
    };
    
    for (size_t i = 0; i < G.size(); ++i) {
        add_pairs_for(i);
    }
    
    auto pair_cmp = [&options](const Pair &a, const Pair &b) {
        if (a.lcm_deg != b.lcm_deg) {
            return a.lcm_deg < b.lcm_deg;
        }
        return compare_monomial_less(a.lcm_mon, b.lcm_mon, options.order);
    };
    std::sort(pairs.begin(), pairs.end(), pair_cmp);
    
    while (!pairs.empty()) {
        if (options.cancellation_token && options.cancellation_token->is_cancelled()) {
            throw SymEngineException("cancelled");
        }
        if (options.max_s_pairs > 0 && s_pairs_processed >= options.max_s_pairs) {
            throw SymEngineException("limit exceeded: max s-pairs");
        }
        
        Pair p = pairs.front();
        pairs.erase(pairs.begin());
        s_pairs_processed++;
        
        GPoly<Coeff, Domain> &f = G[p.i];
        GPoly<Coeff, Domain> &g = G[p.j];
        
        Monomial lm_f = f.leading_monomial();
        Monomial lm_g = g.leading_monomial();
        
        bool coprime = true;
        for (size_t k = 0; k < lm_f.size(); ++k) {
            if (lm_f[k] > 0 && lm_g[k] > 0) {
                coprime = false;
                break;
            }
        }
        if (coprime) {
            reductions_to_zero++;
            continue;
        }
        
        Monomial l = p.lcm_mon;
        Monomial q_f = quotient(l, lm_f);
        Monomial q_g = quotient(l, lm_g);
        
        GPoly<Coeff, Domain> term_f = f;
        term_f.mul_monomial(q_f);
        
        GPoly<Coeff, Domain> term_g = g;
        term_g.mul_monomial(q_g);
        
        GPoly<Coeff, Domain> s_poly = term_f;
        s_poly.sub_poly(term_g);
        
        GPoly<Coeff, Domain> rem = normal_form(s_poly, G, options, reduction_steps);
        if (!rem.is_zero()) {
            rem.make_monic();
            G.push_back(rem);
            max_basis_size = std::max(max_basis_size, (unsigned)G.size());
            
            add_pairs_for(G.size() - 1);
            std::sort(pairs.begin(), pairs.end(), pair_cmp);
        } else {
            reductions_to_zero++;
        }
    }
    
    if (options.reduced) {
        G = interreduce(G, options, reduction_steps);
    }
    
    stats.s_pairs_processed = s_pairs_processed;
    stats.reductions_to_zero = reductions_to_zero;
    stats.max_basis_size = std::max(max_basis_size, (unsigned)G.size());
    stats.input_polys = F.size();
    stats.output_polys = G.size();
    
    return G;
}

GroebnerResult groebner_basis(const vec_basic &polys,
                              const vec_sym &variables,
                              const GroebnerOptions &options) {
    GroebnerResult result;
    result.variables = variables;
    result.order = options.order;
    
    if (polys.empty()) {
        result.status = GroebnerStatus::Success;
        return result;
    }
    
    bool symbolic = has_symbolic_parameters(polys, variables);
    
    try {
        if (symbolic) {
            std::vector<GPoly<Expression, ExpressionCoeffDomain>> F;
            for (const auto &poly : polys) {
                F.push_back(basic_to_gpoly<Expression, ExpressionCoeffDomain>(poly, variables, options.order));
            }
            GroebnerStats stats;
            auto G = buchberger<Expression, ExpressionCoeffDomain>(F, options, stats);
            result.basis.clear();
            if (options.sort_output) {
                std::sort(G.begin(), G.end(), [options](const GPoly<Expression, ExpressionCoeffDomain> &a, const GPoly<Expression, ExpressionCoeffDomain> &b) {
                    return compare_monomial_less(b.leading_monomial(), a.leading_monomial(), options.order);
                });
            }
            for (const auto &g : G) {
                result.basis.push_back(gpoly_to_basic<Expression, ExpressionCoeffDomain>(g, variables));
            }
            result.stats = stats;
            result.status = GroebnerStatus::Success;
        } else {
            std::vector<GPoly<rational_class, RationalCoeffDomain>> F;
            for (const auto &poly : polys) {
                F.push_back(basic_to_gpoly<rational_class, RationalCoeffDomain>(poly, variables, options.order));
            }
            GroebnerStats stats;
            auto G = buchberger<rational_class, RationalCoeffDomain>(F, options, stats);
            result.basis.clear();
            if (options.sort_output) {
                std::sort(G.begin(), G.end(), [options](const GPoly<rational_class, RationalCoeffDomain> &a, const GPoly<rational_class, RationalCoeffDomain> &b) {
                    return compare_monomial_less(b.leading_monomial(), a.leading_monomial(), options.order);
                });
            }
            for (const auto &g : G) {
                result.basis.push_back(gpoly_to_basic<rational_class, RationalCoeffDomain>(g, variables));
            }
            result.stats = stats;
            result.status = GroebnerStatus::Success;
        }
    } catch (const SymEngineException &e) {
        if (std::string(e.what()).find("cancelled") != std::string::npos) {
            result.status = GroebnerStatus::Cancelled;
        } else if (std::string(e.what()).find("limit exceeded") != std::string::npos) {
            result.status = GroebnerStatus::ResourceLimitExceeded;
        } else {
            throw e;
        }
    }
    return result;
}

bool is_zero_dimensional(const std::vector<Monomial> &leading_monomials, size_t n) {
    std::vector<bool> has_pure_power(n, false);
    for (const auto &lm : leading_monomials) {
        int nonzero_var_idx = -1;
        bool is_pure = true;
        for (size_t i = 0; i < n; ++i) {
            if (lm[i] > 0) {
                if (nonzero_var_idx == -1) {
                    nonzero_var_idx = (int)i;
                } else {
                    is_pure = false;
                    break;
                }
            }
        }
        if (is_pure && nonzero_var_idx != -1) {
            has_pure_power[nonzero_var_idx] = true;
        }
    }
    for (bool b : has_pure_power) {
        if (!b) return false;
    }
    return true;
}

void find_standard_monomials(size_t var_idx, Monomial current, const std::vector<Monomial> &leading_monomials, const std::vector<unsigned int> &bounds, std::vector<Monomial> &std_monomials) {
    if (var_idx == bounds.size()) {
        bool is_divisible = false;
        for (const auto &lm : leading_monomials) {
            if (divides(lm, current)) {
                is_divisible = true;
                break;
            }
        }
        if (!is_divisible) {
            std_monomials.push_back(current);
        }
        return;
    }

    for (unsigned int e = 0; e < bounds[var_idx]; ++e) {
        current[var_idx] = e;
        find_standard_monomials(var_idx + 1, current, leading_monomials, bounds, std_monomials);
    }
}

template <typename Coeff, typename Domain>
GroebnerResult fglm_impl(const GroebnerResult &source, MonomialOrder target_order, const GroebnerOptions &options) {
    Domain dom;
    size_t n = source.variables.size();
    
    std::vector<GPoly<Coeff, Domain>> G1;
    std::vector<Monomial> leading_monomials;
    for (const auto &poly : source.basis) {
        auto gp = basic_to_gpoly<Coeff, Domain>(poly, source.variables, source.order);
        if (!gp.is_zero()) {
            G1.push_back(gp);
            leading_monomials.push_back(gp.leading_monomial());
        }
    }
    
    if (!is_zero_dimensional(leading_monomials, n)) {
        throw SymEngineException("not zero-dimensional");
    }
    
    std::vector<unsigned int> bounds(n, 0);
    for (const auto &lm : leading_monomials) {
        int nonzero_var_idx = -1;
        bool is_pure = true;
        for (size_t i = 0; i < n; ++i) {
            if (lm[i] > 0) {
                if (nonzero_var_idx == -1) {
                    nonzero_var_idx = (int)i;
                } else {
                    is_pure = false;
                    break;
                }
            }
        }
        if (is_pure && nonzero_var_idx != -1) {
            if (bounds[nonzero_var_idx] == 0 || lm[nonzero_var_idx] < bounds[nonzero_var_idx]) {
                bounds[nonzero_var_idx] = lm[nonzero_var_idx];
            }
        }
    }
    
    std::vector<Monomial> O;
    Monomial current(n, 0);
    find_standard_monomials(0, current, leading_monomials, bounds, O);
    
    std::sort(O.begin(), O.end(), [target_order](const Monomial &a, const Monomial &b) {
        return compare_monomial_less(a, b, target_order);
    });
    
    size_t d = O.size();
    
    auto get_coordinate_vector = [&](const GPoly<Coeff, Domain> &gp) {
        std::vector<Coeff> vec(d, dom.zero());
        for (const auto &t : gp.terms) {
            auto it = std::find(O.begin(), O.end(), t.monomial);
            if (it != O.end()) {
                size_t idx = std::distance(O.begin(), it);
                vec[idx] = t.coeff;
            }
        }
        return vec;
    };
    
    struct EchelonRow {
        Monomial target_monomial;
        std::vector<Coeff> reduced_vec;
        GPoly<Coeff, Domain> relation;
    };
    std::vector<EchelonRow> echelon_rows;
    std::vector<Monomial> std_basis_target;
    std::vector<GPoly<Coeff, Domain>> G2;
    
    std::vector<Monomial> queue;
    std::vector<Monomial> queue_set;
    queue.push_back(Monomial(n, 0));
    queue_set.push_back(Monomial(n, 0));
    
    unsigned reduction_steps = 0;
    
    while (!queue.empty()) {
        if (options.cancellation_token && options.cancellation_token->is_cancelled()) {
            throw SymEngineException("cancelled");
        }
        
        std::sort(queue.begin(), queue.end(), [target_order](const Monomial &a, const Monomial &b) {
            return compare_monomial_less(a, b, target_order);
        });
        
        Monomial v = queue.front();
        queue.erase(queue.begin());
        
        bool divisible_by_g2 = false;
        for (const auto &g : G2) {
            if (divides(g.leading_monomial(), v)) {
                divisible_by_g2 = true;
                break;
            }
        }
        if (divisible_by_g2) continue;
        
        GPoly<Coeff, Domain> gp_v(source.order);
        gp_v.terms.push_back({v, dom.one()});
        
        GPoly<Coeff, Domain> nf_v = normal_form(gp_v, G1, options, reduction_steps);
        std::vector<Coeff> c_v = get_coordinate_vector(nf_v);
        
        GPoly<Coeff, Domain> rel_v(target_order);
        rel_v.terms.push_back({v, dom.one()});
        
        std::vector<Coeff> red_v = c_v;
        
        for (const auto &row : echelon_rows) {
            int pivot = -1;
            for (size_t i = 0; i < d; ++i) {
                if (!dom.is_zero(row.reduced_vec[i])) {
                    pivot = (int)i;
                    break;
                }
            }
            if (pivot != -1 && !dom.is_zero(red_v[pivot])) {
                Coeff factor = dom.div(red_v[pivot], row.reduced_vec[pivot]);
                for (size_t i = 0; i < d; ++i) {
                    red_v[i] = dom.sub(red_v[i], dom.mul(factor, row.reduced_vec[i]));
                }
                auto scaled_relation = row.relation;
                scaled_relation.scale(factor);
                rel_v.sub_poly(scaled_relation);
            }
        }
        
        bool is_zero = true;
        for (const auto &coeff : red_v) {
            if (!dom.is_zero(coeff)) {
                is_zero = false;
                break;
            }
        }
        
        if (is_zero) {
            rel_v.make_monic();
            G2.push_back(rel_v);
        } else {
            int pivot = -1;
            for (size_t i = 0; i < d; ++i) {
                if (!dom.is_zero(red_v[i])) {
                    pivot = (int)i;
                    break;
                }
            }
            Coeff p_coeff = red_v[pivot];
            for (size_t i = 0; i < d; ++i) {
                red_v[i] = dom.div(red_v[i], p_coeff);
            }
            rel_v.scale(dom.div(dom.one(), p_coeff));
            
            EchelonRow new_row = {v, red_v, rel_v};
            echelon_rows.push_back(new_row);
            std_basis_target.push_back(v);
            
            for (size_t i = 0; i < n; ++i) {
                Monomial next_v = v;
                next_v[i]++;
                
                if (std::find(queue_set.begin(), queue_set.end(), next_v) == queue_set.end() &&
                    std::find(std_basis_target.begin(), std_basis_target.end(), next_v) == std_basis_target.end()) {
                    
                    bool div_g2 = false;
                    for (const auto &g : G2) {
                        if (divides(g.leading_monomial(), next_v)) {
                            div_g2 = true;
                            break;
                        }
                    }
                    if (!div_g2) {
                        queue.push_back(next_v);
                        queue_set.push_back(next_v);
                    }
                }
            }
        }
    }
    
    if (options.sort_output) {
        std::sort(G2.begin(), G2.end(), [target_order](const GPoly<Coeff, Domain> &a, const GPoly<Coeff, Domain> &b) {
            return compare_monomial_less(b.leading_monomial(), a.leading_monomial(), target_order);
        });
    }
    
    GroebnerResult res;
    res.variables = source.variables;
    res.order = target_order;
    res.status = GroebnerStatus::Success;
    res.basis.clear();
    for (const auto &g : G2) {
        res.basis.push_back(gpoly_to_basic<Coeff, Domain>(g, source.variables));
    }
    return res;
}

GroebnerResult fglm_convert(const GroebnerResult &source,
                            MonomialOrder target_order,
                            const GroebnerOptions &options) {
    GroebnerResult result;
    result.variables = source.variables;
    result.order = target_order;
    
    if (source.status != GroebnerStatus::Success) {
        result.status = source.status;
        return result;
    }
    
    if (source.basis.empty()) {
        result.status = GroebnerStatus::Success;
        return result;
    }
    
    bool symbolic = has_symbolic_parameters(source.basis, source.variables);
    
    try {
        if (symbolic) {
            result = fglm_impl<Expression, ExpressionCoeffDomain>(source, target_order, options);
        } else {
            result = fglm_impl<rational_class, RationalCoeffDomain>(source, target_order, options);
        }
    } catch (const SymEngineException &e) {
        if (std::string(e.what()).find("cancelled") != std::string::npos) {
            result.status = GroebnerStatus::Cancelled;
        } else if (std::string(e.what()).find("limit exceeded") != std::string::npos) {
            result.status = GroebnerStatus::ResourceLimitExceeded;
        } else if (std::string(e.what()).find("not zero-dimensional") != std::string::npos) {
            result.status = GroebnerStatus::NotZeroDimensional;
        } else {
            throw e;
        }
    }
    return result;
}

GroebnerResult augmented_groebner_basis(const vec_basic &polys,
                                        const vec_sym &variables,
                                        const RCP<const Symbol> &selected_variable,
                                        const RCP<const Symbol> &auxiliary_variable,
                                        const GroebnerOptions &options) {
    vec_basic augmented_polys = polys;
    augmented_polys.push_back(sub(auxiliary_variable, selected_variable));
    
    vec_sym augmented_variables = variables;
    augmented_variables.push_back(auxiliary_variable);
    
    return groebner_basis(augmented_polys, augmented_variables, options);
}

TriangularUnivariateSolution
extract_univariate_linear_shape(const GroebnerResult &lex_basis,
                                const RCP<const Symbol> &root_variable) {
    TriangularUnivariateSolution sol;
    sol.root_variable = root_variable;
    
    vec_basic eqns;
    RCP<const Basic> upoly; // defaults to null
    
    set_basic other_vars;
    for (const auto &v : lex_basis.variables) {
        if (!eq(*v, *root_variable)) {
            other_vars.insert(v);
        }
    }
    
    for (const auto &poly : lex_basis.basis) {
        set_basic syms = free_symbols(*poly);
        bool has_other = false;
        for (const auto &s : syms) {
            if (other_vars.find(s) != other_vars.end()) {
                has_other = true;
                break;
            }
        }
        if (!has_other) {
            if (!eq(*poly, *zero) && !is_a<Integer>(*poly) && !is_a<Rational>(*poly)) {
                upoly = poly;
            }
        } else {
            eqns.push_back(poly);
        }
    }
    
    if (upoly.is_null()) {
        throw SymEngineException("Univariate polynomial not found");
    }
    
    sol.univariate_polynomial = upoly;
    
    vec_sym linear_vars;
    for (const auto &v : lex_basis.variables) {
        if (!eq(*v, *root_variable)) {
            linear_vars.push_back(v);
        }
    }
    sol.linear_variables = linear_vars;
    
    try {
        auto pair_mat = linear_eqns_to_matrix(eqns, linear_vars);
        DenseMatrix A = pair_mat.first;
        DenseMatrix b = pair_mat.second;
        
        unsigned nrows = A.nrows();
        unsigned ncols = A.ncols();
        DenseMatrix augmented(nrows, ncols + 1);
        for (unsigned i = 0; i < nrows; ++i) {
            for (unsigned j = 0; j < ncols; ++j) {
                augmented.set(i, j, A.get(i, j));
            }
            augmented.set(i, ncols, b.get(i, 0));
        }
        vec_basic solution_vec = linsolve(augmented, linear_vars);
        sol.linear_solution = DenseMatrix(static_cast<unsigned int>(linear_vars.size()), 1, solution_vec);
    } catch (...) {
        throw SymEngineException("System is not linear in the remaining variables");
    }
    
    return sol;
}

RCP<const Set> solve_poly_system_via_univariate_root(const vec_basic &equations,
                                                     const vec_sym &variables,
                                                     const RCP<const Symbol> &selected_variable,
                                                     const RCP<const Symbol> &auxiliary_variable,
                                                     const GroebnerOptions &options) {
    GroebnerResult aug_gb = augmented_groebner_basis(equations, variables, selected_variable, auxiliary_variable, options);
    if (aug_gb.status != GroebnerStatus::Success) {
        return emptyset();
    }
    
    GroebnerResult lex_gb = fglm_convert(aug_gb, MonomialOrder::Lex, options);
    if (lex_gb.status != GroebnerStatus::Success) {
        return emptyset();
    }
    
    try {
        TriangularUnivariateSolution sol = extract_univariate_linear_shape(lex_gb, auxiliary_variable);
        
        RCP<const Set> root_set = solve_poly(sol.univariate_polynomial, sol.root_variable);
        if (is_a<EmptySet>(*root_set)) {
            return emptyset();
        }
        
        if (is_a<FiniteSet>(*root_set)) {
            const FiniteSet &fset = down_cast<const FiniteSet &>(*root_set);
            set_basic points;
            
            for (const auto &root : fset.get_container()) {
                map_basic_basic subs_map;
                subs_map[sol.root_variable] = root;
                
                vec_basic point_vals;
                for (const auto &v : variables) {
                    if (eq(*v, *selected_variable)) {
                        point_vals.push_back(root);
                    } else {
                        auto it = std::find(sol.linear_variables.begin(), sol.linear_variables.end(), v);
                        if (it != sol.linear_variables.end()) {
                            size_t idx = std::distance(sol.linear_variables.begin(), it);
                            RCP<const Basic> sol_expr = sol.linear_solution.get(static_cast<unsigned int>(idx), 0);
                            point_vals.push_back(sol_expr->subs(subs_map));
                        } else {
                            point_vals.push_back(v);
                        }
                    }
                }
                points.insert(tuple(point_vals));
            }
            return make_rcp<const FiniteSet>(points);
        }
    } catch (...) {
    }
    
    return emptyset();
}

RCP<const Set> solve_poly_system(const vec_basic &equations,
                                 const vec_sym &variables,
                                 const GroebnerOptions &options) {
    if (variables.empty()) {
        return emptyset();
    }
    RCP<const Symbol> selected_var = variables[0];
    RCP<const Symbol> aux_var = symbol("zaug");
    return solve_poly_system_via_univariate_root(equations, variables, selected_var, aux_var, options);
}

} // namespace SymEngine
