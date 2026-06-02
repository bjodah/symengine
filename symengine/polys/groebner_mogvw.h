#ifndef SYMENGINE_POLYS_GROEBNER_MOGVW_H
#define SYMENGINE_POLYS_GROEBNER_MOGVW_H

#include <symengine/polys/groebner_internal.h>
#include <symengine/polys/groebner_signature.h>
#include <vector>
#include <unordered_map>
#include <set>
#include <algorithm>

namespace SymEngine {

template <typename Coeff, typename Domain>
struct LabeledPolynomial {
    Signature signature;
    GPoly<Coeff, Domain> polynomial;
    unsigned id = 0;
};

template <typename Coeff, typename Domain>
struct LabeledMonomial {
    PackedMonomial m;                    // the "available" monomial
    unsigned       labeled_poly_id = 0;  // index into ctx.labeled_polys
    Signature      signature;            // (m / lm(f)) * lm(u)
    bool           is_syzygy = false;    // polynomial component reduces to 0
    bool           lifted = false;       // x_i * this has been queued for each i
};

template <typename Coeff, typename Domain>
struct MoGVWContext {
    std::vector<LabeledPolynomial<Coeff, Domain>> labeled_polys;
    std::vector<LabeledMonomial<Coeff, Domain>>   G;                  // the full set
    std::unordered_map<PackedMonomial, size_t, PackedMonomialHash>
        best_by_monomial;                              // m -> index into G
    std::set<Signature, SignatureLess>
        principal_syzygies;                            // sig where polynomial = 0
    LeadingMonomialIndex lm_index;                     // for divisor / cover tests
    GroebnerOptions options;
    GroebnerStats   stats;

    explicit MoGVWContext(const GroebnerOptions& opt) : options(opt) {}
};

template <typename Coeff, typename Domain>
bool is_rejected_by_lcm(const LabeledMonomial<Coeff, Domain>& brm, const MoGVWContext<Coeff, Domain>& ctx) {
    if (brm.is_syzygy) return false;
    PackedMonomial lm_f = packed_divide(brm.m, brm.signature.monomial);
    for (const auto& g : ctx.G) {
        if (g.is_syzygy) continue;
        if (g.signature.generator_index == brm.signature.generator_index &&
            g.signature.monomial.total_degree == 0) {
            continue;
        }
        if (packed_divides(g.m, brm.m)) {
            PackedMonomial lm_g = packed_divide(g.m, g.signature.monomial);
            PackedMonomial L = packed_lcm(lm_f, lm_g);
            if (!(brm.m == L)) {
                return true;
            }
        }
    }
    return false;
}

template <typename Coeff, typename Domain>
bool is_rejected_by_syzygy(const LabeledMonomial<Coeff, Domain>& brm, const MoGVWContext<Coeff, Domain>& ctx) {
    for (const auto& sig : ctx.principal_syzygies) {
        if (sig.generator_index < brm.signature.generator_index &&
            packed_divides(sig.monomial, brm.signature.monomial)) {
            return true;
        }
    }
    return false;
}

template <typename Coeff, typename Domain>
bool is_rejected_by_rewritten(const LabeledMonomial<Coeff, Domain>& brm, const MoGVWContext<Coeff, Domain>& ctx) {
    if (brm.is_syzygy) return false;
    for (const auto& g : ctx.G) {
        if (g.is_syzygy) continue;
        if (g.signature.generator_index == brm.signature.generator_index &&
            packed_divides(g.signature.monomial, brm.signature.monomial)) {
            PackedMonomial t = packed_divide(brm.signature.monomial, g.signature.monomial);
            PackedMonomial tg_m = packed_multiply(g.m, t);
            if (!(tg_m == brm.m)) {
                return true;
            }
        }
    }
    return false;
}

template <typename Coeff, typename Domain>
LabeledMonomial<Coeff, Domain> top_reduce(LabeledMonomial<Coeff, Domain> brm, MoGVWContext<Coeff, Domain>& ctx) {
    if (brm.is_syzygy) return brm;
    
    GPoly<Coeff, Domain> p = ctx.labeled_polys[brm.labeled_poly_id].polynomial;
    Domain dom = p.dom;
    if (p.is_zero()) {
        brm.is_syzygy = true;
        return brm;
    }

    SignatureLess less{p.order};
    bool reduced = true;
    while (reduced && !p.is_zero()) {
        reduced = false;
        PackedMonomial lm_p = pack_monomial(p.leading_monomial());
        
        for (const auto& g : ctx.G) {
            if (g.is_syzygy) continue;
            
            const auto& g_poly = ctx.labeled_polys[g.labeled_poly_id].polynomial;
            PackedMonomial lm_g = pack_monomial(g_poly.leading_monomial());
            
            if (packed_divides(lm_g, lm_p)) {
                PackedMonomial t = packed_divide(lm_p, lm_g);
                Signature th_sig = signature_mul(g.signature, t);
                
                if (less(th_sig, brm.signature)) {
                    Coeff q_coeff = dom.div(p.leading_coeff(), g_poly.leading_coeff());
                    GPoly<Coeff, Domain> term_g = g_poly;
                    term_g.mul_monomial(t.exponents);
                    term_g.scale(q_coeff);
                    p.sub_poly(term_g);
                    reduced = true;
                    break;
                }
            }
        }
    }

    LabeledPolynomial<Coeff, Domain> new_lp;
    new_lp.signature = brm.signature;
    new_lp.polynomial = p;
    new_lp.id = static_cast<unsigned>(ctx.labeled_polys.size());
    ctx.labeled_polys.push_back(new_lp);

    brm.labeled_poly_id = new_lp.id;
    if (p.is_zero()) {
        brm.is_syzygy = true;
    } else {
        brm.m = pack_monomial(p.leading_monomial());
    }
    return brm;
}

template <typename Coeff, typename Domain>
void mutual_reduce(LabeledMonomial<Coeff, Domain> brm, MoGVWContext<Coeff, Domain>& ctx) {
    check_cancellation(ctx.options);

    bool reducible = ctx.lm_index.has_divisor(brm.m);
    if (reducible) {
        if (is_rejected_by_lcm(brm, ctx))       { ++ctx.stats.rejected_by_lcm;       return; }
        if (is_rejected_by_syzygy(brm, ctx))    { ++ctx.stats.rejected_by_syzygy;    return; }
        if (is_rejected_by_rewritten(brm, ctx)) { ++ctx.stats.rejected_by_rewritten; return; }
    }

    LabeledMonomial<Coeff, Domain> brm_pp = top_reduce(brm, ctx);

    if (!brm_pp.is_syzygy) {
        auto it = ctx.best_by_monomial.find(brm_pp.m);
        if (it != ctx.best_by_monomial.end()) {
            LabeledMonomial<Coeff, Domain>& existing = ctx.G[it->second];
            SignatureLess less{brm_pp.m.num_vars > 0 ? ctx.options.order : MonomialOrder::DegRevLex};
            if (less(brm_pp.signature, existing.signature)) {
                LabeledMonomial<Coeff, Domain> displaced = existing;
                existing = brm_pp;
                ++ctx.stats.collisions_resolved;
                mutual_reduce(displaced, ctx);
            }
            return;
        }
    }

    size_t new_id = ctx.G.size();
    ctx.G.push_back(brm_pp);
    if (brm_pp.is_syzygy) {
        ctx.principal_syzygies.insert(brm_pp.signature);
    } else {
        ctx.best_by_monomial[brm_pp.m] = new_id;
        ctx.lm_index.insert(brm_pp.m, new_id);
    }
}

template <typename Coeff, typename Domain>
unsigned max_primitive_degree(const MoGVWContext<Coeff, Domain>& ctx) {
    unsigned max_deg = 0;
    for (const auto& g : ctx.G) {
        if (!g.is_syzygy) {
            max_deg = std::max(max_deg, g.m.total_degree);
        }
    }
    return max_deg;
}

template <typename Coeff, typename Domain>
unsigned max_cp_lcm_degree(const MoGVWContext<Coeff, Domain>& ctx) {
    // In scalar MoGVW, since critical pairs are selected degree by degree,
    // we can return maximum degree + 1 as standard bounds or loop limit.
    return max_primitive_degree(ctx) + 1;
}

template <typename Coeff, typename Domain>
void initialize_G_from_input(MoGVWContext<Coeff, Domain>& ctx, const std::vector<GPoly<Coeff, Domain>>& F) {
    for (size_t i = 0; i < F.size(); ++i) {
        LabeledPolynomial<Coeff, Domain> lp;
        lp.signature = signature_mul({pack_monomial(Monomial(F[i].leading_monomial().size(), 0)), static_cast<unsigned>(i + 1)}, pack_monomial(Monomial(F[i].leading_monomial().size(), 0)));
        lp.polynomial = F[i];
        lp.polynomial.make_monic();
        lp.id = static_cast<unsigned>(ctx.labeled_polys.size());
        ctx.labeled_polys.push_back(lp);

        LabeledMonomial<Coeff, Domain> brm;
        brm.m = pack_monomial(lp.polynomial.leading_monomial());
        brm.labeled_poly_id = lp.id;
        brm.signature = lp.signature;
        brm.is_syzygy = false;
        brm.lifted = false;

        mutual_reduce(brm, ctx);
    }
}

template <typename Coeff, typename Domain>
std::vector<GPoly<Coeff, Domain>> groebner_basis_mogvw_impl(std::vector<GPoly<Coeff, Domain>> F,
                                                             const GroebnerOptions& options,
                                                             GroebnerStats& stats) {
    unsigned reduction_steps = 0;
    F = interreduce(F, options, reduction_steps);

    MoGVWContext<Coeff, Domain> ctx(options);
    ctx.stats = stats;
    initialize_G_from_input(ctx, F);
    unsigned liftdeg = max_primitive_degree(ctx);

    size_t num_vars = F.empty() ? 0 : F[0].leading_monomial().size();

    for (;;) {
        check_cancellation(options);
        bool progressed = false;
        size_t current_size = ctx.G.size();
        for (size_t i = 0; i < current_size; ++i) {
            LabeledMonomial<Coeff, Domain> brm = ctx.G[i];
            if (brm.is_syzygy || brm.lifted) continue;
            if (brm.m.total_degree > liftdeg) continue;

            for (unsigned v = 0; v < num_vars; ++v) {
                // Queue multiplication by variable v:
                std::vector<exponent_t> exps(num_vars, 0);
                exps[v] = 1;
                PackedMonomial var_m = pack_monomial(exps);
                
                LabeledPolynomial<Coeff, Domain> lp_lifted;
                lp_lifted.signature = signature_mul(brm.signature, var_m);
                lp_lifted.polynomial = ctx.labeled_polys[brm.labeled_poly_id].polynomial;
                lp_lifted.polynomial.mul_monomial(exps);
                lp_lifted.id = static_cast<unsigned>(ctx.labeled_polys.size());
                ctx.labeled_polys.push_back(lp_lifted);

                LabeledMonomial<Coeff, Domain> lifted;
                lifted.m = packed_multiply(brm.m, var_m);
                lifted.signature = lp_lifted.signature;
                lifted.labeled_poly_id = lp_lifted.id;
                lifted.is_syzygy = false;
                lifted.lifted = false;

                mutual_reduce(lifted, ctx);
                ++ctx.stats.labeled_monomials_lifted;
            }
            ctx.G[i].lifted = true;
            progressed = true;
            liftdeg = std::max(liftdeg, max_primitive_degree(ctx));
        }

        unsigned maxcpdeg = max_cp_lcm_degree(ctx);
        if (maxcpdeg > liftdeg + 1) {
            liftdeg = maxcpdeg - 1;
            continue;
        }
        if (!progressed) break;
    }

    std::vector<GPoly<Coeff, Domain>> G_out;
    for (const auto& g : ctx.G) {
        if (!g.is_syzygy) {
            G_out.push_back(ctx.labeled_polys[g.labeled_poly_id].polynomial);
        }
    }
    stats = ctx.stats;
    return interreduce(G_out, options, reduction_steps);
}

} // namespace SymEngine

#endif // SYMENGINE_POLYS_GROEBNER_MOGVW_H
