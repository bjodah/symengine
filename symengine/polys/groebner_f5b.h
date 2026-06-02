#ifndef SYMENGINE_POLYS_GROEBNER_F5B_H
#define SYMENGINE_POLYS_GROEBNER_F5B_H

#include <symengine/polys/groebner_internal.h>
#include <symengine/polys/groebner_signature.h>
#include <vector>
#include <algorithm>
#include <memory>

namespace SymEngine {

template <typename Coeff, typename Domain>
struct LabeledPoly {
    Signature   signature;
    GPoly<Coeff, Domain> polynomial;
    unsigned    number = 0;
};

template <typename Coeff, typename Domain>
struct CriticalPair {
    Signature                 sig_um_f;     // signature of um * f
    PackedMonomial            um;
    const LabeledPoly<Coeff, Domain>* f = nullptr;
    Signature                 sig_vm_g;     // signature of vm * g
    PackedMonomial            vm;
    const LabeledPoly<Coeff, Domain>* g = nullptr;
};

template <typename Coeff, typename Domain>
struct CriticalPairDesc {
    MonomialOrder order;
    explicit CriticalPairDesc(MonomialOrder ord) : order(ord) {}

    bool operator()(const CriticalPair<Coeff, Domain>& a, const CriticalPair<Coeff, Domain>& b) const {
        SignatureLess less{order};
        if (!(a.sig_um_f == b.sig_um_f)) {
            return less(b.sig_um_f, a.sig_um_f);
        }
        return less(b.sig_vm_g, a.sig_vm_g);
    }
};

template <typename Coeff, typename Domain>
LabeledPoly<Coeff, Domain> lbp_mul_term(const LabeledPoly<Coeff, Domain>& h, const PackedMonomial& t, const Coeff& c) {
    LabeledPoly<Coeff, Domain> res;
    res.signature = signature_mul(h.signature, t);
    res.polynomial = h.polynomial;
    res.polynomial.mul_monomial(t.exponents);
    res.polynomial.scale(c);
    res.number = h.number;
    return res;
}

template <typename Coeff, typename Domain>
LabeledPoly<Coeff, Domain> lbp_sub(const LabeledPoly<Coeff, Domain>& f, const LabeledPoly<Coeff, Domain>& g) {
    LabeledPoly<Coeff, Domain> res;
    SignatureLess less{f.polynomial.order};
    if (less(f.signature, g.signature)) {
        res.signature = g.signature;
    } else {
        res.signature = f.signature;
    }
    res.polynomial = f.polynomial;
    res.polynomial.sub_poly(g.polynomial);
    res.number = f.number;
    return res;
}

template <typename Coeff, typename Domain>
CriticalPair<Coeff, Domain> make_critical_pair(const LabeledPoly<Coeff, Domain>& f,
                                               const LabeledPoly<Coeff, Domain>& g,
                                               MonomialOrder order) {
    PackedMonomial lm_f = pack_monomial(f.polynomial.leading_monomial());
    PackedMonomial lm_g = pack_monomial(g.polynomial.leading_monomial());
    PackedMonomial L = packed_lcm(lm_f, lm_g);
    PackedMonomial um = packed_divide(L, lm_f);
    PackedMonomial vm = packed_divide(L, lm_g);

    Signature sig_u = signature_mul(f.signature, um);
    Signature sig_v = signature_mul(g.signature, vm);

    CriticalPair<Coeff, Domain> cp;
    if (SignatureLess{order}(sig_u, sig_v)) {
        cp.sig_um_f = sig_v;
        cp.um = vm;
        cp.f = &g;
        cp.sig_vm_g = sig_u;
        cp.vm = um;
        cp.g = &f;
    } else {
        cp.sig_um_f = sig_u;
        cp.um = um;
        cp.f = &f;
        cp.sig_vm_g = sig_v;
        cp.vm = vm;
        cp.g = &g;
    }
    return cp;
}

template <typename Coeff, typename Domain>
LabeledPoly<Coeff, Domain> s_polynomial(const CriticalPair<Coeff, Domain>& cp) {
    Domain dom = cp.f->polynomial.dom;
    Coeff inv_lc_f = dom.inv(cp.f->polynomial.leading_coeff());
    Coeff inv_lc_g = dom.inv(cp.g->polynomial.leading_coeff());

    auto term1 = lbp_mul_term(*cp.f, cp.um, inv_lc_f);
    auto term2 = lbp_mul_term(*cp.g, cp.vm, inv_lc_g);
    return lbp_sub(term1, term2);
}

template <typename Coeff, typename Domain>
bool is_rewritable_or_comparable(const Signature& s,
                                 unsigned num,
                                 const std::vector<const LabeledPoly<Coeff, Domain>*>& B,
                                 GroebnerStats& stats)
{
    for (const auto& h : B) {
        if (h->polynomial.is_zero()) continue;
        if (h->signature.generator_index > s.generator_index &&
            packed_divides(pack_monomial(h->polynomial.leading_monomial()), s.monomial)) {
            ++stats.rejected_by_syzygy;
            return true;
        }
        if (h->signature.generator_index == s.generator_index &&
            h->number > num &&
            packed_divides(h->signature.monomial, s.monomial)) {
            ++stats.rejected_by_rewritten;
            return true;
        }
    }
    return false;
}

template <typename Coeff, typename Domain>
LabeledPoly<Coeff, Domain> f5_reduce(LabeledPoly<Coeff, Domain> f,
                                     const std::vector<const LabeledPoly<Coeff, Domain>*>& B,
                                     const GroebnerOptions& options,
                                     GroebnerStats& stats)
{
    if (f.polynomial.is_zero()) return f;
    Domain dom = f.polynomial.dom;
    for (;;) {
        check_cancellation(options);
        bool reduced = false;
        for (const auto& h : B) {
            if (h->polynomial.is_zero()) continue;
            PackedMonomial lm_h = pack_monomial(h->polynomial.leading_monomial());
            PackedMonomial lm_f = pack_monomial(f.polynomial.leading_monomial());
            if (!packed_divides(lm_h, lm_f)) continue;
            PackedMonomial t = packed_divide(lm_f, lm_h);
            Signature th_sig = signature_mul(h->signature, t);
            if (!SignatureLess{options.order}(th_sig, f.signature)) continue;

            if (options.max_reduction_steps > 0
                && stats.f5_reductions >= options.max_reduction_steps) {
                throw GroebnerLimitExceeded{};
            }
            Coeff q_coeff = dom.div(f.polynomial.leading_coeff(), h->polynomial.leading_coeff());
            f = lbp_sub(f, lbp_mul_term(*h, t, q_coeff));
            ++stats.f5_reductions;
            reduced = true;
            break;
        }
        if (!reduced || f.polynomial.is_zero()) return f;
    }
}

template <typename Coeff, typename Domain>
std::vector<CriticalPair<Coeff, Domain>> build_initial_critical_pairs(const std::vector<const LabeledPoly<Coeff, Domain>*>& B, MonomialOrder order) {
    std::vector<CriticalPair<Coeff, Domain>> CP;
    for (size_t i = 0; i < B.size(); ++i) {
        for (size_t j = i + 1; j < B.size(); ++j) {
            CP.push_back(make_critical_pair(*B[i], *B[j], order));
        }
    }
    return CP;
}

template <typename Coeff, typename Domain>
std::vector<GPoly<Coeff, Domain>> f5b_impl(std::vector<GPoly<Coeff, Domain>> F,
                                           const GroebnerOptions& options,
                                           GroebnerStats& stats)
{
    unsigned reduction_steps = 0;
    F = interreduce(F, options, reduction_steps);

    std::vector<std::unique_ptr<LabeledPoly<Coeff, Domain>>> B_storage;
    std::vector<const LabeledPoly<Coeff, Domain>*> B;
    B.reserve(F.size());

    for (size_t i = 0; i < F.size(); ++i) {
        auto lp = std::make_unique<LabeledPoly<Coeff, Domain>>();
        lp->signature = signature_mul({pack_monomial(Monomial(F[i].leading_monomial().size(), 0)), static_cast<unsigned>(i + 1)}, pack_monomial(Monomial(F[i].leading_monomial().size(), 0)));
        lp->polynomial = F[i];
        lp->number = static_cast<unsigned>(i + 1);
        B_storage.push_back(std::move(lp));
        B.push_back(B_storage.back().get());
    }

    std::sort(B.begin(), B.end(), [](const LabeledPoly<Coeff, Domain>* a, const LabeledPoly<Coeff, Domain>* b) {
        return compare_monomial_less(b->polynomial.leading_monomial(),
                                    a->polynomial.leading_monomial(),
                                    a->polynomial.order);
    });

    std::vector<CriticalPair<Coeff, Domain>> CP = build_initial_critical_pairs(B, options.order);
    std::sort(CP.begin(), CP.end(), CriticalPairDesc<Coeff, Domain>{options.order});

    unsigned k = static_cast<unsigned>(B.size());

    while (!CP.empty()) {
        check_cancellation(options);
        CriticalPair<Coeff, Domain> cp = std::move(CP.back());
        CP.pop_back();
        if (options.max_s_pairs > 0
            && stats.s_pairs_processed >= options.max_s_pairs) {
            throw GroebnerLimitExceeded{};
        }
        ++stats.s_pairs_processed;

        if (is_rewritable_or_comparable(cp.sig_um_f, cp.f->number, B, stats)) continue;
        if (is_rewritable_or_comparable(cp.sig_vm_g, cp.g->number, B, stats)) continue;

        LabeledPoly<Coeff, Domain> s = s_polynomial(cp);
        LabeledPoly<Coeff, Domain> p = f5_reduce(std::move(s), B, options, stats);
        p.polynomial.make_monic();
        p.number = ++k;

        if (p.polynomial.is_zero()) {
            ++stats.reductions_to_zero;
            continue;
        }

        auto lp_new = std::make_unique<LabeledPoly<Coeff, Domain>>(std::move(p));
        const LabeledPoly<Coeff, Domain>* p_ptr = lp_new.get();
        B_storage.push_back(std::move(lp_new));

        CP.erase(std::remove_if(CP.begin(), CP.end(), [&](const CriticalPair<Coeff, Domain>& cpx) {
            std::vector<const LabeledPoly<Coeff, Domain>*> p_vec = {p_ptr};
            return is_rewritable_or_comparable(cpx.sig_um_f, cpx.f->number, p_vec, stats) ||
                   is_rewritable_or_comparable(cpx.sig_vm_g, cpx.g->number, p_vec, stats);
        }), CP.end());

        for (const auto& g : B) {
            if (g->polynomial.is_zero()) continue;
            CriticalPair<Coeff, Domain> cp_new = make_critical_pair(*p_ptr, *g, options.order);
            std::vector<const LabeledPoly<Coeff, Domain>*> p_vec = {p_ptr};
            if (is_rewritable_or_comparable(cp_new.sig_um_f, cp_new.f->number, p_vec, stats)) continue;
            if (is_rewritable_or_comparable(cp_new.sig_vm_g, cp_new.g->number, p_vec, stats)) continue;
            CP.push_back(std::move(cp_new));
        }
        std::sort(CP.begin(), CP.end(), CriticalPairDesc<Coeff, Domain>{options.order});

        B.push_back(p_ptr);
        std::sort(B.begin(), B.end(), [](const LabeledPoly<Coeff, Domain>* a, const LabeledPoly<Coeff, Domain>* b) {
            if (a->polynomial.is_zero()) return false;
            if (b->polynomial.is_zero()) return true;
            return compare_monomial_less(b->polynomial.leading_monomial(),
                                        a->polynomial.leading_monomial(),
                                        a->polynomial.order);
        });

        stats.max_basis_size = std::max(stats.max_basis_size, static_cast<unsigned>(B.size()));
    }

    std::vector<GPoly<Coeff, Domain>> G_out;
    for (const auto& h : B) {
        if (!h->polynomial.is_zero()) {
            G_out.push_back(h->polynomial);
        }
    }
    return interreduce(G_out, options, reduction_steps);
}

} // namespace SymEngine

#endif // SYMENGINE_POLYS_GROEBNER_F5B_H
