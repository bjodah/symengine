#ifndef SYMENGINE_POLYS_GROEBNER_M4GB_H
#define SYMENGINE_POLYS_GROEBNER_M4GB_H

#include <symengine/polys/groebner_internal.h>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <cassert>

namespace SymEngine {

using DensePoly = std::vector<uint64_t>;

struct PolyMatrixEntry {
    DensePoly tail;
    int generation = 0;
    OrdinalMonomialId basis_lm = 0;
};

struct BasisRep {
    std::unordered_set<OrdinalMonomialId> matrix_multiples;
    OrdinalMonomialId next_mul = 0;
    OrdinalMonomialId next_res = 0;
};

struct PolyMatrix {
    GFpCoeffDomain field;
    RuntimeDegRevLexCodec codec;
    MonomialOrder order = MonomialOrder::DegRevLex;
    OrdinalMonomialId upper_bound = 0;

    std::vector<OrdinalMonomialId> dense_index;
    std::unordered_map<OrdinalMonomialId, size_t> dense_invindex;
    std::unordered_map<OrdinalMonomialId, PolyMatrixEntry> matrix;

    std::map<OrdinalMonomialId, unsigned> basis;
    std::unordered_map<OrdinalMonomialId, BasisRep> basis_info;

    int generation = 0;
    GroebnerOptions options;

    PolyMatrix() : codec(0, 0) {}
};

struct M4GBCriticalPair {
    OrdinalMonomialId lcm = 0;
    OrdinalMonomialId lm_f = 0;
    OrdinalMonomialId lm_g = 0;

    bool operator<(const M4GBCriticalPair& other) const {
        return lcm > other.lcm;
    }
};

struct M4GBContext {
    PolyMatrix pm;
    std::vector<M4GBCriticalPair> critical_pairs;
    std::vector<std::pair<OrdinalMonomialId, std::map<OrdinalMonomialId, uint64_t>>> to_add;
    GroebnerOptions options;
    GroebnerStats stats;
};

// Forward declarations of helpers
std::map<OrdinalMonomialId, uint64_t> to_sparse(const DensePoly& tail, const PolyMatrix& pm);
void reduce_sparse(std::map<OrdinalMonomialId, uint64_t>& poly, PolyMatrix& pm);
DensePoly to_dense(const std::map<OrdinalMonomialId, uint64_t>& sparse, const std::unordered_map<OrdinalMonomialId, size_t>& new_dense_invindex);
DensePoly get_u_g(OrdinalMonomialId u, const DensePoly& g_tail, PolyMatrix& pm);
void shrink(PolyMatrix& pm);
void assert_invariants(const PolyMatrix& pm);

inline std::map<OrdinalMonomialId, uint64_t> to_sparse(const DensePoly& tail, const PolyMatrix& pm) {
    std::map<OrdinalMonomialId, uint64_t> res;
    for (size_t i = 0; i < tail.size(); ++i) {
        if (!pm.field.is_zero(tail[i])) {
            res[pm.dense_index[i]] = tail[i];
        }
    }
    return res;
}

inline void reduce_sparse(std::map<OrdinalMonomialId, uint64_t>& poly, PolyMatrix& pm) {
    bool reduced = true;
    while (reduced) {
        reduced = false;
        OrdinalMonomialId m_to_reduce = 0;
        OrdinalMonomialId divisor_lm = 0;
        uint64_t coeff = 0;
        for (auto it = poly.rbegin(); it != poly.rend(); ++it) {
            OrdinalMonomialId m = it->first;
            for (const auto& pair : pm.basis) {
                if (pm.codec.divides(pair.first, m)) {
                    m_to_reduce = m;
                    divisor_lm = pair.first;
                    coeff = it->second;
                    reduced = true;
                    break;
                }
            }
            if (reduced) {
                break;
            }
        }
        if (reduced) {
            poly.erase(m_to_reduce);
            OrdinalMonomialId u = pm.codec.divide(m_to_reduce, divisor_lm);
            const DensePoly& basis_tail = pm.matrix[divisor_lm].tail;
            DensePoly u_g_dense = get_u_g(u, basis_tail, pm);
            std::map<OrdinalMonomialId, uint64_t> u_g_sparse = to_sparse(u_g_dense, pm);
            for (const auto& term : u_g_sparse) {
                uint64_t prod = pm.field.mul(coeff, term.second);
                uint64_t new_val = pm.field.sub(poly[term.first], prod);
                if (pm.field.is_zero(new_val)) {
                    poly.erase(term.first);
                } else {
                    poly[term.first] = new_val;
                }
            }
        }
    }
}


inline DensePoly to_dense(const std::map<OrdinalMonomialId, uint64_t>& sparse, const std::unordered_map<OrdinalMonomialId, size_t>& new_dense_invindex) {
    DensePoly res;
    for (const auto& term : sparse) {
        auto it = new_dense_invindex.find(term.first);
        if (it != new_dense_invindex.end()) {
            if (it->second >= res.size()) res.resize(it->second + 1, 0);
            res[it->second] = term.second;
        }
    }
    return res;
}

inline DensePoly get_u_g(OrdinalMonomialId u, const DensePoly& g_tail, PolyMatrix& pm) {
    DensePoly ret;
    if (g_tail.empty()) return ret;

    OrdinalMonomialId top = pm.codec.multiply(u, pm.dense_index[g_tail.size() - 1]);
    auto ub_it = std::lower_bound(pm.dense_index.begin(), pm.dense_index.end(), top);
    ret.reserve(ub_it - pm.dense_index.begin() + 1);

    for (size_t i = g_tail.size(); i-- > 0; ) {
        if (pm.field.is_zero(g_tail[i])) continue;
        OrdinalMonomialId m = pm.codec.multiply(u, pm.dense_index[i]);

        auto dit = pm.dense_invindex.find(m);
        if (dit != pm.dense_invindex.end()) {
            if (dit->second >= ret.size()) ret.resize(dit->second + 1, 0);
            ret[dit->second] = pm.field.add(ret[dit->second], g_tail[i]);
        } else {
            auto mit = pm.matrix.find(m);
            if (mit == pm.matrix.end())
                throw std::runtime_error("M4GB invariant violated: m neither dense nor reducible");
            const DensePoly& reducer = mit->second.tail;
            for (size_t j = 0; j < reducer.size(); ++j) {
                if (pm.field.is_zero(reducer[j])) continue;
                uint64_t prod = pm.field.mul(g_tail[i], reducer[j]);
                if (j >= ret.size()) ret.resize(j + 1, 0);
                ret[j] = pm.field.sub(ret[j], prod);
            }
        }

        if (should_check_cancellation(static_cast<unsigned>(i + 1), pm.options))
            check_cancellation(pm.options);
    }
    while (!ret.empty() && pm.field.is_zero(ret.back())) {
        ret.pop_back();
    }
    return ret;
}

inline void shrink(PolyMatrix& pm) {
    std::vector<OrdinalMonomialId> reducible_ids;
    for (auto id : pm.dense_index) {
        for (const auto& pair : pm.basis) {
            if (pm.codec.divides(pair.first, id)) {
                reducible_ids.push_back(id);
                break;
            }
        }
    }
    if (reducible_ids.empty()) return;

    for (auto id : reducible_ids) {
        if (pm.matrix.find(id) != pm.matrix.end()) continue;
        OrdinalMonomialId divisor_lm = 0;
        for (const auto& pair : pm.basis) {
            if (pm.codec.divides(pair.first, id)) {
                divisor_lm = pair.first;
                break;
            }
        }
        PolyMatrixEntry entry;
        entry.basis_lm = divisor_lm;
        entry.generation = pm.generation;
        OrdinalMonomialId u = pm.codec.divide(id, divisor_lm);
        const DensePoly& basis_tail = pm.matrix[divisor_lm].tail;
        entry.tail = get_u_g(u, basis_tail, pm);
        pm.matrix[id] = entry;
    }

    std::vector<OrdinalMonomialId> new_dense_index;
    for (auto id : pm.dense_index) {
        bool reducible = false;
        for (const auto& pair : pm.basis) {
            if (pm.codec.divides(pair.first, id)) {
                reducible = true;
                break;
            }
        }
        if (!reducible) {
            new_dense_index.push_back(id);
        }
    }
    std::unordered_map<OrdinalMonomialId, size_t> new_dense_invindex;
    for (size_t i = 0; i < new_dense_index.size(); ++i) {
        new_dense_invindex[new_dense_index[i]] = i;
    }

    for (auto& pair : pm.matrix) {
        auto sparse = to_sparse(pair.second.tail, pm);
        reduce_sparse(sparse, pm);
        pair.second.tail = to_dense(sparse, new_dense_invindex);
    }

    pm.dense_index = std::move(new_dense_index);
    pm.dense_invindex = std::move(new_dense_invindex);
}

inline void increase_upper_bound(OrdinalMonomialId m, PolyMatrix& pm) {
    if (m < pm.upper_bound) return;
    for (OrdinalMonomialId id = pm.upper_bound; id <= m; ++id) {
        bool divisible = false;
        OrdinalMonomialId divisor_lm = 0;
        for (const auto& pair : pm.basis) {
            if (pm.codec.divides(pair.first, id)) {
                divisible = true;
                divisor_lm = pair.first;
                break;
            }
        }
        if (divisible) {
            PolyMatrixEntry entry;
            entry.basis_lm = divisor_lm;
            entry.generation = pm.generation;
            OrdinalMonomialId u = pm.codec.divide(id, divisor_lm);
            const DensePoly& basis_tail = pm.matrix[divisor_lm].tail;
            entry.tail = get_u_g(u, basis_tail, pm);
            pm.matrix[id] = entry;
        } else {
            size_t new_idx = pm.dense_index.size();
            pm.dense_index.push_back(id);
            pm.dense_invindex[id] = new_idx;
        }
    }
    pm.upper_bound = m + 1;
}

inline void insert_basis_polynomial(OrdinalMonomialId lm, const DensePoly& tail, PolyMatrix& pm) {
    PolyMatrixEntry entry;
    entry.basis_lm = lm;
    entry.generation = pm.generation;
    entry.tail = tail;
    pm.matrix[lm] = entry;

    pm.basis[lm] = pm.codec.degree(lm);
    pm.basis_info[lm] = BasisRep();

    shrink(pm);
}

inline void assert_invariants(const PolyMatrix& pm) {
#ifndef NDEBUG
    for (auto id : pm.dense_index) {
        assert(id < pm.upper_bound);
        for (const auto& pair : pm.basis) {
            assert(!pm.codec.divides(pair.first, id));
        }
    }
    for (const auto& pair : pm.matrix) {
        assert(pair.first < pm.upper_bound);
        bool divisible = false;
        for (const auto& b : pm.basis) {
            if (pm.codec.divides(b.first, pair.first)) {
                divisible = true;
                break;
            }
        }
        assert(divisible);
    }
#endif
}

inline void selection(M4GBContext& ctx) {
    if (ctx.critical_pairs.empty()) return;

    std::vector<M4GBCriticalPair> batch;
    unsigned min_deg = ctx.pm.codec.degree(ctx.critical_pairs.back().lcm);
    
    unsigned batch_size = ctx.options.matrix_batch_size > 0 ? ctx.options.matrix_batch_size : 512;
    while (!ctx.critical_pairs.empty() && batch.size() < batch_size) {
        if (ctx.pm.codec.degree(ctx.critical_pairs.back().lcm) == min_deg) {
            batch.push_back(std::move(ctx.critical_pairs.back()));
            ctx.critical_pairs.pop_back();
        } else {
            break;
        }
    }

    OrdinalMonomialId max_lcm = 0;
    for (const auto& cp : batch) {
        if (cp.lcm > max_lcm) max_lcm = cp.lcm;
    }
    increase_upper_bound(max_lcm, ctx.pm);

    std::vector<MatrixRow<uint64_t>> submatrix;
    for (const auto& cp : batch) {
        OrdinalMonomialId u_f = ctx.pm.codec.divide(cp.lcm, cp.lm_f);
        OrdinalMonomialId u_g = ctx.pm.codec.divide(cp.lcm, cp.lm_g);

        DensePoly row_f = get_u_g(u_f, ctx.pm.matrix[cp.lm_f].tail, ctx.pm);
        DensePoly row_g = get_u_g(u_g, ctx.pm.matrix[cp.lm_g].tail, ctx.pm);

        MatrixRow<uint64_t> row;
        std::map<OrdinalMonomialId, uint64_t> diff_sparse;
        std::map<OrdinalMonomialId, uint64_t> sparse_f = to_sparse(row_f, ctx.pm);
        std::map<OrdinalMonomialId, uint64_t> sparse_g = to_sparse(row_g, ctx.pm);
        for (const auto& term : sparse_f) {
            diff_sparse[term.first] = term.second;
        }
        for (const auto& term : sparse_g) {
            diff_sparse[term.first] = ctx.pm.field.sub(diff_sparse[term.first], term.second);
            if (ctx.pm.field.is_zero(diff_sparse[term.first])) {
                diff_sparse.erase(term.first);
            }
        }
        for (const auto& term : diff_sparse) {
            row.columns.push_back(term.first);
            row.values.push_back(term.second);
        }
        if (!row.columns.empty()) {
            submatrix.push_back(std::move(row));
        }
    }

    if (submatrix.empty()) return;

    GroebnerMatrixReducer<uint64_t, GFpCoeffDomain> reducer;
    reducer.reduce(submatrix, ctx.options, ctx.stats, ctx.pm.field);

    for (const auto& row : submatrix) {
        if (row.columns.empty()) continue;
        OrdinalMonomialId lm = row.columns[0];
        std::map<OrdinalMonomialId, uint64_t> tail_sparse;
        for (size_t i = 1; i < row.columns.size(); ++i) {
            tail_sparse[row.columns[i]] = row.values[i];
        }
        ctx.to_add.push_back({lm, std::move(tail_sparse)});
    }
}

inline void update(M4GBContext& ctx) {
    std::sort(ctx.to_add.begin(), ctx.to_add.end(), [](const std::pair<OrdinalMonomialId, std::map<OrdinalMonomialId, uint64_t>>& a, const std::pair<OrdinalMonomialId, std::map<OrdinalMonomialId, uint64_t>>& b) {
        return a.first < b.first;
    });

    for (size_t i = 0; i < ctx.to_add.size(); ++i) {
        OrdinalMonomialId lm = ctx.to_add[i].first;
        std::map<OrdinalMonomialId, uint64_t> tail_sparse = ctx.to_add[i].second;

        bool reducible = false;
        for (const auto& pair : ctx.pm.basis) {
            if (ctx.pm.codec.divides(pair.first, lm)) {
                reducible = true;
                break;
            }
        }
        if (reducible) {
            tail_sparse[lm] = 1;
            reduce_sparse(tail_sparse, ctx.pm);
            if (tail_sparse.empty()) continue;
            
            OrdinalMonomialId new_lm = tail_sparse.rbegin()->first;
            uint64_t lc = tail_sparse.rbegin()->second;
            uint64_t inv_lc = ctx.pm.field.inv(lc);
            std::map<OrdinalMonomialId, uint64_t> new_tail_sparse;
            for (const auto& term : tail_sparse) {
                if (term.first != new_lm) {
                    new_tail_sparse[term.first] = ctx.pm.field.mul(term.second, inv_lc);
                }
            }
            lm = new_lm;
            tail_sparse = std::move(new_tail_sparse);
        }

        DensePoly tail = to_dense(tail_sparse, ctx.pm.dense_invindex);

        insert_basis_polynomial(lm, tail, ctx.pm);

        std::vector<M4GBCriticalPair> new_pairs;
        for (const auto& pair : ctx.pm.basis) {
            if (pair.first == lm) continue;
            OrdinalMonomialId L = ctx.pm.codec.lcm(lm, pair.first);
            
            std::vector<exponent_t> exps_lm, exps_g;
            ctx.pm.codec.from_index(lm, exps_lm);
            ctx.pm.codec.from_index(pair.first, exps_g);
            bool coprime = true;
            for (size_t k = 0; k < exps_lm.size(); ++k) {
                if (exps_lm[k] > 0 && exps_g[k] > 0) {
                    coprime = false;
                    break;
                }
            }
            if (coprime) continue;

            M4GBCriticalPair cp;
            cp.lcm = L;
            cp.lm_f = lm;
            cp.lm_g = pair.first;
            new_pairs.push_back(cp);
        }

        ctx.critical_pairs.erase(std::remove_if(ctx.critical_pairs.begin(), ctx.critical_pairs.end(), [&](const M4GBCriticalPair& cp) {
            return ctx.pm.codec.divides(lm, cp.lcm) &&
                   !(cp.lm_f == lm) && !(cp.lm_g == lm) &&
                   ctx.pm.codec.lcm(cp.lm_f, lm) < cp.lcm &&
                   ctx.pm.codec.lcm(cp.lm_g, lm) < cp.lcm;
        }), ctx.critical_pairs.end());

        for (auto& cp : new_pairs) {
            ctx.critical_pairs.push_back(std::move(cp));
        }
        std::sort(ctx.critical_pairs.begin(), ctx.critical_pairs.end());

        std::vector<OrdinalMonomialId> to_remove;
        for (const auto& pair : ctx.pm.basis) {
            if (pair.first != lm && ctx.pm.codec.divides(lm, pair.first)) {
                to_remove.push_back(pair.first);
            }
        }
        for (auto id : to_remove) {
            ctx.pm.basis.erase(id);
            ctx.pm.basis_info.erase(id);
            
            DensePoly old_tail = ctx.pm.matrix[id].tail;
            PolyMatrixEntry entry;
            entry.basis_lm = lm;
            entry.generation = ctx.pm.generation;
            OrdinalMonomialId u = ctx.pm.codec.divide(id, lm);
            entry.tail = get_u_g(u, ctx.pm.matrix[lm].tail, ctx.pm);
            ctx.pm.matrix[id] = entry;
            
            auto sparse = to_sparse(old_tail, ctx.pm);
            sparse[id] = 1;
            reduce_sparse(sparse, ctx.pm);
            if (sparse.empty()) continue;
            
            OrdinalMonomialId new_lm = sparse.rbegin()->first;
            uint64_t lc = sparse.rbegin()->second;
            uint64_t inv_lc = ctx.pm.field.inv(lc);
            std::map<OrdinalMonomialId, uint64_t> new_tail_sparse;
            for (const auto& term : sparse) {
                if (term.first != new_lm) {
                    new_tail_sparse[term.first] = ctx.pm.field.mul(term.second, inv_lc);
                }
            }
            
            ctx.to_add.push_back({new_lm, std::move(new_tail_sparse)});
        }
    }
    ctx.to_add.clear();
}

inline void initialize_from_input(M4GBContext& ctx, const std::vector<GPoly<uint64_t, GFpCoeffDomain>>& F) {
    unsigned max_deg = 0;
    for (const auto& f : F) {
        if (f.is_zero()) continue;
        unsigned d = 0;
        for (auto e : f.leading_monomial()) {
            d += e;
        }
        if (d > max_deg) max_deg = d;
    }
    unsigned max_codec_deg = max_deg + 30;
    
    size_t num_vars = F.empty() ? 0 : F[0].leading_monomial().size();
    ctx.pm.codec = RuntimeDegRevLexCodec(static_cast<unsigned>(num_vars),
                                         max_codec_deg);

    for (const auto& f : F) {
        if (f.is_zero()) continue;
        GPoly<uint64_t, GFpCoeffDomain> f_monic = f;
        f_monic.make_monic();
        
        OrdinalMonomialId lm = ctx.pm.codec.to_index(f_monic.leading_monomial());
        
        std::map<OrdinalMonomialId, uint64_t> tail_sparse;
        for (size_t i = 1; i < f_monic.terms.size(); ++i) {
            OrdinalMonomialId id = ctx.pm.codec.to_index(f_monic.terms[i].monomial);
            tail_sparse[id] = f_monic.terms[i].coeff;
        }
        
        increase_upper_bound(lm, ctx.pm);
        for (const auto& term : tail_sparse) {
            increase_upper_bound(term.first, ctx.pm);
        }
        
        ctx.to_add.push_back({lm, std::move(tail_sparse)});
    }
}

inline void extract_basis(M4GBContext& ctx, std::vector<GPoly<uint64_t, GFpCoeffDomain>>& out) {
    out.clear();
    for (const auto& pair : ctx.pm.basis) {
        OrdinalMonomialId lm = pair.first;
        const DensePoly& tail = ctx.pm.matrix[lm].tail;
        
        GPoly<uint64_t, GFpCoeffDomain> g(ctx.pm.order);
        g.dom = ctx.pm.field;
        
        std::vector<exponent_t> lm_exps;
        ctx.pm.codec.from_index(lm, lm_exps);
        g.terms.push_back({lm_exps, 1});
        
        std::map<OrdinalMonomialId, uint64_t> sparse = to_sparse(tail, ctx.pm);
        for (auto it = sparse.rbegin(); it != sparse.rend(); ++it) {
            std::vector<exponent_t> exps;
            ctx.pm.codec.from_index(it->first, exps);
            g.terms.push_back({exps, it->second});
        }
        out.push_back(std::move(g));
    }
}

inline std::vector<GPoly<uint64_t, GFpCoeffDomain>> groebner_basis_m4gb_impl(
    std::vector<GPoly<uint64_t, GFpCoeffDomain>> F,
    const GroebnerOptions& options,
    GroebnerStats& stats)
{
    M4GBContext ctx;
    ctx.pm.field = F[0].dom;
    ctx.pm.order = options.order;
    ctx.pm.options = options;
    ctx.options = options;
    ctx.stats = stats;

    initialize_from_input(ctx, F);
    update(ctx);

    while (!ctx.critical_pairs.empty() || !ctx.to_add.empty()) {
        check_cancellation(options);
        if (ctx.critical_pairs.empty() && ctx.to_add.empty()) break;
        selection(ctx);
        update(ctx);
    }

    std::vector<GPoly<uint64_t, GFpCoeffDomain>> G_out;
    extract_basis(ctx, G_out);
    stats = ctx.stats;
    stats.polymatrix_entries_built = static_cast<unsigned>(ctx.pm.matrix.size());
    stats.polymatrix_generation = static_cast<unsigned>(ctx.pm.generation);
    unsigned reduction_steps = 0;
    return interreduce(G_out, options, reduction_steps);
}

} // namespace SymEngine

#endif // SYMENGINE_POLYS_GROEBNER_M4GB_H
