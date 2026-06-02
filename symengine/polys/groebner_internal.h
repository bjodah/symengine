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
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <chrono>

namespace SymEngine {

struct GroebnerCancelled {};
struct GroebnerLimitExceeded {};
struct GroebnerNotZeroDimensional {};
struct GroebnerUnsupportedDomain {};

extern thread_local std::chrono::steady_clock::time_point active_groebner_start;

inline void check_cancellation(const GroebnerOptions &options) {
    if (options.cancellation_token && options.cancellation_token->is_cancelled()) {
        throw GroebnerCancelled{};
    }
    if (options.max_milliseconds > 0) {
        auto now = std::chrono::steady_clock::now();
        auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(now - active_groebner_start).count();
        if (elapsed_us >= static_cast<int64_t>(options.max_milliseconds) * 1000) {
            throw GroebnerLimitExceeded{};
        }
    }
}

inline void enforce_max_degree(unsigned degree, const GroebnerOptions &options) {
    if (options.max_degree > 0 && degree > options.max_degree) {
        throw GroebnerLimitExceeded{};
    }
}

inline bool should_check_cancellation(unsigned counter,
                                      const GroebnerOptions &options) {
    unsigned interval = options.cancellation_check_interval;
    return interval != 0 && (counter % interval) == 0;
}

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

struct PackedMonomial {
    uint64_t bits = 0;       // bit-packed exponents
    uint32_t total_degree = 0;
    uint32_t num_vars = 0;   // for the comparator + decoding
    std::vector<exponent_t> exponents;

    bool operator==(const PackedMonomial& other) const {
        if (num_vars != other.num_vars) return false;
        if (total_degree != other.total_degree) return false;
        if (bits != other.bits) return false;
        return exponents == other.exponents;
    }
};

struct PackedMonomialHash {
    size_t operator()(const PackedMonomial& m) const {
        size_t h = m.bits ^ (static_cast<size_t>(m.total_degree) << 32);
        for (auto e : m.exponents) {
            h ^= std::hash<exponent_t>{}(e) + 0x9e3779b9 + (h << 6) + (h >> 2);
        }
        return h;
    }
};

inline PackedMonomial pack_monomial(const std::vector<exponent_t>& exps) {
    PackedMonomial pm;
    pm.exponents = exps;
    pm.num_vars = static_cast<uint32_t>(exps.size());
    pm.total_degree = 0;
    for (auto e : exps) {
        pm.total_degree += e;
    }
    if (pm.num_vars <= 8) {
        bool fits = true;
        uint64_t b = 0;
        for (size_t i = 0; i < exps.size(); ++i) {
            if (exps[i] > 255) {
                fits = false;
                break;
            }
            b |= (static_cast<uint64_t>(exps[i]) << (i * 8));
        }
        if (fits) {
            pm.bits = b;
        }
    } else if (pm.num_vars <= 16) {
        bool fits = true;
        uint64_t b = 0;
        for (size_t i = 0; i < exps.size(); ++i) {
            if (exps[i] > 15) {
                fits = false;
                break;
            }
            b |= (static_cast<uint64_t>(exps[i]) << (i * 4));
        }
        if (fits) {
            pm.bits = b;
        }
    }
    return pm;
}

struct PackedMonomialLess {
    MonomialOrder order;
    explicit PackedMonomialLess(MonomialOrder ord = MonomialOrder::DegRevLex) : order(ord) {}

    bool operator()(const PackedMonomial& a, const PackedMonomial& b) const {
        return compare_monomial_less(a.exponents, b.exponents, order);
    }
};

inline bool packed_divides(const PackedMonomial& a, const PackedMonomial& b) {
    if (a.num_vars != b.num_vars) return false;
    for (size_t i = 0; i < a.exponents.size(); ++i) {
        if (a.exponents[i] > b.exponents[i]) return false;
    }
    return true;
}

inline PackedMonomial packed_divide(const PackedMonomial& b, const PackedMonomial& a) {
    std::vector<exponent_t> res(b.num_vars);
    for (size_t i = 0; i < b.num_vars; ++i) {
        res[i] = b.exponents[i] - a.exponents[i];
    }
    return pack_monomial(res);
}

inline PackedMonomial packed_multiply(const PackedMonomial& a, const PackedMonomial& b) {
    std::vector<exponent_t> res(a.num_vars);
    for (size_t i = 0; i < a.num_vars; ++i) {
        res[i] = a.exponents[i] + b.exponents[i];
    }
    return pack_monomial(res);
}

inline PackedMonomial packed_lcm(const PackedMonomial& a, const PackedMonomial& b) {
    std::vector<exponent_t> res(a.num_vars);
    for (size_t i = 0; i < a.num_vars; ++i) {
        res[i] = std::max(a.exponents[i], b.exponents[i]);
    }
    return pack_monomial(res);
}

class LeadingMonomialIndex {
    std::vector<std::pair<PackedMonomial, size_t>> entries;
public:
    void insert(const PackedMonomial& m, size_t val) {
        entries.push_back({m, val});
    }
    void clear() {
        entries.clear();
    }
    bool has_divisor(const PackedMonomial& m) const {
        for (const auto& entry : entries) {
            if (packed_divides(entry.first, m)) {
                return true;
            }
        }
        return false;
    }
};

class GFpCoeffDomain {
    uint64_t p_;
public:
    using Coeff = uint64_t;

    GFpCoeffDomain() : p_(0) {}
    explicit GFpCoeffDomain(uint64_t p) : p_(p) {}

    uint64_t modulus() const { return p_; }
    uint64_t zero() const { return 0; }
    uint64_t one() const { return 1; }
    bool is_zero(const uint64_t &c) const { return c == 0; }
    bool is_one(const uint64_t &c) const { return c == 1; }
    uint64_t neg(const uint64_t &c) const { return c == 0 ? 0 : p_ - c; }
    uint64_t add(const uint64_t &a, const uint64_t &b) const {
        uint64_t s = a + b;
        return s >= p_ ? s - p_ : s;
    }
    uint64_t sub(const uint64_t &a, const uint64_t &b) const {
        return a >= b ? a - b : a + p_ - b;
    }
    uint64_t mul(const uint64_t &a, const uint64_t &b) const {
        return (a * b) % p_;
    }
    uint64_t inv(const uint64_t &a) const {
        int64_t t = 0, newt = 1;
        int64_t r = p_, newr = a;
        while (newr != 0) {
            int64_t quotient = r / newr;
            int64_t temp = t - quotient * newt;
            t = newt;
            newt = temp;
            temp = r - quotient * newr;
            r = newr;
            newr = temp;
        }
        if (r > 1) {
            throw std::runtime_error("modular inverse does not exist");
        }
        if (t < 0) {
            t += p_;
        }
        return static_cast<uint64_t>(t);
    }
    uint64_t div(const uint64_t &a, const uint64_t &b) const {
        return mul(a, inv(b));
    }
    uint64_t normalize(const uint64_t &c) const {
        return c;
    }
};

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
    Coeff inv(const Coeff &c) const { return div(one(), c); }
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
    Coeff inv(const Coeff &c) const { return div(one(), c); }
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

#include <symengine/polys/groebner_signature.h>

namespace SymEngine {
template <typename Coeff>
struct MatrixRow {
    std::vector<uint64_t> columns;  // monomial ids
    std::vector<Coeff>    values;   // parallel to columns
    Signature             signature; // ignored when the matrix is not signature-driven
};

template <typename Coeff, typename Domain>
class GroebnerMatrixReducer {
public:
    struct Result {
        unsigned rows_in = 0;
        unsigned rows_out = 0;
        unsigned rows_zero = 0;
    };

    Result reduce(std::vector<MatrixRow<Coeff>>& rows,
                  const GroebnerOptions& options,
                  GroebnerStats& stats,
                  const Domain& dom) {
        Result res;
        res.rows_in = static_cast<unsigned>(rows.size());
        ++stats.matrices_built;
        stats.max_matrix_rows_seen = std::max(
            stats.max_matrix_rows_seen, static_cast<unsigned>(rows.size()));
        if (rows.empty()) {
            return res;
        }

        if (options.max_matrix_rows > 0 && rows.size() > options.max_matrix_rows) {
            throw GroebnerLimitExceeded{};
        }

        std::unordered_set<uint64_t> unique_cols;
        for (const auto& r : rows) {
            for (auto col : r.columns) {
                unique_cols.insert(col);
            }
        }

        if (options.max_matrix_columns > 0 && unique_cols.size() > options.max_matrix_columns) {
            throw GroebnerLimitExceeded{};
        }

        std::vector<uint64_t> cols(unique_cols.begin(), unique_cols.end());
        std::sort(cols.begin(), cols.end(), std::greater<uint64_t>());
        stats.max_matrix_columns_seen = std::max(
            stats.max_matrix_columns_seen, static_cast<unsigned>(cols.size()));

        size_t num_cols = cols.size();
        std::unordered_map<uint64_t, size_t> col_to_dense;
        col_to_dense.reserve(num_cols);
        for (size_t i = 0; i < num_cols; ++i) {
            col_to_dense[cols[i]] = i;
        }

        size_t num_rows = rows.size();
        std::vector<std::vector<Coeff>> mat(num_rows, std::vector<Coeff>(num_cols, dom.zero()));
        for (size_t i = 0; i < num_rows; ++i) {
            for (size_t j = 0; j < rows[i].columns.size(); ++j) {
                uint64_t col_id = rows[i].columns[j];
                mat[i][col_to_dense[col_id]] = rows[i].values[j];
            }
        }

        size_t r = 0;
        unsigned pivot_steps = 0;
        for (size_t c = 0; c < num_cols; ++c) {
            if (r >= num_rows) break;

            size_t pivot = r;
            while (pivot < num_rows && dom.is_zero(mat[pivot][c])) {
                pivot++;
            }

            if (pivot == num_rows) {
                continue;
            }

            if (pivot != r) {
                std::swap(mat[r], mat[pivot]);
                std::swap(rows[r].signature, rows[pivot].signature);
            }

            pivot_steps++;
            if (should_check_cancellation(pivot_steps, options)) {
                check_cancellation(options);
            }

            Coeff p_val = mat[r][c];
            if (!dom.is_one(p_val)) {
                Coeff inv_p = dom.inv(p_val);
                for (size_t j = c; j < num_cols; ++j) {
                    if (!dom.is_zero(mat[r][j])) {
                        mat[r][j] = dom.mul(mat[r][j], inv_p);
                    }
                }
            }

            for (size_t i = 0; i < num_rows; ++i) {
                if (i != r && !dom.is_zero(mat[i][c])) {
                    Coeff factor = mat[i][c];
                    for (size_t j = c; j < num_cols; ++j) {
                        if (!dom.is_zero(mat[r][j])) {
                            Coeff prod = dom.mul(factor, mat[r][j]);
                            mat[i][j] = dom.sub(mat[i][j], prod);
                        }
                    }
                }
            }
            r++;
        }

        std::vector<MatrixRow<Coeff>> reduced_rows;
        reduced_rows.reserve(num_rows);
        for (size_t i = 0; i < num_rows; ++i) {
            MatrixRow<Coeff> new_row;
            new_row.signature = rows[i].signature;
            new_row.columns.clear();
            new_row.values.clear();
            for (size_t j = 0; j < num_cols; ++j) {
                if (!dom.is_zero(mat[i][j])) {
                    new_row.columns.push_back(cols[j]);
                    new_row.values.push_back(mat[i][j]);
                }
            }
            if (new_row.columns.empty()) {
                res.rows_zero++;
                stats.reductions_to_zero++;
                stats.rows_reduced_to_zero++;
            } else {
                res.rows_out++;
                reduced_rows.push_back(std::move(new_row));
            }
        }
        rows = std::move(reduced_rows);
        return res;
    }
};

struct CodecWorkspace {
    std::vector<exponent_t> exps_a;
    std::vector<exponent_t> exps_b;
    std::vector<exponent_t> res;
};

inline CodecWorkspace& get_codec_workspace() {
    thread_local CodecWorkspace ws;
    return ws;
}

using OrdinalMonomialId = uint64_t;

class RuntimeDegRevLexCodec {
    unsigned num_vars_;
    unsigned max_degree_;
    std::vector<OrdinalMonomialId> L_;
    std::vector<std::vector<std::vector<OrdinalMonomialId>>> T2_;
    std::vector<std::vector<std::vector<OrdinalMonomialId>>> T3_;
    OrdinalMonomialId max_value_;

    uint64_t binomial_coefficient(size_t n, size_t k) const {
        if (k > n) return 0;
        if (k == 0 || k == n) return 1;
        if (k > n - k) k = n - k;
        
        std::vector<std::vector<uint64_t>> table(n + 1);
        for (size_t i = 0; i <= n; ++i) {
            table[i].resize(i + 1, 0);
            table[i][0] = 1;
            table[i][i] = 1;
            for (size_t j = 1; j < i; ++j) {
                uint64_t val = table[i - 1][j - 1] + table[i - 1][j];
                if (val < table[i - 1][j - 1]) {
                    throw GroebnerUnsupportedDomain{};
                }
                table[i][j] = val;
            }
        }
        return table[n][k];
    }

    uint64_t multiset_coefficient(size_t n, size_t k) const {
        if (n == 0) {
            return (k == 0) ? 1 : 0;
        }
        return binomial_coefficient(n + k - 1, k);
    }

    unsigned decode_i(unsigned d, unsigned& l, OrdinalMonomialId& v) const {
        unsigned i = l;
        while (i > 0 && v >= T3_[d][l][i - 1]) {
            --i;
        }
        v -= T3_[d][l][i];
        l = i - 1;
        return i;
    }

    unsigned decode_e(unsigned& d, unsigned i, OrdinalMonomialId& v) const {
        unsigned e = 1;
        while (v < T2_[d][i][e]) {
            ++e;
        }
        v -= T2_[d][i][e];
        d -= e;
        return e;
    }

public:
    RuntimeDegRevLexCodec(unsigned num_vars, unsigned max_degree)
        : num_vars_(num_vars), max_degree_(max_degree) {
        
        unsigned N = num_vars_;
        unsigned D = max_degree_;

        std::vector<std::vector<uint64_t>> multisetcoef(N + 1, std::vector<uint64_t>(D + 1, 0));
        for (unsigned i = 0; i <= N; ++i) {
            for (unsigned d = 0; d <= D; ++d) {
                multisetcoef[i][d] = multiset_coefficient(i, d);
            }
        }

        L_.resize(D + 2, 0);
        L_[0] = 0;
        L_[1] = 1;
        for (unsigned i = 2; i <= D + 1; ++i) {
            L_[i] = L_[i - 1] + multisetcoef[N][i - 1];
            if (L_[i] < L_[i - 1]) {
                throw GroebnerUnsupportedDomain{};
            }
        }
        max_value_ = L_[D + 1] - 1;

        T2_.resize(D + 1, std::vector<std::vector<OrdinalMonomialId>>(N, std::vector<OrdinalMonomialId>(D + 1, 0)));
        T3_.resize(D + 1, std::vector<std::vector<OrdinalMonomialId>>(N, std::vector<OrdinalMonomialId>(N, 0)));

        for (unsigned d = 0; d <= D; ++d) {
            for (unsigned i = 0; i < N; ++i) {
                for (unsigned e = 0; e <= d; ++e) {
                    uint64_t c = 0;
                    for (unsigned j = e + 1; j <= d; ++j) {
                        c += multisetcoef[i][d - j];
                    }
                    T2_[d][i][e] = c;
                }

                for (unsigned l = i; l < N; ++l) {
                    uint64_t c = 0;
                    for (unsigned j = 1; j <= d; ++j) {
                        if (i < l) {
                            uint64_t term = multisetcoef[l - i][j] * multisetcoef[i + 1][d - j];
                            c += term;
                        }
                    }
                    T3_[d][l][i] = c;
                }
            }
        }
    }

    unsigned degree(OrdinalMonomialId index) const {
        unsigned d = 0;
        while (d <= max_degree_ && L_[d + 1] <= index) {
            ++d;
        }
        return d;
    }

    OrdinalMonomialId max_id() const {
        return max_value_;
    }

    OrdinalMonomialId next(OrdinalMonomialId id) const { return id + 1; }

    OrdinalMonomialId to_index(const std::vector<exponent_t>& exps) const {
        unsigned d = 0;
        for (auto e : exps) {
            d += e;
        }
        if (d > max_degree_) {
            throw GroebnerLimitExceeded{};
        }

        OrdinalMonomialId v = L_[d];
        unsigned cur_d = d;
        unsigned l = num_vars_;

        for (size_t i = num_vars_; i > 0; --i) {
            size_t idx = i - 1;
            if (exps[idx] > 0) {
                v += T3_[cur_d][l - 1][idx] + T2_[cur_d][idx][exps[idx]];
                cur_d -= exps[idx];
                l = static_cast<unsigned>(idx);
            }
        }
        return v;
    }

    void from_index(OrdinalMonomialId id, std::vector<exponent_t>& out) const {
        out.assign(num_vars_, 0);
        unsigned d = degree(id);
        OrdinalMonomialId v = id - L_[d];
        unsigned li = num_vars_ - 1;
        unsigned cur_d = d;
        while (cur_d > 0) {
            unsigned idx = decode_i(cur_d, li, v);
            unsigned exp = decode_e(cur_d, idx, v);
            out[idx] = exp;
        }
    }

    OrdinalMonomialId multiply(OrdinalMonomialId a, OrdinalMonomialId b) const {
        auto& ws = get_codec_workspace();
        from_index(a, ws.exps_a);
        from_index(b, ws.exps_b);
        ws.res.resize(num_vars_);
        for (size_t i = 0; i < num_vars_; ++i) {
            ws.res[i] = ws.exps_a[i] + ws.exps_b[i];
        }
        return to_index(ws.res);
    }

    bool divides(OrdinalMonomialId a, OrdinalMonomialId b) const {
        auto& ws = get_codec_workspace();
        from_index(a, ws.exps_a);
        from_index(b, ws.exps_b);
        for (size_t i = 0; i < num_vars_; ++i) {
            if (ws.exps_a[i] > ws.exps_b[i]) return false;
        }
        return true;
    }

    OrdinalMonomialId divide(OrdinalMonomialId a, OrdinalMonomialId b) const {
        auto& ws = get_codec_workspace();
        from_index(a, ws.exps_a);
        from_index(b, ws.exps_b);
        ws.res.resize(num_vars_);
        for (size_t i = 0; i < num_vars_; ++i) {
            ws.res[i] = ws.exps_a[i] - ws.exps_b[i];
        }
        return to_index(ws.res);
    }

    OrdinalMonomialId lcm(OrdinalMonomialId a, OrdinalMonomialId b) const {
        auto& ws = get_codec_workspace();
        from_index(a, ws.exps_a);
        from_index(b, ws.exps_b);
        ws.res.resize(num_vars_);
        for (size_t i = 0; i < num_vars_; ++i) {
            ws.res[i] = std::max(ws.exps_a[i], ws.exps_b[i]);
        }
        return to_index(ws.res);
    }
};

template <typename Coeff, typename Domain>
GPoly<Coeff, Domain> normal_form_impl(GPoly<Coeff, Domain> p, const std::vector<GPoly<Coeff, Domain>> &G, const GroebnerOptions &options, unsigned &reduction_steps) {
    GPoly<Coeff, Domain> remainder(p.order);
    remainder.dom = p.dom;
    Domain dom = p.dom;
    while (!p.is_zero()) {
        check_cancellation(options);
        if (options.max_reduction_steps > 0 && reduction_steps >= options.max_reduction_steps) {
            throw GroebnerLimitExceeded{};
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
    auto polys_equal = [](const GPoly<Coeff, Domain> &a, const GPoly<Coeff, Domain> &b) {
        if (a.terms.size() != b.terms.size()) return false;
        Domain dom = a.dom;
        for (size_t k = 0; k < a.terms.size(); ++k) {
            if (a.terms[k].monomial != b.terms[k].monomial) return false;
            if (!dom.is_zero(dom.sub(a.terms[k].coeff, b.terms[k].coeff))) return false;
        }
        return true;
    };

    bool changed = true;
    while (changed) {
        changed = false;
        for (size_t i = 0; i < G.size(); ++i) {
            GPoly<Coeff, Domain> g = G[i];
            G.erase(G.begin() + i);
            GPoly<Coeff, Domain> reduced_g = normal_form_impl(g, G, options, reduction_steps);
            if (!reduced_g.is_zero()) {
                reduced_g.make_monic();
                G.insert(G.begin() + i, reduced_g);
                if (!polys_equal(reduced_g, g)) {
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

} // namespace SymEngine

#endif // SYMENGINE_POLYS_GROEBNER_INTERNAL_H
