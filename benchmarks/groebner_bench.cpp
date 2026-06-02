#include <symengine/polys/groebner.h>
#include <symengine/add.h>
#include <symengine/integer.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/rational.h>
#include <symengine/symbol.h>

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace SymEngine;

namespace {

RCP<const Integer> Z(long n)
{
    return integer(n);
}

RCP<const Basic> P(const RCP<const Basic> &x, unsigned n)
{
    return pow(x, integer(n));
}

std::string status_name(GroebnerStatus status)
{
    switch (status) {
        case GroebnerStatus::Success:
            return "success";
        case GroebnerStatus::Cancelled:
            return "cancelled";
        case GroebnerStatus::ResourceLimitExceeded:
            return "limit";
        case GroebnerStatus::NotZeroDimensional:
            return "not_zero_dimensional";
        case GroebnerStatus::UnsupportedCoefficientDomain:
            return "unsupported_domain";
    }
    return "unknown";
}

std::string order_name(MonomialOrder order)
{
    switch (order) {
        case MonomialOrder::Lex:
            return "lex";
        case MonomialOrder::GrLex:
            return "grlex";
        case MonomialOrder::DegRevLex:
            return "degrevlex";
    }
    return "unknown";
}

struct BenchResult {
    GroebnerResult result;
    double seconds;
};

template <typename F>
BenchResult run_repeated(unsigned iterations, F f)
{
    GroebnerResult last;
    auto start = std::chrono::high_resolution_clock::now();
    for (unsigned i = 0; i < iterations; ++i) {
        last = f();
    }
    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = stop - start;
    return BenchResult{last, elapsed.count()};
}

std::string algorithm_name(GroebnerAlgorithm alg)
{
    switch (alg) {
        case GroebnerAlgorithm::Buchberger:
            return "buchberger";
        case GroebnerAlgorithm::F5B:
            return "f5b";
        case GroebnerAlgorithm::MoGVW:
            return "mogvw";
        case GroebnerAlgorithm::M4GB:
            return "m4gb";
        case GroebnerAlgorithm::Auto:
            return "auto";
    }
    return "unknown";
}

void print_result(const std::string &name, unsigned iterations, const BenchResult &res)
{
    const GroebnerStats &s = res.result.stats;
    std::cout << std::left << std::setw(34) << name
              << " total_s=" << std::right << std::setw(10) << std::fixed << std::setprecision(6) << res.seconds
              << " avg_ms=" << std::setw(10) << std::fixed << std::setprecision(3) << (1000.0 * res.seconds / iterations)
              << " status=" << std::setw(20) << std::left << status_name(res.result.status)
              << " alg=" << std::setw(12) << algorithm_name(res.result.selected_algorithm)
              << " order=" << std::setw(10) << order_name(res.result.order)
              << " basis=" << std::setw(4) << res.result.basis.size()
              << " spairs=" << std::setw(6) << s.s_pairs_processed
              << " zero=" << std::setw(6) << s.reductions_to_zero
              << " max_basis=" << std::setw(4) << s.max_basis_size
              << " syz=" << std::setw(6) << s.rejected_by_syzygy
              << " rew=" << std::setw(6) << s.rejected_by_rewritten
              << " f5red=" << std::setw(6) << s.f5_reductions
              << " lift=" << std::setw(6) << s.labeled_monomials_lifted
              << " lcmrej=" << std::setw(6) << s.rejected_by_lcm
              << " mat=" << std::setw(4) << s.matrices_built
              << " max_rows=" << std::setw(5) << s.max_matrix_rows_seen
              << " max_cols=" << std::setw(5) << s.max_matrix_columns_seen
              << " pm_entries=" << s.polymatrix_entries_built
              << std::endl;
}

GroebnerResult compute_basis(const vec_basic &polys,
                             const vec_sym &vars,
                             MonomialOrder order,
                             GroebnerAlgorithm algorithm = GroebnerAlgorithm::Buchberger,
                             uint64_t modulus = 0,
                             bool verify = false)
{
    GroebnerOptions opt;
    opt.order = order;
    opt.algorithm = algorithm;
    opt.modulus = modulus;
    opt.verify_with_buchberger = verify;
    return groebner_basis(polys, vars, opt);
}

GroebnerResult compute_degrevlex_then_lex(const vec_basic &polys, const vec_sym &vars)
{
    GroebnerOptions opt;
    opt.order = MonomialOrder::DegRevLex;
    auto drl = groebner_basis(polys, vars, opt);
    if (drl.status != GroebnerStatus::Success) {
        return drl;
    }
    auto lex = fglm_convert(drl, MonomialOrder::Lex, opt);
    // fglm_convert currently reports conversion status but not source basis
    // stats. Carry those through here so direct-lex and drl+FGLM benchmark
    // rows remain comparable.
    lex.stats = drl.stats;
    return lex;
}

} // namespace

int main(int argc, char *argv[])
{
    unsigned iterations = 1;
    if (argc > 1) {
        iterations = std::max(1, std::atoi(argv[1]));
    }

    auto x = symbol("x");
    auto y = symbol("y");
    auto z = symbol("z");
    auto w = symbol("w");
    auto a = symbol("a");
    auto b = symbol("b");
    auto c = symbol("c");
    auto d = symbol("d");
    auto x0 = symbol("x0");
    auto x1 = symbol("x1");
    auto x2 = symbol("x2");
    auto r1 = symbol("r1");
    auto r2 = symbol("r2");
    auto xs = symbol("xs");

    std::vector<std::pair<std::string, std::function<GroebnerResult()>>> cases;

    // SymPy test_groebnertools: small elimination, lex is much slower than a
    // favorable graded order in many implementations.
    {
        vec_basic polys = {
            add(P(x, 2), mul(Z(2), mul(x, P(y, 2)))),
            sub(add(mul(x, y), mul(Z(2), P(y, 3))), Z(1)),
        };
        vec_sym vars = {x, y};
        cases.push_back({"sympy_small_elim_lex", [=]() { return compute_basis(polys, vars, MonomialOrder::Lex); }});
        cases.push_back({"sympy_small_elim_lex_f5b", [=]() { return compute_basis(polys, vars, MonomialOrder::Lex, GroebnerAlgorithm::F5B); }});
        cases.push_back({"sympy_small_elim_degrevlex", [=]() { return compute_basis(polys, vars, MonomialOrder::DegRevLex); }});
        cases.push_back({"sympy_small_elim_drl_f5b", [=]() { return compute_basis(polys, vars, MonomialOrder::DegRevLex, GroebnerAlgorithm::F5B); }});
    }

    // SymPy minpoly benchmark from test_groebnertools.
    {
        vec_basic polys = {
            add(add(P(x, 3), x), Z(1)),
            add(add(P(y, 2), y), Z(1)),
            sub(mul(add(x, y), z), add(P(x, 2), y)),
        };
        vec_sym vars = {x, y, z};
        cases.push_back({"sympy_minpoly_lex", [=]() { return compute_basis(polys, vars, MonomialOrder::Lex); }});
        cases.push_back({"sympy_minpoly_lex_f5b", [=]() { return compute_basis(polys, vars, MonomialOrder::Lex, GroebnerAlgorithm::F5B); }});
        cases.push_back({"sympy_minpoly_drl_fglm", [=]() { return compute_degrevlex_then_lex(polys, vars); }});
    }

    // SymPy FGLM test: elementary symmetric equations with a degree-12 lex
    // eliminant in d.
    {
        vec_basic polys = {
            add(add(add(a, b), c), d),
            add(add(add(mul(a, b), mul(a, d)), mul(b, c)), mul(b, d)),
            add(add(add(mul(mul(a, b), c), mul(mul(a, b), d)), mul(mul(a, c), d)), mul(mul(b, c), d)),
            sub(mul(mul(mul(a, b), c), d), Z(1)),
        };
        vec_sym vars = {a, b, c, d};
        cases.push_back({"sympy_fglm_symmetric_lex", [=]() { return compute_basis(polys, vars, MonomialOrder::Lex); }});
        cases.push_back({"sympy_fglm_symmetric_drl_fglm", [=]() { return compute_degrevlex_then_lex(polys, vars); }});
    }

    // Katsura-3 from SymPy: useful small structured benchmark with moderately
    // growing coefficients under lex order.
    {
        vec_basic polys = {
            sub(add(add(x0, mul(Z(2), x1)), mul(Z(2), x2)), Z(1)),
            sub(add(add(P(x0, 2), mul(Z(2), P(x1, 2))), mul(Z(2), P(x2, 2))), x0),
            sub(add(mul(Z(2), mul(x0, x1)), mul(Z(2), mul(x1, x2))), x1),
        };
        vec_sym vars = {x0, x1, x2};
        cases.push_back({"katsura3_lex", [=]() { return compute_basis(polys, vars, MonomialOrder::Lex); }});
        cases.push_back({"katsura3_lex_f5b", [=]() { return compute_basis(polys, vars, MonomialOrder::Lex, GroebnerAlgorithm::F5B); }});
        cases.push_back({"katsura3_drl", [=]() { return compute_basis(polys, vars, MonomialOrder::DegRevLex); }});
        cases.push_back({"katsura3_drl_f5b", [=]() { return compute_basis(polys, vars, MonomialOrder::DegRevLex, GroebnerAlgorithm::F5B); }});
        cases.push_back({"katsura3_grlex_mogvw", [=]() { return compute_basis(polys, vars, MonomialOrder::GrLex, GroebnerAlgorithm::MoGVW); }});
        cases.push_back({"katsura3_gf31_m4gb", [=]() { return compute_basis(polys, vars, MonomialOrder::DegRevLex, GroebnerAlgorithm::M4GB, 31); }});
        cases.push_back({"katsura3_gf31_auto", [=]() { return compute_basis(polys, vars, MonomialOrder::DegRevLex, GroebnerAlgorithm::Auto, 31); }});
        cases.push_back({"katsura3_gf31_m4gb_verify", [=]() { return compute_basis(polys, vars, MonomialOrder::DegRevLex, GroebnerAlgorithm::M4GB, 31, true); }});
    }

    // Singular-derived toric/binomial case from 01-IMPL-PLAN-GROEBNER.md.
    {
        vec_basic polys = {
            sub(P(w, 2), mul(x, z)),
            sub(mul(w, x), mul(y, z)),
            sub(P(x, 2), mul(w, y)),
            sub(mul(x, y), P(z, 2)),
            sub(P(y, 2), mul(w, z)),
        };
        vec_sym vars = {w, x, y, z};
        cases.push_back({"toric_binomial_drl", [=]() { return compute_basis(polys, vars, MonomialOrder::DegRevLex); }});
        cases.push_back({"toric_binomial_drl_f5b", [=]() { return compute_basis(polys, vars, MonomialOrder::DegRevLex, GroebnerAlgorithm::F5B); }});
    }

    // Inspired by demo_groebner_bench_sympy.py: eliminate two algebraic roots
    // of the same cubic from an expression in xs. This avoids CRootOf in C++
    // while preserving the "composition through root variables" shape.
    {
        auto cubic = [](const RCP<const Basic> &v) {
            return add(add(add(mul(Z(4000000), P(v, 3)),
                               mul(Z(-239960000), P(v, 2))),
                           mul(Z(4782399900), v)),
                       Z(-31663998001));
        };
        vec_basic polys = {
            cubic(r1),
            cubic(r2),
            sub(xs, add(add(mul(Z(17), r1), mul(Z(-29), r2)), mul(r1, r2))),
        };
        vec_sym vars = {r1, r2, xs};
        cases.push_back({"algebraic_composition_lex", [=]() { return compute_basis(polys, vars, MonomialOrder::Lex); }});
        cases.push_back({"algebraic_composition_drl_fglm", [=]() { return compute_degrevlex_then_lex(polys, vars); }});
    }

    std::cout << "Groebner benchmark iterations: " << iterations << std::endl;
    for (const auto &entry : cases) {
        print_result(entry.first, iterations, run_repeated(iterations, entry.second));
    }

    return 0;
}
