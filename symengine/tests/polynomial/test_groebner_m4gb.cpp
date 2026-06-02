#include "catch.hpp"
#include <symengine/polys/groebner.h>
#include <symengine/polys/groebner_internal.h>
#include <symengine/symbol.h>
#include <symengine/integer.h>
#include <symengine/rational.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/symengine_exception.h>
#include <string>
#include <vector>

using namespace SymEngine;

namespace {

bool basis_is_canonical(const GroebnerResult &res)
{
    return res.status == GroebnerStatus::Success
           && is_groebner(res.basis, res.variables, res.order)
           && is_reduced_basis(res.basis, res.variables, res.order);
}

RCP<const Integer> Z(long n)
{
    return integer(n);
}

RCP<const Basic> P(const RCP<const Basic> &x, unsigned n)
{
    return pow(x, integer(n));
}

std::vector<std::string> basis_strings(const GroebnerResult &res)
{
    std::vector<std::string> strings;
    for (const auto &poly : res.basis) {
        strings.push_back(poly->__str__());
    }
    return strings;
}

} // namespace

TEST_CASE("M4GB basic tests", "[groebner_m4gb]")
{
    auto x = symbol("x");
    auto y = symbol("y");

    // {x^3 - 2*x*y, x^2*y + x - 2*y^2} mod 5
    vec_basic polys = {
        sub(pow(x, integer(3)), mul(integer(2), mul(x, y))),
        sub(add(mul(pow(x, integer(2)), y), x), mul(integer(2), pow(y, integer(2))))
    };
    vec_sym vars = {x, y};
    GroebnerOptions opt;
    opt.order = MonomialOrder::DegRevLex;
    opt.algorithm = GroebnerAlgorithm::M4GB;
    opt.modulus = 5;

    auto res = groebner_basis(polys, vars, opt);
    REQUIRE(basis_is_canonical(res));
    REQUIRE(res.basis.size() == 3);
    REQUIRE(res.basis[0]->__str__() == "x**2");
    REQUIRE(res.basis[1]->__str__() == "x*y");
    REQUIRE(res.basis[2]->__str__() == "2*x + y**2");

    GroebnerOptions buchberger = opt;
    buchberger.algorithm = GroebnerAlgorithm::Buchberger;
    auto ref = groebner_basis(polys, vars, buchberger);
    REQUIRE(ref.status == GroebnerStatus::Success);
    REQUIRE(basis_strings(res) == basis_strings(ref));
}

TEST_CASE("M4GB Katsura-3 test", "[groebner_m4gb]")
{
    auto x0 = symbol("x0");
    auto x1 = symbol("x1");
    auto x2 = symbol("x2");
    vec_basic polys = {
        sub(add(add(x0, mul(Z(2), x1)), mul(Z(2), x2)), Z(1)),
        sub(add(add(P(x0, 2), mul(Z(2), P(x1, 2))), mul(Z(2), P(x2, 2))), x0),
        sub(add(mul(Z(2), mul(x0, x1)), mul(Z(2), mul(x1, x2))), x1)
    };
    vec_sym vars = {x0, x1, x2};
    GroebnerOptions opt;
    opt.order = MonomialOrder::DegRevLex;
    opt.algorithm = GroebnerAlgorithm::M4GB;
    opt.modulus = 31;

    auto res = groebner_basis(polys, vars, opt);
    REQUIRE(res.status == GroebnerStatus::Success);
    REQUIRE(res.basis.size() == 4);

    // Cross-verify with F5B under same options (modulo 31)
    GroebnerOptions opt_f5b = opt;
    opt_f5b.algorithm = GroebnerAlgorithm::F5B;
    auto res_f5b = groebner_basis(polys, vars, opt_f5b);
    REQUIRE(res_f5b.status == GroebnerStatus::Success);
    REQUIRE(res_f5b.basis.size() == 4);

    for (size_t i = 0; i < 4; ++i) {
        REQUIRE(res.basis[i]->__str__() == res_f5b.basis[i]->__str__());
    }

    GroebnerOptions opt_buchberger = opt;
    opt_buchberger.algorithm = GroebnerAlgorithm::Buchberger;
    auto res_buchberger = groebner_basis(polys, vars, opt_buchberger);
    REQUIRE(res_buchberger.status == GroebnerStatus::Success);
    REQUIRE(basis_strings(res) == basis_strings(res_buchberger));
}

TEST_CASE("M4GB preconditions and cancellation interval edge case", "[groebner_m4gb]")
{
    auto x = symbol("x");
    auto y = symbol("y");

    vec_basic polys = {
        sub(pow(x, integer(3)), mul(integer(2), mul(x, y))),
        sub(add(mul(pow(x, integer(2)), y), x), mul(integer(2), pow(y, integer(2))))
    };

    GroebnerOptions no_modulus;
    no_modulus.order = MonomialOrder::DegRevLex;
    no_modulus.algorithm = GroebnerAlgorithm::M4GB;
    auto no_modulus_res = groebner_basis(polys, {x, y}, no_modulus);
    REQUIRE(no_modulus_res.status == GroebnerStatus::UnsupportedCoefficientDomain);

    GroebnerOptions wrong_order = no_modulus;
    wrong_order.order = MonomialOrder::Lex;
    wrong_order.modulus = 5;
    auto wrong_order_res = groebner_basis(polys, {x, y}, wrong_order);
    REQUIRE(wrong_order_res.status == GroebnerStatus::UnsupportedCoefficientDomain);

    GroebnerOptions zero_interval = no_modulus;
    zero_interval.modulus = 5;
    zero_interval.cancellation_check_interval = 0;
    auto zero_interval_res = groebner_basis(polys, {x, y}, zero_interval);
    REQUIRE(zero_interval_res.status == GroebnerStatus::Success);
}

TEST_CASE("M4GB Codec Round-trips and Monomial Operations", "[groebner_m4gb]")
{
    unsigned num_vars = 3;
    unsigned max_degree = 6;
    RuntimeDegRevLexCodec codec(num_vars, max_degree);

    std::vector<std::vector<unsigned>> test_monomials = {
        {0, 0, 0},
        {1, 0, 0},
        {0, 2, 0},
        {0, 0, 3},
        {1, 1, 1},
        {2, 0, 4},
        {0, 5, 1}
    };

    for (const auto& mon : test_monomials) {
        auto idx = codec.to_index(mon);
        std::vector<unsigned> round_trip;
        codec.from_index(idx, round_trip);
        REQUIRE(mon == round_trip);
    }

    std::vector<unsigned> m1 = {1, 1, 0};
    std::vector<unsigned> m2 = {0, 1, 2};
    auto id1 = codec.to_index(m1);
    auto id2 = codec.to_index(m2);

    auto id_mul = codec.multiply(id1, id2);
    std::vector<unsigned> res_mul;
    codec.from_index(id_mul, res_mul);
    REQUIRE(res_mul == std::vector<unsigned>{1, 2, 2});

    REQUIRE(codec.divides(id1, id_mul));
    REQUIRE(codec.divides(id2, id_mul));
    REQUIRE(!codec.divides(id_mul, id1));

    auto id_div = codec.divide(id_mul, id1);
    std::vector<unsigned> res_div;
    codec.from_index(id_div, res_div);
    REQUIRE(res_div == m2);

    auto id_lcm = codec.lcm(id1, id2);
    std::vector<unsigned> res_lcm;
    codec.from_index(id_lcm, res_lcm);
    REQUIRE(res_lcm == std::vector<unsigned>{1, 1, 2});
}
