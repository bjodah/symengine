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

} // namespace

TEST_CASE("F5B basic QQ tests", "[groebner_f5b]")
{
    auto x = symbol("x");
    auto y = symbol("y");

    // {x^2 + 2*x*y^2, x*y + 2*y^3 - 1} -> {x, y^3 - 1/2}
    vec_basic polys1 = {add(pow(x, integer(2)), mul(integer(2), mul(x, pow(y, integer(2))))),
                        sub(add(mul(x, y), mul(integer(2), pow(y, integer(3)))), integer(1))};
    vec_sym vars1 = {x, y};
    GroebnerOptions opt1;
    opt1.order = MonomialOrder::Lex;
    opt1.algorithm = GroebnerAlgorithm::F5B;
    
    auto res1 = groebner_basis(polys1, vars1, opt1);
    REQUIRE(basis_is_canonical(res1));
    REQUIRE(res1.basis.size() == 2);
    REQUIRE(res1.basis[0]->__str__() == "x");
    REQUIRE(res1.basis[1]->__str__() == "-1/2 + y**3");

    // Katsura-3 lex
    auto x0 = symbol("x0");
    auto x1 = symbol("x1");
    auto x2 = symbol("x2");
    vec_basic polys2 = {
        sub(add(add(x0, mul(Z(2), x1)), mul(Z(2), x2)), Z(1)),
        sub(add(add(P(x0, 2), mul(Z(2), P(x1, 2))), mul(Z(2), P(x2, 2))), x0),
        sub(add(mul(Z(2), mul(x0, x1)), mul(Z(2), mul(x1, x2))), x1)
    };
    vec_sym vars2 = {x0, x1, x2};
    GroebnerOptions opt2;
    opt2.order = MonomialOrder::Lex;
    opt2.algorithm = GroebnerAlgorithm::F5B;
    auto res2 = groebner_basis(polys2, vars2, opt2);
    REQUIRE(basis_is_canonical(res2));
    REQUIRE(res2.stats.rejected_by_syzygy + res2.stats.rejected_by_rewritten > 0);
}

TEST_CASE("F5B benchmark minpoly test", "[groebner_f5b]")
{
    auto x = symbol("x");
    auto y = symbol("y");
    auto z = symbol("z");

    vec_basic polys = {
        add(add(P(x, 3), x), Z(1)),
        add(add(P(y, 2), y), Z(1)),
        sub(mul(add(x, y), z), add(P(x, 2), y))
    };
    vec_sym vars = {x, y, z};
    GroebnerOptions opt;
    opt.order = MonomialOrder::Lex;
    opt.algorithm = GroebnerAlgorithm::F5B;
    
    auto res = groebner_basis(polys, vars, opt);
    REQUIRE(basis_is_canonical(res));
    REQUIRE(res.basis.size() == 3);
}

TEST_CASE("F5B modular prime test", "[groebner_f5b]")
{
    auto x = symbol("x");
    auto y = symbol("y");

    vec_basic polys = {
        add(pow(x, integer(2)), mul(integer(2), mul(x, pow(y, integer(2))))),
        sub(add(mul(x, y), mul(integer(2), pow(y, integer(3)))), integer(1))
    };
    vec_sym vars = {x, y};
    GroebnerOptions opt;
    opt.order = MonomialOrder::Lex;
    opt.algorithm = GroebnerAlgorithm::F5B;
    opt.modulus = 5; // GF(5)
    
    auto res = groebner_basis(polys, vars, opt);
    REQUIRE(res.status == GroebnerStatus::Success);
    REQUIRE(res.basis.size() == 2);
    // modulo 5, -1/2 is -3, which is 2. So y^3 - 1/2 is y^3 + 2!
    REQUIRE(res.basis[0]->__str__() == "x");
    REQUIRE(res.basis[1]->__str__() == "2 + y**3");
}

TEST_CASE("F5B honors public limits and modular domain preconditions", "[groebner_f5b]")
{
    auto w = symbol("w");
    auto x = symbol("x");
    auto y = symbol("y");
    auto z = symbol("z");
    auto a = symbol("a");

    vec_basic toric = {
        sub(pow(w, integer(2)), mul(x, z)),
        sub(mul(w, x), mul(y, z)),
        sub(pow(x, integer(2)), mul(w, y)),
        sub(mul(x, y), pow(z, integer(2))),
        sub(pow(y, integer(2)), mul(w, z)),
    };

    GroebnerOptions limited;
    limited.order = MonomialOrder::DegRevLex;
    limited.algorithm = GroebnerAlgorithm::F5B;
    limited.max_s_pairs = 1;
    auto limited_res = groebner_basis(toric, {w, x, y, z}, limited);
    REQUIRE(limited_res.status == GroebnerStatus::ResourceLimitExceeded);

    GroebnerOptions composite;
    composite.algorithm = GroebnerAlgorithm::F5B;
    composite.modulus = 6;
    auto composite_res = groebner_basis({add(pow(x, integer(2)), y)}, {x, y}, composite);
    REQUIRE(composite_res.status == GroebnerStatus::UnsupportedCoefficientDomain);

    GroebnerOptions symbolic_mod;
    symbolic_mod.algorithm = GroebnerAlgorithm::F5B;
    symbolic_mod.modulus = 5;
    auto symbolic_res = groebner_basis({add(mul(a, x), integer(1))}, {x}, symbolic_mod);
    REQUIRE(symbolic_res.status == GroebnerStatus::UnsupportedCoefficientDomain);

    GroebnerOptions bad_denominator;
    bad_denominator.algorithm = GroebnerAlgorithm::F5B;
    bad_denominator.modulus = 5;
    auto rational_coeff = Rational::from_two_ints(1, 5);
    auto denominator_res = groebner_basis({add(mul(rational_coeff, x), integer(1))},
                                          {x},
                                          bad_denominator);
    REQUIRE(denominator_res.status == GroebnerStatus::UnsupportedCoefficientDomain);
}
