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

TEST_CASE("MoGVW basic tests", "[groebner_mogvw]")
{
    auto x = symbol("x");
    auto y = symbol("y");

    // {x^3 - 2*x*y, x^2*y + x - 2*y^2} grlex -> {x^2, x*y, -1/2*x + y^2}
    vec_basic polys3 = {sub(pow(x, integer(3)), mul(integer(2), mul(x, y))),
                        sub(add(mul(pow(x, integer(2)), y), x), mul(integer(2), pow(y, integer(2))))};
    vec_sym vars3 = {x, y};
    GroebnerOptions opt3;
    opt3.order = MonomialOrder::GrLex;
    opt3.algorithm = GroebnerAlgorithm::MoGVW;

    auto res3 = groebner_basis(polys3, vars3, opt3);
    REQUIRE(basis_is_canonical(res3));
    REQUIRE(res3.basis.size() == 3);
    REQUIRE(res3.basis[0]->__str__() == "x**2");
    REQUIRE(res3.basis[1]->__str__() == "x*y");
    REQUIRE(res3.basis[2]->__str__() == "(-1/2)*x + y**2");
}

TEST_CASE("MoGVW Katsura-3 grlex test", "[groebner_mogvw]")
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
    opt.order = MonomialOrder::GrLex;
    opt.algorithm = GroebnerAlgorithm::MoGVW;

    auto res = groebner_basis(polys, vars, opt);
    REQUIRE(basis_is_canonical(res));
    REQUIRE(res.basis.size() == 4);
}

TEST_CASE("MoGVW modular GF(p) test", "[groebner_mogvw]")
{
    auto x = symbol("x");
    auto y = symbol("y");

    vec_basic polys = {
        sub(pow(x, integer(3)), mul(integer(2), mul(x, y))),
        sub(add(mul(pow(x, integer(2)), y), x), mul(integer(2), pow(y, integer(2))))
    };
    vec_sym vars = {x, y};
    GroebnerOptions opt;
    opt.order = MonomialOrder::GrLex;
    opt.algorithm = GroebnerAlgorithm::MoGVW;
    opt.modulus = 5;

    auto res = groebner_basis(polys, vars, opt);
    REQUIRE(res.status == GroebnerStatus::Success);
    REQUIRE(res.basis.size() == 3);
    REQUIRE(res.basis[0]->__str__() == "x**2");
    REQUIRE(res.basis[1]->__str__() == "x*y");
    REQUIRE(res.basis[2]->__str__() == "2*x + y**2"); // -1/2 modulo 5 is 2!
}
