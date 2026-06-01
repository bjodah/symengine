#include "catch.hpp"
#include <symengine/polys/groebner.h>
#include <symengine/polys/groebner_internal.h>
#include <symengine/symbol.h>
#include <symengine/integer.h>
#include <symengine/rational.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/sets.h>
#include <symengine/tuple.h>
#include <symengine/symengine_exception.h>

using namespace SymEngine;

TEST_CASE("Monomial Ordering Tests", "[groebner]")
{
    // lex((1,2,3)) == (1,2,3)
    Monomial m1 = {1, 2, 3};
    Monomial m2 = {1, 2, 3};
    REQUIRE(!compare_monomial_less(m1, m2, MonomialOrder::Lex));
    REQUIRE(!compare_monomial_less(m2, m1, MonomialOrder::Lex));

    // grlex((1,2,3))
    Monomial m3 = {0, 0, 2}; // deg 2
    Monomial m4 = {0, 1, 1}; // deg 2
    // Equal-degree comparison checks:
    // (0, 1, 1) > (0, 0, 2) under grlex and degrevlex
    REQUIRE(compare_monomial_less(m3, m4, MonomialOrder::GrLex));
    REQUIRE(compare_monomial_less(m3, m4, MonomialOrder::DegRevLex));

    // (0, 3, 1) < (2, 2, 1) under grlex and degrevlex
    Monomial m5 = {0, 3, 1}; // deg 4
    Monomial m6 = {2, 2, 1}; // deg 4
    REQUIRE(compare_monomial_less(m5, m6, MonomialOrder::GrLex));
    REQUIRE(compare_monomial_less(m5, m6, MonomialOrder::DegRevLex));

    // With generator order [z, y, x], sort:
    // {x, x^2*z^2, x*y, x^2, 1, y^2, x^3, y, z, x*y^2*z, x^2*y^2}
    // Exponent vectors under [z, y, x]:
    // 1: (0, 0, 0)
    // x: (0, 0, 1)
    // x^2: (0, 0, 2)
    // x^3: (0, 0, 3)
    // y: (0, 1, 0)
    // x*y: (0, 1, 1)
    // y^2: (0, 2, 0)
    // x^2*y^2: (0, 2, 2)
    // z: (1, 0, 0)
    // x*y^2*z: (1, 2, 1)
    // x^2*z^2: (2, 0, 2)

    std::vector<Monomial> list = {
        {0, 0, 1}, // x
        {2, 0, 2}, // x^2*z^2
        {0, 1, 1}, // x*y
        {0, 0, 2}, // x^2
        {0, 0, 0}, // 1
        {0, 2, 0}, // y^2
        {0, 0, 3}, // x^3
        {0, 1, 0}, // y
        {1, 0, 0}, // z
        {1, 2, 1}, // x*y^2*z
        {0, 2, 2}  // x^2*y^2
    };

    // lex sorted (ascending): {1, x, x^2, x^3, y, x*y, y^2, x^2*y^2, z, x*y^2*z, x^2*z^2}
    std::vector<Monomial> lex_expected = {
        {0, 0, 0}, {0, 0, 1}, {0, 0, 2}, {0, 0, 3}, {0, 1, 0},
        {0, 1, 1}, {0, 2, 0}, {0, 2, 2}, {1, 0, 0}, {1, 2, 1}, {2, 0, 2}
    };
    std::vector<Monomial> lex_actual = list;
    std::sort(lex_actual.begin(), lex_actual.end(), [](const Monomial &a, const Monomial &b) {
        return compare_monomial_less(a, b, MonomialOrder::Lex);
    });
    REQUIRE(lex_actual == lex_expected);

    // grlex sorted (ascending): {1, x, y, z, x^2, x*y, y^2, x^3, x^2*y^2, x*y^2*z, x^2*z^2}
    std::vector<Monomial> grlex_expected = {
        {0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {1, 0, 0}, {0, 0, 2},
        {0, 1, 1}, {0, 2, 0}, {0, 0, 3}, {0, 2, 2}, {1, 2, 1}, {2, 0, 2}
    };
    std::vector<Monomial> grlex_actual = list;
    std::sort(grlex_actual.begin(), grlex_actual.end(), [](const Monomial &a, const Monomial &b) {
        return compare_monomial_less(a, b, MonomialOrder::GrLex);
    });
    REQUIRE(grlex_actual == grlex_expected);

    // degrevlex sorted (ascending): {1, x, y, z, x^2, x*y, y^2, x^3, x^2*y^2, x^2*z^2, x*y^2*z}
    std::vector<Monomial> degrevlex_expected = {
        {0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {1, 0, 0}, {0, 0, 2},
        {0, 1, 1}, {0, 2, 0}, {0, 0, 3}, {0, 2, 2}, {2, 0, 2}, {1, 2, 1}
    };
    std::vector<Monomial> degrevlex_actual = list;
    std::sort(degrevlex_actual.begin(), degrevlex_actual.end(), [](const Monomial &a, const Monomial &b) {
        return compare_monomial_less(a, b, MonomialOrder::DegRevLex);
    });
    REQUIRE(degrevlex_actual == degrevlex_expected);
}

TEST_CASE("Rational Buchberger Exact Groebner Basis Tests", "[groebner]")
{
    RCP<const Symbol> x = symbol("x");
    RCP<const Symbol> y = symbol("y");
    RCP<const Symbol> z = symbol("z");

    // Over QQ[x, y], lex:
    // input {x^2 + 2*x*y^2, x*y + 2*y^3 - 1};
    // expected {x, y^3 - 1/2}
    vec_basic polys1 = {add(pow(x, integer(2)), mul(integer(2), mul(x, pow(y, integer(2))))),
                        sub(add(mul(x, y), mul(integer(2), pow(y, integer(3)))), integer(1))};
    vec_sym vars1 = {x, y};
    GroebnerOptions opt1;
    opt1.order = MonomialOrder::Lex;
    auto res1 = groebner_basis(polys1, vars1, opt1);
    REQUIRE(res1.status == GroebnerStatus::Success);
    REQUIRE(res1.basis.size() == 2);
    REQUIRE(res1.basis[0]->__str__() == "x");
    REQUIRE(res1.basis[1]->__str__() == "-1/2 + y**3");

    // Over QQ[y, x], lex:
    // input {2*x^2*y + y^2, 2*x^3 + x*y - 1};
    // expected {y, x^3 - 1/2}
    vec_basic polys2 = {add(mul(integer(2), mul(pow(x, integer(2)), y)), pow(y, integer(2))),
                        sub(add(mul(integer(2), pow(x, integer(3))), mul(x, y)), integer(1))};
    vec_sym vars2 = {y, x};
    GroebnerOptions opt2;
    opt2.order = MonomialOrder::Lex;
    auto res2 = groebner_basis(polys2, vars2, opt2);
    REQUIRE(res2.status == GroebnerStatus::Success);
    REQUIRE(res2.basis.size() == 2);
    REQUIRE(res2.basis[0]->__str__() == "y");
    REQUIRE(res2.basis[1]->__str__() == "-1/2 + x**3");

    // Over QQ[x, y], grlex:
    // input {x^3 - 2*x*y, x^2*y + x - 2*y^2};
    // expected {x^2, x*y, -1/2*x + y^2}
    vec_basic polys3 = {sub(pow(x, integer(3)), mul(integer(2), mul(x, y))),
                        sub(add(mul(pow(x, integer(2)), y), x), mul(integer(2), pow(y, integer(2))))};
    vec_sym vars3 = {x, y};
    GroebnerOptions opt3;
    opt3.order = MonomialOrder::GrLex;
    auto res3 = groebner_basis(polys3, vars3, opt3);
    REQUIRE(res3.status == GroebnerStatus::Success);
    REQUIRE(res3.basis.size() == 3);
    REQUIRE(res3.basis[0]->__str__() == "x**2");
    REQUIRE(res3.basis[1]->__str__() == "x*y");
    REQUIRE(res3.basis[2]->__str__() == "(-1/2)*x + y**2");

    // Over QQ[x, y, z], lex:
    // input {-x^2 + y, -x^3 + z};
    // expected {x^2 - y, x*y - z, x*z - y^2, y^3 - z^2}
    vec_basic polys4 = {sub(y, pow(x, integer(2))), sub(z, pow(x, integer(3)))};
    vec_sym vars4 = {x, y, z};
    GroebnerOptions opt4;
    opt4.order = MonomialOrder::Lex;
    auto res4 = groebner_basis(polys4, vars4, opt4);
    REQUIRE(res4.status == GroebnerStatus::Success);
    REQUIRE(res4.basis.size() == 4);
    REQUIRE(eq(*res4.basis[0], *sub(pow(x, integer(2)), y)));
    REQUIRE(eq(*res4.basis[1], *sub(mul(x, y), z)));
    REQUIRE(eq(*res4.basis[2], *sub(mul(x, z), pow(y, integer(2)))));
    REQUIRE(eq(*res4.basis[3], *sub(pow(y, integer(3)), pow(z, integer(2)))));

    // Empty and unit ideals:
    auto res_empty = groebner_basis({}, vars4);
    REQUIRE(res_empty.basis.empty());
    REQUIRE(res_empty.status == GroebnerStatus::Success);

    auto res_unit = groebner_basis({integer(1)}, vars4);
    REQUIRE(res_unit.basis.size() == 1);
    REQUIRE(res_unit.basis[0]->__str__() == "1");
}

TEST_CASE("Symbolic Parameters and Sage/symcse Acceptance Tests", "[groebner]")
{
    RCP<const Symbol> x = symbol("x");
    RCP<const Symbol> y = symbol("y");
    RCP<const Symbol> K1 = symbol("K1");
    RCP<const Symbol> K2 = symbol("K2");

    // Variables {x, y}, parameters {K1, K2}:
    // input {x^2 + K1*x*y, x*y + 2*y^3 - K2}
    vec_basic polys = {add(pow(x, integer(2)), mul(K1, mul(x, y))),
                       sub(add(mul(x, y), mul(integer(2), pow(y, integer(3)))), K2)};
    vec_sym vars = {x, y};
    auto res = groebner_basis(polys, vars);
    REQUIRE(res.status == GroebnerStatus::Success);
    // Should compute a successful basis over symbolic coefficient domain!
    REQUIRE(!res.basis.empty());
}

TEST_CASE("FGLM and Zero-Dimensional Conversion", "[groebner]")
{
    RCP<const Symbol> x = symbol("x");
    RCP<const Symbol> y = symbol("y");
    RCP<const Symbol> z = symbol("z");

    // For F = {x^2 - x - 3*y + 1, -2*x + y^2 + y - 1}
    vec_basic polys = {sub(add(sub(pow(x, integer(2)), x), mul(integer(-3), y)), integer(-1)),
                       sub(add(mul(integer(-2), x), pow(y, integer(2))), sub(integer(1), y))};
    vec_sym vars = {x, y};
    
    // First, compute degrevlex basis
    GroebnerOptions opt;
    opt.order = MonomialOrder::DegRevLex;
    auto degrevlex_res = groebner_basis(polys, vars, opt);
    REQUIRE(degrevlex_res.status == GroebnerStatus::Success);

    // Convert to lex using FGLM
    auto lex_res = fglm_convert(degrevlex_res, MonomialOrder::Lex, opt);
    REQUIRE(lex_res.status == GroebnerStatus::Success);
    REQUIRE(!lex_res.basis.empty());
}

TEST_CASE("Cooperative Cancellation Token and Resources", "[groebner]")
{
    RCP<const Symbol> x = symbol("x");
    RCP<const Symbol> y = symbol("y");
    vec_basic polys = {sub(pow(x, integer(2)), y), sub(pow(y, integer(2)), x)};
    vec_sym vars = {x, y};

    GroebnerCancellationToken token;
    token.cancel();

    GroebnerOptions opt;
    opt.cancellation_token = &token;

    auto res = groebner_basis(polys, vars, opt);
    REQUIRE(res.status == GroebnerStatus::Cancelled);
}

TEST_CASE("Polynomial System Solving via Univariate Root", "[groebner]")
{
    RCP<const Symbol> x = symbol("x");
    RCP<const Symbol> y = symbol("y");

    // {y - x^2, y + x^2} gives {(0, 0)}
    vec_basic polys = {sub(y, pow(x, integer(2))), add(y, pow(x, integer(2)))};
    vec_sym vars = {x, y};

    auto sol_set = solve_poly_system_via_univariate_root(polys, vars, x, symbol("zaug"));
    REQUIRE(is_a<FiniteSet>(*sol_set));
    const FiniteSet &fset = down_cast<const FiniteSet &>(*sol_set);
    REQUIRE(fset.get_container().size() == 1);
    
    auto point = *fset.get_container().begin();
    REQUIRE(is_a<Tuple>(*point));
    const Tuple &t = down_cast<const Tuple &>(*point);
    REQUIRE(t.get_args().size() == 2);
    REQUIRE(t.get_args()[0]->__str__() == "0");
    REQUIRE(t.get_args()[1]->__str__() == "0");
}
