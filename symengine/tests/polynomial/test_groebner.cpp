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

// ---------------------------------------------------------------------------
// Phase 0 — Monomial ordering and arithmetic
// ---------------------------------------------------------------------------

TEST_CASE("Monomial Ordering Tests", "[groebner]")
{
    // lex((1,2,3)) == (1,2,3)
    Monomial m1 = {1, 2, 3};
    Monomial m2 = {1, 2, 3};
    REQUIRE(!compare_monomial_less(m1, m2, MonomialOrder::Lex));
    REQUIRE(!compare_monomial_less(m2, m1, MonomialOrder::Lex));

    // (0, 1, 1) > (0, 0, 2) under grlex and degrevlex
    Monomial m3 = {0, 0, 2};
    Monomial m4 = {0, 1, 1};
    REQUIRE(compare_monomial_less(m3, m4, MonomialOrder::GrLex));
    REQUIRE(compare_monomial_less(m3, m4, MonomialOrder::DegRevLex));

    // (0, 3, 1) < (2, 2, 1) under grlex and degrevlex
    Monomial m5 = {0, 3, 1};
    Monomial m6 = {2, 2, 1};
    REQUIRE(compare_monomial_less(m5, m6, MonomialOrder::GrLex));
    REQUIRE(compare_monomial_less(m5, m6, MonomialOrder::DegRevLex));

    // With generator order [z, y, x], sort:
    // {x, x^2*z^2, x*y, x^2, 1, y^2, x^3, y, z, x*y^2*z, x^2*y^2}
    std::vector<Monomial> list = {
        {0, 0, 1}, {2, 0, 2}, {0, 1, 1}, {0, 0, 2}, {0, 0, 0},
        {0, 2, 0}, {0, 0, 3}, {0, 1, 0}, {1, 0, 0}, {1, 2, 1}, {0, 2, 2}
    };

    std::vector<Monomial> lex_expected = {
        {0, 0, 0}, {0, 0, 1}, {0, 0, 2}, {0, 0, 3}, {0, 1, 0},
        {0, 1, 1}, {0, 2, 0}, {0, 2, 2}, {1, 0, 0}, {1, 2, 1}, {2, 0, 2}
    };
    std::vector<Monomial> lex_actual = list;
    std::sort(lex_actual.begin(), lex_actual.end(), [](const Monomial &a, const Monomial &b) {
        return compare_monomial_less(a, b, MonomialOrder::Lex);
    });
    REQUIRE(lex_actual == lex_expected);

    std::vector<Monomial> grlex_expected = {
        {0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {1, 0, 0}, {0, 0, 2},
        {0, 1, 1}, {0, 2, 0}, {0, 0, 3}, {0, 2, 2}, {1, 2, 1}, {2, 0, 2}
    };
    std::vector<Monomial> grlex_actual = list;
    std::sort(grlex_actual.begin(), grlex_actual.end(), [](const Monomial &a, const Monomial &b) {
        return compare_monomial_less(a, b, MonomialOrder::GrLex);
    });
    REQUIRE(grlex_actual == grlex_expected);

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

TEST_CASE("Monomial Arithmetic Helpers", "[groebner]")
{
    Monomial a = {1, 2, 0};
    Monomial b = {2, 2, 1};
    Monomial c = {0, 1, 0};

    REQUIRE(degree(a) == 3);
    REQUIRE(degree(Monomial{0, 0, 0}) == 0);

    REQUIRE(divides(a, b));
    REQUIRE(!divides(b, a));
    REQUIRE(divides(Monomial{0, 0, 0}, a));

    Monomial q = quotient(b, a);
    REQUIRE(q == Monomial{1, 0, 1});

    Monomial p = multiply(a, c);
    REQUIRE(p == Monomial{1, 3, 0});

    Monomial l = lcm(a, b);
    REQUIRE(l == Monomial{2, 2, 1});
    Monomial l2 = lcm(Monomial{1, 0, 2}, Monomial{0, 3, 1});
    REQUIRE(l2 == Monomial{1, 3, 2});
}

// ---------------------------------------------------------------------------
// Phase 1 — Buchberger over QQ
// ---------------------------------------------------------------------------

namespace {

bool basis_is_canonical(const GroebnerResult &res)
{
    return res.status == GroebnerStatus::Success
           && is_groebner(res.basis, res.variables, res.order)
           && is_reduced_basis(res.basis, res.variables, res.order);
}

bool every_input_in_ideal(const vec_basic &inputs, const GroebnerResult &res)
{
    for (const auto &p : inputs) {
        auto nf = normal_form(p, res.basis, res.variables, res.order);
        if (!eq(*nf, *zero)) return false;
    }
    return true;
}

// Check ideal equality by mutual containment, where A is a Groebner basis
// under order_a and B under order_b — each direction reduces with respect to
// the basis using that basis's own order.
bool same_ideal(const vec_basic &A, MonomialOrder order_a,
                const vec_basic &B, MonomialOrder order_b,
                const vec_sym &vars)
{
    for (const auto &p : A) {
        if (!eq(*normal_form(p, B, vars, order_b), *zero)) return false;
    }
    for (const auto &p : B) {
        if (!eq(*normal_form(p, A, vars, order_a), *zero)) return false;
    }
    return true;
}

} // namespace

TEST_CASE("Rational Buchberger — small exact bases", "[groebner]")
{
    auto x = symbol("x");
    auto y = symbol("y");
    auto z = symbol("z");

    // {x^2 + 2*x*y^2, x*y + 2*y^3 - 1} -> {x, y^3 - 1/2}
    vec_basic polys1 = {add(pow(x, integer(2)), mul(integer(2), mul(x, pow(y, integer(2))))),
                        sub(add(mul(x, y), mul(integer(2), pow(y, integer(3)))), integer(1))};
    vec_sym vars1 = {x, y};
    GroebnerOptions opt1;
    opt1.order = MonomialOrder::Lex;
    auto res1 = groebner_basis(polys1, vars1, opt1);
    REQUIRE(basis_is_canonical(res1));
    REQUIRE(res1.basis.size() == 2);
    REQUIRE(res1.basis[0]->__str__() == "x");
    REQUIRE(res1.basis[1]->__str__() == "-1/2 + y**3");
    REQUIRE(every_input_in_ideal(polys1, res1));

    // {2*x^2*y + y^2, 2*x^3 + x*y - 1} over {y, x} lex -> {y, x^3 - 1/2}
    vec_basic polys2 = {add(mul(integer(2), mul(pow(x, integer(2)), y)), pow(y, integer(2))),
                        sub(add(mul(integer(2), pow(x, integer(3))), mul(x, y)), integer(1))};
    vec_sym vars2 = {y, x};
    auto res2 = groebner_basis(polys2, vars2, opt1);
    REQUIRE(basis_is_canonical(res2));
    REQUIRE(res2.basis.size() == 2);
    REQUIRE(res2.basis[0]->__str__() == "y");
    REQUIRE(res2.basis[1]->__str__() == "-1/2 + x**3");

    // {x^3 - 2*x*y, x^2*y + x - 2*y^2} grlex -> {x^2, x*y, -1/2*x + y^2}
    vec_basic polys3 = {sub(pow(x, integer(3)), mul(integer(2), mul(x, y))),
                        sub(add(mul(pow(x, integer(2)), y), x), mul(integer(2), pow(y, integer(2))))};
    vec_sym vars3 = {x, y};
    GroebnerOptions opt3; opt3.order = MonomialOrder::GrLex;
    auto res3 = groebner_basis(polys3, vars3, opt3);
    REQUIRE(basis_is_canonical(res3));
    REQUIRE(res3.basis.size() == 3);
    REQUIRE(res3.basis[0]->__str__() == "x**2");
    REQUIRE(res3.basis[1]->__str__() == "x*y");
    REQUIRE(res3.basis[2]->__str__() == "(-1/2)*x + y**2");
    REQUIRE(every_input_in_ideal(polys3, res3));

    // {-x^2+y, -x^3+z} lex -> {x^2-y, x*y-z, x*z-y^2, y^3-z^2}
    vec_basic polys4 = {sub(y, pow(x, integer(2))), sub(z, pow(x, integer(3)))};
    vec_sym vars4 = {x, y, z};
    auto res4 = groebner_basis(polys4, vars4, opt1);
    REQUIRE(basis_is_canonical(res4));
    REQUIRE(res4.basis.size() == 4);
    REQUIRE(eq(*res4.basis[0], *sub(pow(x, integer(2)), y)));
    REQUIRE(eq(*res4.basis[1], *sub(mul(x, y), z)));
    REQUIRE(eq(*res4.basis[2], *sub(mul(x, z), pow(y, integer(2)))));
    REQUIRE(eq(*res4.basis[3], *sub(pow(y, integer(3)), pow(z, integer(2)))));
    REQUIRE(every_input_in_ideal(polys4, res4));

    // {-x^2+y, -x^3+z} grlex — verify property (size and identical ideal)
    auto res4g = groebner_basis(polys4, vars4, opt3);
    REQUIRE(basis_is_canonical(res4g));
    REQUIRE(every_input_in_ideal(polys4, res4g));
    REQUIRE(same_ideal(res4.basis, MonomialOrder::Lex,
                       res4g.basis, MonomialOrder::GrLex, vars4));
}

TEST_CASE("Rational Buchberger — additional small exact bases", "[groebner]")
{
    auto x = symbol("x");
    auto y = symbol("y");
    auto z = symbol("z");
    GroebnerOptions lex; lex.order = MonomialOrder::Lex;
    GroebnerOptions grl; grl.order = MonomialOrder::GrLex;

    // {x - y^2, -y^3 + z} lex -> {x - y^2, y^3 - z}
    vec_basic A = {sub(x, pow(y, integer(2))), sub(z, pow(y, integer(3)))};
    auto rA = groebner_basis(A, {x, y, z}, lex);
    REQUIRE(basis_is_canonical(rA));
    REQUIRE(rA.basis.size() == 2);
    REQUIRE(every_input_in_ideal(A, rA));

    // same grlex -> {x^2 - y*z, x*y - z, -x + y^2}
    auto rAg = groebner_basis(A, {x, y, z}, grl);
    REQUIRE(basis_is_canonical(rAg));
    REQUIRE(rAg.basis.size() == 3);
    REQUIRE(same_ideal(rA.basis, MonomialOrder::Lex,
                       rAg.basis, MonomialOrder::GrLex, {x, y, z}));

    // {x - z^2, y - z^3} lex -> {x - z^2, y - z^3}
    vec_basic B = {sub(x, pow(z, integer(2))), sub(y, pow(z, integer(3)))};
    auto rB = groebner_basis(B, {x, y, z}, lex);
    REQUIRE(basis_is_canonical(rB));
    REQUIRE(rB.basis.size() == 2);

    // same grlex -> three generators
    auto rBg = groebner_basis(B, {x, y, z}, grl);
    REQUIRE(basis_is_canonical(rBg));
    REQUIRE(rBg.basis.size() == 3);
    REQUIRE(same_ideal(rB.basis, MonomialOrder::Lex,
                       rBg.basis, MonomialOrder::GrLex, {x, y, z}));

    // {-y^2 + z, x - y^3} lex -> {x - y*z, y^2 - z}
    vec_basic C = {sub(z, pow(y, integer(2))), sub(x, pow(y, integer(3)))};
    auto rC = groebner_basis(C, {x, y, z}, lex);
    REQUIRE(basis_is_canonical(rC));
    REQUIRE(rC.basis.size() == 2);

    // same grlex -> 4 generators
    auto rCg = groebner_basis(C, {x, y, z}, grl);
    REQUIRE(basis_is_canonical(rCg));
    REQUIRE(rCg.basis.size() == 4);
}

TEST_CASE("Rational Buchberger — univariate reductions", "[groebner]")
{
    auto x = symbol("x");
    GroebnerOptions opt; opt.order = MonomialOrder::Lex;

    // {x^3 - 1, x^2 - 1} -> {x - 1}
    vec_basic p1 = {sub(pow(x, integer(3)), integer(1)),
                    sub(pow(x, integer(2)), integer(1))};
    auto r1 = groebner_basis(p1, {x}, opt);
    REQUIRE(basis_is_canonical(r1));
    REQUIRE(r1.basis.size() == 1);
    REQUIRE(eq(*r1.basis[0], *sub(x, integer(1))));

    // {x^2 - 1, x^3 + 1} -> {x + 1}
    vec_basic p2 = {sub(pow(x, integer(2)), integer(1)),
                    add(pow(x, integer(3)), integer(1))};
    auto r2 = groebner_basis(p2, {x}, opt);
    REQUIRE(basis_is_canonical(r2));
    REQUIRE(r2.basis.size() == 1);
    REQUIRE(eq(*r2.basis[0], *add(x, integer(1))));
}

TEST_CASE("Rational Buchberger — empty, unit and zero generator", "[groebner]")
{
    auto x = symbol("x");
    auto y = symbol("y");

    // Empty -> empty
    auto rE = groebner_basis({}, {x, y});
    REQUIRE(rE.status == GroebnerStatus::Success);
    REQUIRE(rE.basis.empty());

    // {1} -> {1}
    auto rU = groebner_basis({integer(1)}, {x, y});
    REQUIRE(rU.status == GroebnerStatus::Success);
    REQUIRE(rU.basis.size() == 1);
    REQUIRE(eq(*rU.basis[0], *integer(1)));

    // Inconsistent {y - x, y - x - 1} -> {1}
    vec_basic incon = {sub(y, x), sub(sub(y, x), integer(1))};
    auto rI = groebner_basis(incon, {x, y});
    REQUIRE(rI.status == GroebnerStatus::Success);
    REQUIRE(rI.basis.size() == 1);
    REQUIRE(eq(*rI.basis[0], *integer(1)));

    // Adding the zero polynomial does not change the basis.
    vec_basic ps = {sub(pow(x, integer(2)), y), sub(pow(y, integer(2)), x)};
    auto r_ref = groebner_basis(ps, {x, y});
    vec_basic ps_with_zero = ps;
    ps_with_zero.push_back(integer(0));
    auto r_zero = groebner_basis(ps_with_zero, {x, y});
    REQUIRE(same_ideal(r_ref.basis, MonomialOrder::DegRevLex,
                       r_zero.basis, MonomialOrder::DegRevLex, {x, y}));
}

// ---------------------------------------------------------------------------
// Property tests: normal_form, is_groebner, is_reduced_basis
// ---------------------------------------------------------------------------

TEST_CASE("normal_form and ideal membership", "[groebner]")
{
    auto x = symbol("x");
    auto y = symbol("y");

    // normal_form(2*x^4 + y^2 - x^2 + y^3, {x^3 - x, y^3 - y}) == x^2 + y^2 + y
    auto p = add(add(mul(integer(2), pow(x, integer(4))), pow(y, integer(2))),
                 sub(pow(y, integer(3)), pow(x, integer(2))));
    vec_basic G = {sub(pow(x, integer(3)), x), sub(pow(y, integer(3)), y)};
    auto r = normal_form(p, G, {x, y}, MonomialOrder::GrLex);
    auto expected = add(add(pow(x, integer(2)), pow(y, integer(2))), y);
    REQUIRE(eq(*r, *expected));

    // 2*x^3 + y^3 + 3*y is in ideal of {x^2 + y^2 - 1, x*y - 2}.
    auto e = add(add(mul(integer(2), pow(x, integer(3))), pow(y, integer(3))),
                 mul(integer(3), y));
    vec_basic G2 = {sub(add(pow(x, integer(2)), pow(y, integer(2))), integer(1)),
                    sub(mul(x, y), integer(2))};
    auto gb = groebner_basis(G2, {x, y});
    REQUIRE(eq(*normal_form(e, gb.basis, {x, y}, gb.order), *zero));

    // e + 1 is not in the ideal.
    auto e1 = add(e, integer(1));
    REQUIRE(!eq(*normal_form(e1, gb.basis, {x, y}, gb.order), *zero));
}

TEST_CASE("is_groebner and is_reduced_basis", "[groebner]")
{
    auto x = symbol("x");
    auto y = symbol("y");
    auto half = Rational::from_two_ints(-1, 2);

    // {x^2, x*y, -1/2*x + y^2} is a Groebner basis (grlex).
    vec_basic G_ok = {pow(x, integer(2)), mul(x, y),
                      add(mul(half, x), pow(y, integer(2)))};
    REQUIRE(is_groebner(G_ok, {x, y}, MonomialOrder::GrLex));
    REQUIRE(is_reduced_basis(G_ok, {x, y}, MonomialOrder::GrLex));

    // {x^3, x*y, -1/2*x + y^2} is NOT a Groebner basis: S(x*y, -x/2 + y^2)
    // reduces to -y^3 which is not divisible by x^3 or x*y.
    vec_basic G_bad = {pow(x, integer(3)), mul(x, y),
                       add(mul(half, x), pow(y, integer(2)))};
    REQUIRE(!is_groebner(G_bad, {x, y}, MonomialOrder::GrLex));

    // {2*x^2 + y^2, x^2 + y^2} is not reduced (not monic, and LM of one
    // divides a term of the other).
    vec_basic G_not_reduced = {
        add(mul(integer(2), pow(x, integer(2))), pow(y, integer(2))),
        add(pow(x, integer(2)), pow(y, integer(2)))
    };
    REQUIRE(!is_reduced_basis(G_not_reduced, {x, y}, MonomialOrder::GrLex));

    // The reduced basis of the same ideal is reduced.
    auto reduced = groebner_basis(G_not_reduced, {x, y});
    REQUIRE(is_reduced_basis(reduced.basis, {x, y}, MonomialOrder::DegRevLex));
}

TEST_CASE("Toric/binomial ideal — Buchberger property tests", "[groebner]")
{
    auto w = symbol("w");
    auto x = symbol("x");
    auto y = symbol("y");
    auto z = symbol("z");

    vec_basic input = {
        sub(pow(w, integer(2)), mul(x, z)),
        sub(mul(w, x), mul(y, z)),
        sub(pow(x, integer(2)), mul(w, y)),
        sub(mul(x, y), pow(z, integer(2))),
        sub(pow(y, integer(2)), mul(w, z)),
    };
    vec_sym vars = {w, x, y, z};
    GroebnerOptions opt; opt.order = MonomialOrder::DegRevLex;
    auto res = groebner_basis(input, vars, opt);

    REQUIRE(res.status == GroebnerStatus::Success);
    REQUIRE(is_groebner(res.basis, vars, MonomialOrder::DegRevLex));
    REQUIRE(is_reduced_basis(res.basis, vars, MonomialOrder::DegRevLex));
    REQUIRE(every_input_in_ideal(input, res));
}

// ---------------------------------------------------------------------------
// Symbolic parameters
// ---------------------------------------------------------------------------

TEST_CASE("Symbolic Parameters — small parameter system", "[groebner]")
{
    auto x = symbol("x");
    auto y = symbol("y");
    auto K1 = symbol("K1");
    auto K2 = symbol("K2");
    auto C1 = symbol("C1");

    // {C1*x + y, x - 1} over {x, y}: -> {x - 1, y + C1}
    vec_basic ps_lin = {add(mul(C1, x), y), sub(x, integer(1))};
    GroebnerOptions opt; opt.order = MonomialOrder::Lex;
    auto rL = groebner_basis(ps_lin, {x, y}, opt);
    REQUIRE(rL.status == GroebnerStatus::Success);
    REQUIRE(rL.basis.size() == 2);
    // Verify each input reduces to zero.
    REQUIRE(every_input_in_ideal(ps_lin, rL));

    // {x^2 + K1*x*y, x*y + 2*y^3 - K2} over {x, y}
    vec_basic ps = {add(pow(x, integer(2)), mul(K1, mul(x, y))),
                    sub(add(mul(x, y), mul(integer(2), pow(y, integer(3)))), K2)};
    auto r = groebner_basis(ps, {x, y});
    REQUIRE(r.status == GroebnerStatus::Success);
    REQUIRE(!r.basis.empty());
    REQUIRE(every_input_in_ideal(ps, r));
}

// ---------------------------------------------------------------------------
// FGLM and zero-dimensional handling
// ---------------------------------------------------------------------------

TEST_CASE("FGLM — degrevlex to lex round-trip", "[groebner]")
{
    auto x = symbol("x");
    auto y = symbol("y");

    // F = {x^2 - x - 3y + 1, -2x + y^2 + y - 1}
    vec_basic polys = {sub(add(sub(pow(x, integer(2)), x), mul(integer(-3), y)), integer(-1)),
                       sub(add(mul(integer(-2), x), pow(y, integer(2))), sub(integer(1), y))};
    vec_sym vars = {x, y};

    GroebnerOptions opt; opt.order = MonomialOrder::DegRevLex;
    auto degrevlex_res = groebner_basis(polys, vars, opt);
    REQUIRE(basis_is_canonical(degrevlex_res));

    auto lex_res = fglm_convert(degrevlex_res, MonomialOrder::Lex, opt);
    REQUIRE(lex_res.status == GroebnerStatus::Success);
    REQUIRE(!lex_res.basis.empty());
    REQUIRE(is_groebner(lex_res.basis, vars, MonomialOrder::Lex));
    REQUIRE(is_reduced_basis(lex_res.basis, vars, MonomialOrder::Lex));

    // Round-trip: convert back to degrevlex — should be the same ideal.
    GroebnerOptions opt2; opt2.order = MonomialOrder::Lex;
    auto back = fglm_convert(lex_res, MonomialOrder::DegRevLex, opt2);
    REQUIRE(back.status == GroebnerStatus::Success);
    REQUIRE(same_ideal(degrevlex_res.basis, MonomialOrder::DegRevLex,
                       back.basis, MonomialOrder::DegRevLex, vars));
}

TEST_CASE("FGLM — non-zero-dimensional rejection", "[groebner]")
{
    auto x = symbol("x");
    auto y = symbol("y");
    auto z = symbol("z");

    // {x*y - 1} over {x, y} is not zero-dimensional.
    vec_basic p1 = {sub(mul(x, y), integer(1))};
    GroebnerOptions opt; opt.order = MonomialOrder::DegRevLex;
    auto g1 = groebner_basis(p1, {x, y}, opt);
    auto c1 = fglm_convert(g1, MonomialOrder::Lex, opt);
    REQUIRE(c1.status == GroebnerStatus::NotZeroDimensional);

    // {x^3 + y^2} over {x, y} is not zero-dimensional.
    vec_basic p2 = {add(pow(x, integer(3)), pow(y, integer(2)))};
    auto g2 = groebner_basis(p2, {x, y}, opt);
    auto c2 = fglm_convert(g2, MonomialOrder::Lex, opt);
    REQUIRE(c2.status == GroebnerStatus::NotZeroDimensional);

    // {x, y, z} over {x, y, z, t}: not zero-dimensional (t is unconstrained).
    auto t = symbol("t");
    vec_basic p3 = {x, y, z};
    auto g3 = groebner_basis(p3, {x, y, z, t}, opt);
    auto c3 = fglm_convert(g3, MonomialOrder::Lex, opt);
    REQUIRE(c3.status == GroebnerStatus::NotZeroDimensional);

    // {x, y, z} over {x, y, z}: zero-dimensional. Trivial conversion.
    auto g4 = groebner_basis({x, y, z}, {x, y, z}, opt);
    auto c4 = fglm_convert(g4, MonomialOrder::Lex, opt);
    REQUIRE(c4.status == GroebnerStatus::Success);
}

TEST_CASE("FGLM — cancellation", "[groebner]")
{
    auto x = symbol("x");
    auto y = symbol("y");
    vec_basic ps = {sub(pow(x, integer(2)), y), sub(pow(y, integer(2)), x)};
    GroebnerOptions opt; opt.order = MonomialOrder::DegRevLex;
    auto g = groebner_basis(ps, {x, y}, opt);

    GroebnerCancellationToken token;
    token.cancel();
    GroebnerOptions opt2 = opt;
    opt2.cancellation_token = &token;
    auto c = fglm_convert(g, MonomialOrder::Lex, opt2);
    REQUIRE(c.status == GroebnerStatus::Cancelled);
}

// ---------------------------------------------------------------------------
// Cancellation and resource limits
// ---------------------------------------------------------------------------

TEST_CASE("Cancellation token short-circuits Buchberger", "[groebner]")
{
    auto x = symbol("x");
    auto y = symbol("y");
    vec_basic polys = {sub(pow(x, integer(2)), y), sub(pow(y, integer(2)), x)};

    GroebnerCancellationToken token;
    token.cancel();

    GroebnerOptions opt;
    opt.cancellation_token = &token;
    auto res = groebner_basis(polys, {x, y}, opt);
    REQUIRE(res.status == GroebnerStatus::Cancelled);
}

TEST_CASE("Resource limits exit cleanly", "[groebner]")
{
    auto w = symbol("w");
    auto x = symbol("x");
    auto y = symbol("y");
    auto z = symbol("z");
    // Toric/binomial system has many overlapping S-pairs.
    vec_basic polys = {
        sub(pow(w, integer(2)), mul(x, z)),
        sub(mul(w, x), mul(y, z)),
        sub(pow(x, integer(2)), mul(w, y)),
        sub(mul(x, y), pow(z, integer(2))),
        sub(pow(y, integer(2)), mul(w, z)),
    };
    {
        GroebnerOptions opt; opt.max_s_pairs = 1;
        auto res = groebner_basis(polys, {w, x, y, z}, opt);
        REQUIRE(res.status == GroebnerStatus::ResourceLimitExceeded);
    }
    {
        GroebnerOptions opt; opt.max_reduction_steps = 1;
        auto res = groebner_basis(polys, {w, x, y, z}, opt);
        REQUIRE(res.status == GroebnerStatus::ResourceLimitExceeded);
    }
}

// ---------------------------------------------------------------------------
// Polynomial system solving — first-shape solver
// ---------------------------------------------------------------------------

TEST_CASE("Solver — augmented zaug shape, simple cases", "[groebner]")
{
    auto x = symbol("x");
    auto y = symbol("y");
    auto zaug = symbol("zaug");

    // {y - x^2, y + x^2} -> {(0, 0)}
    vec_basic polys = {sub(y, pow(x, integer(2))), add(y, pow(x, integer(2)))};
    auto sol = solve_poly_system_via_univariate_root(polys, {x, y}, x, zaug);
    REQUIRE(is_a<FiniteSet>(*sol));
    REQUIRE(down_cast<const FiniteSet &>(*sol).get_container().size() == 1);

    // Single-variable: {x - 1} over {x} -> {(1,)}
    vec_basic pp = {sub(x, integer(1))};
    auto sol_uni = solve_poly_system_via_univariate_root(pp, {x}, x, zaug);
    REQUIRE(is_a<FiniteSet>(*sol_uni));
    REQUIRE(down_cast<const FiniteSet &>(*sol_uni).get_container().size() == 1);
}

TEST_CASE("Solver — linear systems", "[groebner]")
{
    auto x = symbol("x");
    auto y = symbol("y");
    auto z = symbol("z");

    // {2x-3, 3y/2 - 2x, z-5y} -> (3/2, 2, 10)
    vec_basic lin = {sub(mul(integer(2), x), integer(3)),
                     sub(mul(Rational::from_two_ints(3, 2), y), mul(integer(2), x)),
                     sub(z, mul(integer(5), y))};
    auto sol = solve_poly_system(lin, {x, y, z});
    REQUIRE(is_a<FiniteSet>(*sol));
    const auto &cont = down_cast<const FiniteSet &>(*sol).get_container();
    REQUIRE(cont.size() == 1);
    auto pt = *cont.begin();
    REQUIRE(is_a<Tuple>(*pt));
    const auto &args = down_cast<const Tuple &>(*pt).get_args();
    REQUIRE(args.size() == 3);
    REQUIRE(eq(*args[0], *Rational::from_two_ints(3, 2)));
    REQUIRE(eq(*args[1], *integer(2)));
    REQUIRE(eq(*args[2], *integer(10)));
}

TEST_CASE("Solver — inconsistent system", "[groebner]")
{
    auto x = symbol("x");
    auto y = symbol("y");
    auto zaug = symbol("zaug");
    vec_basic incon = {sub(y, x), sub(sub(y, x), integer(1))};
    auto sol = solve_poly_system_via_univariate_root(incon, {x, y}, x, zaug);
    REQUIRE(is_a<EmptySet>(*sol));
}

TEST_CASE("Solver — augmented zaug, x^2-y / x^3-x", "[groebner]")
{
    auto x = symbol("x");
    auto y = symbol("y");
    auto zaug = symbol("zaug");

    // {x^2 - y, x^3 - x} has solutions (0,0), (1,1), (-1,1)
    vec_basic ps = {sub(pow(x, integer(2)), y), sub(pow(x, integer(3)), x)};
    auto sol = solve_poly_system_via_univariate_root(ps, {x, y}, x, zaug);
    REQUIRE(is_a<FiniteSet>(*sol));
    REQUIRE(down_cast<const FiniteSet &>(*sol).get_container().size() == 3);
}

// ---------------------------------------------------------------------------
// Statistics
// ---------------------------------------------------------------------------

TEST_CASE("Stats — input_polys, output_polys, reductions_to_zero are populated", "[groebner]")
{
    auto x = symbol("x");
    auto y = symbol("y");
    vec_basic polys = {sub(pow(x, integer(2)), y), sub(pow(y, integer(2)), x)};
    auto r = groebner_basis(polys, {x, y});
    REQUIRE(r.status == GroebnerStatus::Success);
    REQUIRE(r.stats.input_polys == 2);
    REQUIRE(r.stats.output_polys == r.basis.size());
    REQUIRE(r.stats.max_basis_size >= r.stats.output_polys);
}
