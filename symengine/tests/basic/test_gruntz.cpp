#include "catch.hpp"
#include <symengine/basic.h>
#include <symengine/symbol.h>
#include <symengine/integer.h>
#include <symengine/rational.h>
#include <symengine/pow.h>
#include <symengine/functions.h>
#include <symengine/constants.h>
#include <symengine/infinity.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/gruntz.h>

using SymEngine::Basic;
using SymEngine::Symbol;
using SymEngine::symbol;
using SymEngine::Integer;
using SymEngine::integer;
using SymEngine::Rational;
using SymEngine::one;
using SymEngine::zero;
using SymEngine::minus_one;
using SymEngine::Inf;
using SymEngine::NegInf;
using SymEngine::exp;
using SymEngine::log;
using SymEngine::sin;
using SymEngine::cos;
using SymEngine::add;
using SymEngine::sub;
using SymEngine::mul;
using SymEngine::div;
using SymEngine::pow;
using SymEngine::eq;
using SymEngine::gruntz;
using SymEngine::RCP;

TEST_CASE("Gruntz limits at infinity", "[gruntz]")
{
    RCP<const Symbol> x = symbol("x");

    // 1. Polynomial/Rational limits
    // lim_{x -> oo} x = oo
    REQUIRE(eq(*gruntz(x, x, Inf), *Inf));
    // lim_{x -> oo} x^2 = oo
    REQUIRE(eq(*gruntz(pow(x, integer(2)), x, Inf), *Inf));
    // lim_{x -> oo} 1/x = 0
    REQUIRE(eq(*gruntz(div(one, x), x, Inf), *zero));

    // 2. Exponential limits
    // lim_{x -> oo} exp(x) = oo
    REQUIRE(eq(*gruntz(exp(x), x, Inf), *Inf));
    // lim_{x -> oo} exp(-x) = 0
    REQUIRE(eq(*gruntz(exp(neg(x)), x, Inf), *zero));

    // 3. Logarithmic limits
    // lim_{x -> oo} log(x) = oo
    REQUIRE(eq(*gruntz(log(x), x, Inf), *Inf));
    // lim_{x -> oo} log(x)/x = 0
    REQUIRE(eq(*gruntz(div(log(x), x), x, Inf), *zero));

    // 4. Mixed limits
    // lim_{x -> oo} x * exp(-x) = 0
    std::cout << "DEBUG x*exp(-x) limit: " << gruntz(mul(x, exp(neg(x))), x, Inf)->__str__() << std::endl;
    REQUIRE(eq(*gruntz(mul(x, exp(neg(x))), x, Inf), *zero));
    // lim_{x -> oo} exp(x)/x = oo
    REQUIRE(eq(*gruntz(div(exp(x), x), x, Inf), *Inf));
    // lim_{x -> oo} (exp(x) - exp(-x)) / (exp(x) + exp(-x)) = 1
    RCP<const Basic> tanh_exp = div(sub(exp(x), exp(neg(x))), add(exp(x), exp(neg(x))));
    REQUIRE(eq(*gruntz(tanh_exp, x, Inf), *one));

    // 5. Nested exponentials
    // lim_{x -> oo} exp(exp(x)) = oo
    REQUIRE(eq(*gruntz(exp(exp(x)), x, Inf), *Inf));

    // 6. Rational exponents
    // lim_{x -> oo} x^(1/2) = oo
    // In SymEngine, 1/2 is constructed as Rational::from_two_ints(1, 2)
    RCP<const Basic> half = Rational::from_two_ints(1, 2);
    REQUIRE(eq(*gruntz(pow(x, half), x, Inf), *Inf));
}

TEST_CASE("Gruntz limits at -infinity and finite points", "[gruntz]")
{
    RCP<const Symbol> x = symbol("x");

    // 7. Limits at -oo
    // lim_{x -> -oo} exp(x) = 0
    REQUIRE(eq(*gruntz(exp(x), x, NegInf), *zero));
    // lim_{x -> -oo} exp(-x) = oo
    REQUIRE(eq(*gruntz(exp(neg(x)), x, NegInf), *Inf));

    // 8. Limits at finite points
    // lim_{x -> 0+} 1/x = oo
    REQUIRE(eq(*gruntz(div(one, x), x, zero, "+"), *Inf));
    // lim_{x -> 0-} 1/x = -oo
    REQUIRE(eq(*gruntz(div(one, x), x, zero, "-"), *NegInf));

    // lim_{x -> 2+} (x - 2) = 0
    REQUIRE(eq(*gruntz(sub(x, integer(2)), x, integer(2), "+"), *zero));
}

TEST_CASE("Gruntz with trigonometric functions", "[gruntz]")
{
    RCP<const Symbol> x = symbol("x");

    // 9. Trigonometric limit
    // lim_{x -> 0+} sin(x)/x = 1
    REQUIRE(eq(*gruntz(div(sin(x), x), x, zero, "+"), *one));
}
