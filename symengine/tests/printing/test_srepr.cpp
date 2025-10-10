#include "catch.hpp"

#include <symengine/printers/srepr.h>
#include <symengine/parser.h>

namespace se = SymEngine;

using se::SreprPrinter;
using se::Inf;
using se::E;
using se::pi;

TEST_CASE("Simple expressions", "[SreprPrinter]")
{
    SreprPrinter srp;
    REQUIRE(srp.apply(Inf) == "Infty(Integer(1))");
    REQUIRE(srp.apply(E) == "Constant.E");
    REQUIRE(srp.apply(pi) == "Constant.pi");
    REQUIRE(srp.apply(se::parse("sin(x**y)")) == "Sin(Pow(Symbol(\"x\"), Symbol(\"y\")))");
    REQUIRE(srp.apply(se::parse("x*y+z")) == "Add(Symbol(\"z\"), Mul(Symbol(\"x\"), Symbol(\"y\")))");
}
