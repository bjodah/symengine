#include "catch.hpp"

#include <symengine/printers/srepr.h>
#include <symengine/parser.h>

namespace se = SymEngine;

using se::E;
using se::Inf;
using se::pi;
using se::SreprPrinter;

class NamedFunctionWrapper : public se::FunctionWrapper
{
public:
    NamedFunctionWrapper(std::string name, const se::vec_basic &args)
        : FunctionWrapper(std::move(name), args)
    {
    }

    se::RCP<const se::Number> eval(long bits) const override
    {
        return se::zero;
    }

    se::RCP<const se::Basic> create(const se::vec_basic &args) const override
    {
        return se::make_rcp<NamedFunctionWrapper>(get_name(), args);
    }

    se::RCP<const se::Basic>
    diff_impl(const se::RCP<const se::Symbol> &symbol) const override
    {
        return se::zero;
    }
};

TEST_CASE("Simple expressions", "[SreprPrinter]")
{
    SreprPrinter srp;
    REQUIRE(srp.apply(Inf) == "Infty(Integer(1))");
    REQUIRE(srp.apply(E) == "Constant.E");
    REQUIRE(srp.apply(pi) == "Constant.pi");
    REQUIRE(srp.apply(se::parse("sin(x**y)"))
            == "Sin(Pow(Symbol(\"x\"), Symbol(\"y\")))");
    REQUIRE(srp.apply(se::parse("x*y+z"))
            == "Add(Symbol(\"z\"), Mul(Symbol(\"x\"), Symbol(\"y\")))");
}

TEST_CASE("Function identity and escaped names", "[SreprPrinter]")
{
    SreprPrinter srp;
    auto x = se::symbol("x\"\\\n");
    auto y = se::symbol("y");
    auto function = se::function_symbol("f\"\\\n", {x, y});
    auto first = se::make_rcp<NamedFunctionWrapper>("First", se::vec_basic{x});
    auto second
        = se::make_rcp<NamedFunctionWrapper>("Second", se::vec_basic{x});
    auto dummy_first = se::dummy("same", 41);
    auto dummy_second = se::dummy("same", 42);

    REQUIRE(srp.apply(x) == "Symbol(\"x\\\"\\\\\\n\")");
    REQUIRE(srp.apply(function)
            == "FunctionSymbol(\"f\\\"\\\\\\n\", "
               "Symbol(\"x\\\"\\\\\\n\"), Symbol(\"y\"))");
    REQUIRE(srp.apply(first)
            == "FunctionWrapper(\"First\", Symbol(\"x\\\"\\\\\\n\"))");
    REQUIRE(srp.apply(second)
            == "FunctionWrapper(\"Second\", Symbol(\"x\\\"\\\\\\n\"))");
    REQUIRE(srp.apply(first) != srp.apply(second));
    REQUIRE(srp.apply(first) == srp.apply(first));
    REQUIRE(srp.apply(dummy_first) == "Dummy(\"same\", 41)");
    REQUIRE(srp.apply(dummy_second) == "Dummy(\"same\", 42)");
    REQUIRE(srp.apply(dummy_first) != srp.apply(dummy_second));
}

TEST_CASE("Commutative argument order is stable", "[SreprPrinter]")
{
    SreprPrinter srp;
    auto x = se::symbol("x");
    auto y = se::symbol("y");

    REQUIRE(srp.apply(se::add(x, y)) == srp.apply(se::add(y, x)));
    REQUIRE(srp.apply(se::mul(x, y)) == srp.apply(se::mul(y, x)));
}
