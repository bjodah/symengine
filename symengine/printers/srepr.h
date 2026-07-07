#pragma once
// -*- mode: c++ -*-
#include <symengine/visitor.h>
#include <symengine/printers.h>

namespace SymEngine
{

class SreprPrinter : public BaseVisitor<SreprPrinter>
{
protected:
    std::string str_;

public:
    void bvisit(const Basic &x);
    void bvisit(const Symbol &x);
    void bvisit(const Dummy &x);
    void bvisit(const FunctionSymbol &x);
    void bvisit(const FunctionWrapper &x);

    std::string apply(const RCP<const Basic> &x);

private:
    void print_function(const char *type_name, const FunctionSymbol &x);
};
} // namespace SymEngine
