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

    std::string apply(const RCP<const Basic> &x);
};
}
