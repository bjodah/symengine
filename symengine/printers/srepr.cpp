#include <symengine/printers/srepr.h>

namespace SymEngine
{

//     namespace /* anonymous */ {
//         std::string get_name_of_type(const Basic &x) {
//             TypeID type_code = x.get_type_code();
//             switch (type_code) {
// #define SYMENGINE_ENUM(type, Class)             \
//                 case type:                      \
//                     return # Class ;
// #include "symengine/type_codes.inc"
// #undef SYMENGINE_ENUM
//             default:
//                 throw SymEngineException("Not supported: " + typeName<Basic>(x));
//             }
//             //std::unreachable();
//         }
//     }

void SreprPrinter::bvisit(const Basic &x)
{
//    throw SymEngineException("Not supported: " + typeName<Basic>(x));
    std::ostringstream s;
    //s << get_name_of_type(x);
    s << type_code_name(x.get_type_code());
    auto const args = x.get_args();
    if (is_a_sub<Constant>(x)) {
        SYMENGINE_ASSERT(args.size() == 0);
        s << '.' << down_cast<Constant const&>(x).get_name();
    } else {
        s << '(';
        if (is_a<Symbol>(x)) {
            s << '"' << down_cast<Symbol const&>(x).get_name() << '"';
        } else if (args.size() == 0 && is_a_sub<Number>(x)) {
            s << x.__str__();
        } else {
            bool first {true};
            for (auto const& arg : args) {
                if (!first) {
                    s << ", ";
                } else {
                    first = false;
                }
                s << apply(arg);
            }
        }
        s << ')';
    }
    str_ = s.str();
}

std::string SreprPrinter::apply(const RCP<const Basic> &b)
{
    b->accept(*this);
    return str_;
}

}
