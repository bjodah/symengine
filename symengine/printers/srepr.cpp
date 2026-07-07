#include <symengine/printers/srepr.h>

#include <algorithm>

namespace SymEngine
{

namespace
{
std::string quoted(const std::string &value)
{
    static const char hex[] = "0123456789abcdef";
    std::string result{"\""};
    for (unsigned char character : value) {
        switch (character) {
            case '\\':
                result += "\\\\";
                break;
            case '"':
                result += "\\\"";
                break;
            case '\b':
                result += "\\b";
                break;
            case '\f':
                result += "\\f";
                break;
            case '\n':
                result += "\\n";
                break;
            case '\r':
                result += "\\r";
                break;
            case '\t':
                result += "\\t";
                break;
            default:
                if (character < 0x20 or character == 0x7f) {
                    result += "\\u00";
                    result += hex[character >> 4];
                    result += hex[character & 0x0f];
                } else {
                    result += static_cast<char>(character);
                }
        }
    }
    result += '"';
    return result;
}
} // namespace

void SreprPrinter::bvisit(const Basic &x)
{
    std::ostringstream s;
    s << type_code_name(x.get_type_code());
    auto args = x.get_args();
    // Add and Mul store their canonical dictionaries in hash maps, so
    // get_args() iteration order is not a stable structural order. Sorting is
    // valid only for these commutative nodes; all other node argument order is
    // semantically significant and must be preserved.
    if (is_a<Add>(x) or is_a<Mul>(x)) {
        std::sort(args.begin(), args.end(), RCPBasicKeyLess());
    }
    if (is_a_sub<Constant>(x)) {
        SYMENGINE_ASSERT(args.size() == 0);
        s << '.' << down_cast<Constant const &>(x).get_name();
    } else {
        s << '(';
        if (args.size() == 0 && is_a_sub<Number>(x)) {
            s << x.__str__();
        } else {
            bool first{true};
            for (auto const &arg : args) {
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

void SreprPrinter::bvisit(const Symbol &x)
{
    str_ = "Symbol(" + quoted(x.get_name()) + ")";
}

void SreprPrinter::bvisit(const Dummy &x)
{
    str_ = "Dummy(" + quoted(x.get_name()) + ", "
           + std::to_string(x.get_index()) + ")";
}

void SreprPrinter::print_function(const char *type_name,
                                  const FunctionSymbol &x)
{
    std::ostringstream s;
    s << type_name << '(' << quoted(x.get_name());
    for (auto const &arg : x.get_vec()) {
        s << ", " << apply(arg);
    }
    s << ')';
    str_ = s.str();
}

void SreprPrinter::bvisit(const FunctionSymbol &x)
{
    print_function("FunctionSymbol", x);
}

void SreprPrinter::bvisit(const FunctionWrapper &x)
{
    print_function("FunctionWrapper", x);
}

std::string SreprPrinter::apply(const RCP<const Basic> &b)
{
    b->accept(*this);
    return str_;
}

} // namespace SymEngine
