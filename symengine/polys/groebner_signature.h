#ifndef SYMENGINE_POLYS_GROEBNER_SIGNATURE_H
#define SYMENGINE_POLYS_GROEBNER_SIGNATURE_H

#include <symengine/polys/groebner.h>

namespace SymEngine {

struct Signature {
    PackedMonomial monomial;
    unsigned generator_index = 0;   // 1-indexed: 1 .. input_polys

    bool operator==(const Signature& other) const {
        return generator_index == other.generator_index && monomial == other.monomial;
    }
};

struct SignatureLess {
    MonomialOrder order;
    explicit SignatureLess(MonomialOrder ord = MonomialOrder::DegRevLex) : order(ord) {}

    bool operator()(const Signature& a, const Signature& b) const {
        if (a.generator_index != b.generator_index) {
            return a.generator_index > b.generator_index;
        }
        return PackedMonomialLess{order}(a.monomial, b.monomial);
    }
};

inline Signature signature_mul(const Signature& s, const PackedMonomial& m) {
    Signature res;
    res.generator_index = s.generator_index;
    res.monomial = packed_multiply(s.monomial, m);
    return res;
}

} // namespace SymEngine

#endif // SYMENGINE_POLYS_GROEBNER_SIGNATURE_H
