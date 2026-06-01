#ifndef SYMENGINE_POLYS_GROEBNER_H
#define SYMENGINE_POLYS_GROEBNER_H

#include <symengine/basic.h>
#include <symengine/symbol.h>
#include <symengine/expression.h>
#include <symengine/sets.h>
#include <symengine/matrix.h>
#include <atomic>
#include <vector>

namespace SymEngine {

enum class MonomialOrder {
    Lex,
    GrLex,
    DegRevLex
};

enum class GroebnerAlgorithm {
    Buchberger,
    F5B,
    MoGVW,
    Auto
};

enum class GroebnerStatus {
    Success,
    Cancelled,
    ResourceLimitExceeded,
    NotZeroDimensional,
    UnsupportedCoefficientDomain
};

class GroebnerCancellationToken {
private:
    std::atomic<bool> cancelled_;

public:
    GroebnerCancellationToken() : cancelled_(false) {}
    GroebnerCancellationToken(const GroebnerCancellationToken &) = delete;
    GroebnerCancellationToken &operator=(const GroebnerCancellationToken &) = delete;

    void cancel() { cancelled_ = true; }
    bool is_cancelled() const { return cancelled_; }
};

struct GroebnerOptions {
    MonomialOrder order = MonomialOrder::DegRevLex;
    GroebnerAlgorithm algorithm = GroebnerAlgorithm::Auto;
    bool reduced = true;
    bool interreduce_input = true;
    bool sort_output = true;
    unsigned cancellation_check_interval = 1024;
    const GroebnerCancellationToken * cancellation_token = nullptr;
    unsigned max_s_pairs = 0;
    unsigned max_reduction_steps = 0;
    unsigned max_milliseconds = 0;
};

struct GroebnerStats {
    unsigned input_polys = 0;
    unsigned output_polys = 0;
    unsigned s_pairs_processed = 0;
    unsigned reductions_to_zero = 0;
    unsigned max_basis_size = 0;
};

struct GroebnerResult {
    vec_basic basis;
    vec_sym variables;
    MonomialOrder order;
    GroebnerStatus status;
    GroebnerStats stats;
};

struct TriangularUnivariateSolution {
    RCP<const Basic> univariate_polynomial;
    RCP<const Symbol> root_variable;
    vec_sym linear_variables;
    DenseMatrix linear_solution;
};

GroebnerResult groebner_basis(const vec_basic &polys,
                              const vec_sym &variables,
                              const GroebnerOptions &options = GroebnerOptions());

GroebnerResult fglm_convert(const GroebnerResult &source,
                            MonomialOrder target_order,
                            const GroebnerOptions &options = GroebnerOptions());

GroebnerResult augmented_groebner_basis(const vec_basic &polys,
                                        const vec_sym &variables,
                                        const RCP<const Symbol> &selected_variable,
                                        const RCP<const Symbol> &auxiliary_variable,
                                        const GroebnerOptions &options = GroebnerOptions());

TriangularUnivariateSolution
extract_univariate_linear_shape(const GroebnerResult &lex_basis,
                                const RCP<const Symbol> &root_variable);

RCP<const Set> solve_poly_system(const vec_basic &equations,
                                 const vec_sym &variables,
                                 const GroebnerOptions &options = GroebnerOptions());

RCP<const Set> solve_poly_system_via_univariate_root(const vec_basic &equations,
                                                     const vec_sym &variables,
                                                     const RCP<const Symbol> &selected_variable,
                                                     const RCP<const Symbol> &auxiliary_variable,
                                                     const GroebnerOptions &options = GroebnerOptions());

} // namespace SymEngine

#endif // SYMENGINE_POLYS_GROEBNER_H
