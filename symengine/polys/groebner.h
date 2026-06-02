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
    M4GB,
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

    // Finite-field coefficient domain modulus. 0 -> rational/Expression as today.
    // Any prime p with 2 <= p < 2^31 is valid.
    uint64_t modulus = 0;

    // Matrix / degree controls.
    unsigned matrix_batch_size = 0;     // 0 = algorithm default
    unsigned max_matrix_rows = 0;       // 0 = unlimited
    unsigned max_matrix_columns = 0;    // 0 = unlimited
    unsigned max_degree = 0;            // 0 = unlimited

    // Verification toggle for tests / debugging.
    bool verify_with_buchberger = false;
};

struct GroebnerStats {
    unsigned input_polys = 0;
    unsigned output_polys = 0;
    unsigned s_pairs_processed = 0;
    unsigned reductions_to_zero = 0;
    unsigned max_basis_size = 0;

    // F5B + MoGVW signature-based counters
    unsigned rejected_by_syzygy = 0;
    unsigned rejected_by_rewritten = 0;
    unsigned f5_reductions = 0;
    unsigned labeled_monomials_lifted = 0;
    unsigned rejected_by_lcm = 0;
    unsigned collisions_resolved = 0;

    // M4GB / F4 matrix counters
    unsigned matrices_built = 0;
    unsigned max_matrix_rows_seen = 0;
    unsigned max_matrix_columns_seen = 0;
    unsigned rows_reduced_to_zero = 0;
    unsigned polymatrix_entries_built = 0;
    unsigned polymatrix_generation = 0;

    // Parametric assumptions
    std::vector<RCP<const Basic>> genericity_assumptions;
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

// Compute the normal form (remainder) of poly modulo the basis G under the
// supplied monomial order. variables fixes the indexing of exponent vectors,
// the same convention used by groebner_basis.
RCP<const Basic> normal_form(const RCP<const Basic> &poly,
                             const vec_basic &G,
                             const vec_sym &variables,
                             MonomialOrder order = MonomialOrder::DegRevLex);

// Returns true if G is a Groebner basis under the given order: every
// S-polynomial of pairs of G reduces to zero modulo G.
bool is_groebner(const vec_basic &G,
                 const vec_sym &variables,
                 MonomialOrder order = MonomialOrder::DegRevLex);

// Returns true if G is a reduced Groebner basis: it is a Groebner basis, every
// element is monic, and no term of any element is divisible by the leading
// monomial of any other element.
bool is_reduced_basis(const vec_basic &G,
                      const vec_sym &variables,
                      MonomialOrder order = MonomialOrder::DegRevLex);

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
