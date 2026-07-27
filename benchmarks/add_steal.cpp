#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <symengine/add.h>
#include <symengine/constants.h>
#include <symengine/dict.h>
#include <symengine/integer.h>
#include <symengine/mul.h>
#include <symengine/symbol.h>

using SymEngine::Add;
using SymEngine::Basic;
using SymEngine::integer;
using SymEngine::map_basic_basic;
using SymEngine::Mul;
using SymEngine::one;
using SymEngine::RCP;
using SymEngine::symbol;
using SymEngine::umap_basic_num;
using SymEngine::zero;

int main(int argc, char *argv[])
{
    const std::size_t iterations
        = argc > 1 ? std::strtoull(argv[1], nullptr, 10) : 200000;
    const std::size_t factors
        = argc > 2 ? std::strtoull(argv[2], nullptr, 10) : 16;

    if (iterations == 0 || factors < 2) {
        std::cerr << "iterations must be positive and factors must be at least 2\n";
        return 1;
    }

    std::vector<RCP<const Basic>> symbols;
    symbols.reserve(factors);
    for (std::size_t i = 0; i < factors; ++i)
        symbols.push_back(symbol("x" + std::to_string(i)));

    std::uint64_t checksum = 0;
    const auto start = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < iterations; ++i) {
        map_basic_basic mul_dict;
        for (const auto &s : symbols)
            SymEngine::insert(mul_dict, s, one);

        // Move the only RCP into the Add dictionary. This is the exact
        // one-entry Add::from_dict shape that can move Mul::dict_ instead of
        // allocating and copying a second map.
        RCP<const Basic> term = Mul::from_dict(one, std::move(mul_dict));
        umap_basic_num add_dict;
        add_dict.emplace(std::move(term), integer(2));
        RCP<const Basic> result = Add::from_dict(zero, std::move(add_dict));
        checksum += result->hash() + i;
    }
    const auto stop = std::chrono::steady_clock::now();

    const double elapsed_ns
        = std::chrono::duration<double, std::nano>(stop - start).count();
    std::cout << "iterations: " << iterations << '\n';
    std::cout << "factors: " << factors << '\n';
    std::cout << "elapsed_ns: " << elapsed_ns << '\n';
    std::cout << "ns_per_iteration: " << elapsed_ns / iterations << '\n';
    std::cout << "checksum: " << checksum << '\n';
}
