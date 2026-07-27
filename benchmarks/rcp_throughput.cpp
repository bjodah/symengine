#include <iostream>
#include <chrono>
#include <iomanip>
#include <vector>
#include <string>

#include <symengine/basic.h>
#include <symengine/add.h>
#include <symengine/symbol.h>
#include <symengine/integer.h>
#include <symengine/mul.h>
#include <symengine/pow.h>

using SymEngine::add;
using SymEngine::Basic;
using SymEngine::integer;
using SymEngine::Integer;
using SymEngine::mul;
using SymEngine::Number;
using SymEngine::pow;
using SymEngine::RCP;
using SymEngine::rcp_static_cast;
using SymEngine::Symbol;
using SymEngine::symbol;
using SymEngine::vec_basic;

static const int N_SYMBOLS = 50000;
static const int N_INTEGERS = 50000;
static const int N_ADD_CHAIN = 5000;
static const int N_MUL_EXPRS = 1000;
static const int N_DROP_CYCLES = 10000;

double bench_symbol_creation()
{
    auto t1 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N_SYMBOLS; i++) {
        RCP<const Basic> s = symbol("x" + std::to_string(i));
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double>(t2 - t1).count();
}

double bench_integer_creation()
{
    auto t1 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N_INTEGERS; i++) {
        RCP<const Basic> n = integer(i);
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double>(t2 - t1).count();
}

double bench_add_chain()
{
    RCP<const Basic> x = symbol("x");
    RCP<const Basic> a = x;

    auto t1 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N_ADD_CHAIN; i++) {
        a = add(a, integer(i));
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double>(t2 - t1).count();
}

double bench_mul_build()
{
    RCP<const Basic> x = symbol("x");
    RCP<const Basic> y = symbol("y");
    std::vector<RCP<const Basic>> exprs;
    exprs.reserve(N_MUL_EXPRS);

    auto t1 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N_MUL_EXPRS; i++) {
        exprs.push_back(mul(integer(i + 1), pow(x, integer(i % 20 + 1))));
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double>(t2 - t1).count();
}

double bench_expr_build_add()
{
    RCP<const Basic> x = symbol("x");
    RCP<const Basic> y = symbol("y");
    RCP<const Basic> z = symbol("z");
    std::vector<RCP<const Basic>> exprs;
    exprs.reserve(N_MUL_EXPRS);

    auto t1 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N_MUL_EXPRS; i++) {
        RCP<const Basic> e
            = add(add(pow(x, integer(i % 10 + 2)),
                      pow(y, integer(i % 8 + 1))),
                  pow(z, integer(i % 6 + 1)));
        exprs.push_back(e);
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double>(t2 - t1).count();
}

double bench_drop_expressions()
{
    auto t1 = std::chrono::high_resolution_clock::now();
    for (int cycle = 0; cycle < N_DROP_CYCLES; cycle++) {
        std::vector<RCP<const Basic>> batch;
        batch.reserve(10);
        for (int j = 0; j < 10; j++) {
            RCP<const Basic> s = symbol("s" + std::to_string(j));
            RCP<const Basic> e = add(pow(s, integer(3)),
                                     mul(integer(j + 1), s));
            batch.push_back(e);
        }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double>(t2 - t1).count();
}

int main(int argc, char *argv[])
{
    SymEngine::print_stack_on_segfault();

    std::cout << "RCP Throughput Benchmark" << std::endl;
    std::cout << "========================" << std::endl;

    struct Result {
        const char *name;
        int count;
        double time;
    };

    Result results[] = {
        {"symbol_creation", N_SYMBOLS, bench_symbol_creation()},
        {"integer_creation", N_INTEGERS, bench_integer_creation()},
        {"add_chain", N_ADD_CHAIN, bench_add_chain()},
        {"mul_build", N_MUL_EXPRS, bench_mul_build()},
        {"expr_build_add", N_MUL_EXPRS, bench_expr_build_add()},
        {"drop_expressions", N_DROP_CYCLES, bench_drop_expressions()},
    };

    std::cout << std::endl;
    std::cout << std::left << std::setw(20) << "benchmark"
              << std::right << std::setw(12) << "count"
              << std::setw(15) << "time(s)"
              << std::setw(18) << "ops/sec" << std::endl;
    std::cout << std::string(65, '-') << std::endl;

    for (const auto &r : results) {
        double ops = (r.time > 0.0) ? (r.count / r.time) : 0.0;
        std::cout << std::left << std::setw(20) << r.name
                  << std::right << std::setw(12) << r.count
                  << std::setw(15) << std::fixed << std::setprecision(6)
                  << r.time
                  << std::setw(18) << std::fixed << std::setprecision(0)
                  << ops << std::endl;
    }

    std::cout << std::endl;
    std::cout << "parse_begin" << std::endl;
    for (const auto &r : results) {
        double ops = (r.time > 0.0) ? (r.count / r.time) : 0.0;
        std::cout << r.name << "," << r.count << ","
                  << std::fixed << std::setprecision(6) << r.time << ","
                  << std::fixed << std::setprecision(0) << ops << std::endl;
    }
    std::cout << "parse_end" << std::endl;

    return 0;
}
