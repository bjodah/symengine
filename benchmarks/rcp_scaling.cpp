/*
 * rcp_scaling -- how well does expression building scale across threads?
 *
 * The workload is embarrassingly parallel by construction: every thread owns
 * its own symbols and its own expressions and never publishes them, so two
 * threads share nothing except the library's process-global singletons
 * (`zero`, `one`, `minus_one`, `two`, `pi`, ... -- constants.cpp).  Those are
 * ordinary heap `Basic` objects with an intrusive reference count, and every
 * `Mul` that is built stores a reference to `one` in its coefficient slot,
 * every `Add` a reference to `zero`, every `sub()` a reference to
 * `minus_one`.  So each thread writes the same handful of cache lines
 * millions of times a second and the hardware serialises them.
 *
 * If scaling is linear, the singleton reference counts are not a bottleneck.
 * If it flattens out while `nproc` is nowhere near exhausted, they are.
 *
 * Usage:  rcp_scaling [threads...]     (default: 1 2 4 8 16)
 *         SYMENGINE_RCP_SCALING_ITERS=<n>   per-thread iterations
 *         SYMENGINE_RCP_SCALING_REPEAT=<n>  timed repeats, best is reported
 */

#include <atomic>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

#include <symengine/add.h>
#include <symengine/basic.h>
#include <symengine/constants.h>
#include <symengine/integer.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/symbol.h>

using SymEngine::add;
using SymEngine::Basic;
using SymEngine::div;
using SymEngine::hash_t;
using SymEngine::integer;
using SymEngine::mul;
using SymEngine::pow;
using SymEngine::RCP;
using SymEngine::sub;
using SymEngine::symbol;

// Enough work that thread start-up is noise, small enough that the whole
// sweep runs in a couple of minutes.
static long g_iters = 60000;
static int g_repeat = 3;

//! One thread's share.  Everything it touches is private except the
//! singletons that `add`/`mul`/`sub`/`div`/`pow` reach for internally.
static hash_t worker(int id)
{
    // Per-thread symbols: distinct names, so no two threads ever hold a
    // reference to the same Symbol.
    RCP<const Basic> x = symbol("x_" + std::to_string(id));
    RCP<const Basic> y = symbol("y_" + std::to_string(id));
    RCP<const Basic> z = symbol("z_" + std::to_string(id));

    hash_t acc = 0;
    for (long i = 0; i < g_iters; i++) {
        // Mul with an implicit coefficient of `one`.
        RCP<const Basic> a = mul(x, y);
        // Add with an implicit coefficient of `zero`.
        RCP<const Basic> b = add(a, z);
        // sub() multiplies by `minus_one`.
        RCP<const Basic> c = sub(b, x);
        // Small-integer exponent; Pow canonicalisation compares against
        // `zero`/`one`.
        RCP<const Basic> d = pow(c, integer(2 + (i & 1)));
        // div() is mul-by-a-negative-power.
        RCP<const Basic> e = div(d, y);
        acc ^= e->hash();
        // a..e die here: five expression trees torn down, every intrusive
        // count in them decremented, including the singletons'.
    }
    return acc;
}

static double run(int nthreads)
{
    std::vector<std::thread> threads;
    std::vector<hash_t> results(static_cast<size_t>(nthreads), 0);
    std::atomic<int> ready{0};
    std::atomic<bool> go{false};

    threads.reserve(static_cast<size_t>(nthreads));
    for (int t = 0; t < nthreads; t++) {
        threads.emplace_back([t, &results, &ready, &go]() {
            ready.fetch_add(1, std::memory_order_relaxed);
            while (!go.load(std::memory_order_acquire)) {
                std::this_thread::yield();
            }
            results[static_cast<size_t>(t)] = worker(t);
        });
    }
    while (ready.load(std::memory_order_relaxed) < nthreads) {
        std::this_thread::yield();
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    go.store(true, std::memory_order_release);
    for (auto &th : threads) {
        th.join();
    }
    auto t2 = std::chrono::high_resolution_clock::now();

    hash_t sink = 0;
    for (auto r : results) {
        sink ^= r;
    }
    if (sink == 1) { // never true in practice; keeps `results` alive
        std::cout << "";
    }
    return std::chrono::duration<double>(t2 - t1).count();
}

static long env_long(const char *name, long fallback)
{
    if (const char *s = std::getenv(name)) {
        char *end = nullptr;
        long v = std::strtol(s, &end, 10);
        if (end != s && v > 0) {
            return v;
        }
    }
    return fallback;
}

int main(int argc, char *argv[])
{
    g_iters = env_long("SYMENGINE_RCP_SCALING_ITERS", g_iters);
    g_repeat = static_cast<int>(env_long("SYMENGINE_RCP_SCALING_REPEAT", 3));

    std::vector<int> sweep;
    for (int i = 1; i < argc; i++) {
        int n = std::atoi(argv[i]);
        if (n > 0) {
            sweep.push_back(n);
        }
    }
    if (sweep.empty()) {
        sweep = {1, 2, 4, 8, 16};
    }

    std::cout << "rcp_scaling: " << g_iters << " iterations/thread, best of "
              << g_repeat << ", hardware_concurrency="
              << std::thread::hardware_concurrency() << "\n";
    std::cout << std::setw(8) << "threads" << std::setw(12) << "seconds"
              << std::setw(14) << "Mops/s" << std::setw(10) << "speedup"
              << std::setw(12) << "efficiency"
              << "\n";

    double base = 0.0;
    for (int n : sweep) {
        double best = 0.0;
        for (int r = 0; r < g_repeat; r++) {
            double t = run(n);
            if (r == 0 || t < best) {
                best = t;
            }
        }
        double total_ops = static_cast<double>(g_iters) * n;
        if (base == 0.0) {
            base = best / (static_cast<double>(g_iters) * sweep.front());
        }
        double per_op_1 = base;
        double speedup = per_op_1 * total_ops / best;
        std::cout << std::setw(8) << n << std::setw(12) << std::fixed
                  << std::setprecision(3) << best << std::setw(14)
                  << std::setprecision(3) << (total_ops / best / 1e6)
                  << std::setw(10) << std::setprecision(2) << speedup << "x"
                  << std::setw(11) << std::setprecision(0)
                  << (100.0 * speedup / n) << "%"
                  << "\n";
    }
    return 0;
}
