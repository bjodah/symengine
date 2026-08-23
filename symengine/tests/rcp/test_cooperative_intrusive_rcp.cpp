#include "catch.hpp"

#include <atomic>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <thread>

#include <symengine/add.h>
#include <symengine/constants.h>
#include <symengine/dict.h>
#include <symengine/integer.h>
#include <symengine/mul.h>
#include <symengine/symengine_rcp.h>
#include <symengine/symbol.h>

// Guard so this file is empty (no test cases) for non-cooperative_intrusive
// backends
#if defined(WITH_SYMENGINE_COOPERATIVE_INTRUSIVE_RCP)

// The cooperative counter lives in symengine_rcp_cooperative.cpp.
// We only need to install fake external-runtime hooks here.

using SymEngine::EnableRCPFromThis;
using SymEngine::make_rcp;
using SymEngine::RCP;

class Node : public EnableRCPFromThis<Node>
{
public:
    int v = 0;
};

#if defined(WITH_SYMENGINE_THREAD_SAFE)
static std::atomic<unsigned> threaded_node_destructions(0);

class ThreadedNode : public EnableRCPFromThis<ThreadedNode>
{
public:
    ~ThreadedNode()
    {
        ++threaded_node_destructions;
    }
};
#endif

// Counting hooks for testing external-owned mode without a real runtime
static int g_inc_count = 0;
static int g_dec_count = 0;

static void counting_inc_hook(void *) noexcept
{
    ++g_inc_count;
}
static void counting_dec_hook(void *) noexcept
{
    ++g_dec_count;
}

// Fake foreign-runtime pointer. The two low bits of the state word are the
// union's tags -- bit 0 count-or-pointer, bit 1 immortal -- so a wrapper
// pointer must have both clear, and set_self_external() aborts otherwise.
alignas(4) static char fake_runtime_storage;
static void *fake_runtime = &fake_runtime_storage;

// Install hooks once at process start
struct HookInstaller {
    HookInstaller()
    {
        SymEngine::cooperative_intrusive_init(counting_inc_hook,
                                              counting_dec_hook);
    }
};
static HookInstaller installer;

TEST_CASE("C++-owned counting", "[cooperative_intrusive_rcp]")
{
    g_inc_count = 0;
    g_dec_count = 0;

    // Fresh object before any RCP
    Node *raw = new Node();
    REQUIRE(raw->use_count() == 0);
    REQUIRE(not raw->is_external_owned());
    delete raw;

    // With RCP
    RCP<Node> a = make_rcp<Node>();
    REQUIRE(a->use_count() == 1);
    REQUIRE(a->is_uniquely_owned_by_cpp());
    REQUIRE(a->is_uniquely_owned());
    REQUIRE(not a->is_external_owned());

    {
        RCP<Node> b = a;
        REQUIRE(a->use_count() == 2);
        REQUIRE(not a->is_uniquely_owned_by_cpp());
        REQUIRE(not a->is_uniquely_owned());
        REQUIRE(not a->is_external_owned());

        {
            RCP<Node> c = b;
            REQUIRE(a->use_count() == 3);
        }
        REQUIRE(a->use_count() == 2);
    }
    REQUIRE(a->use_count() == 1);
    REQUIRE(a->is_uniquely_owned_by_cpp());
    REQUIRE(a->is_uniquely_owned());
}

#if defined(WITH_SYMENGINE_THREAD_SAFE)
TEST_CASE("C++-owned last release crosses threads",
          "[cooperative_intrusive_rcp][threading]")
{
    threaded_node_destructions = 0;
    RCP<ThreadedNode> main_ref = make_rcp<ThreadedNode>();
    RCP<ThreadedNode> worker_ref = main_ref;

    std::mutex mutex;
    std::condition_variable ready;
    bool worker_ready = false;
    bool release = false;
    std::thread worker([ref = std::move(worker_ref), &mutex, &ready,
                        &worker_ready, &release]() mutable {
        {
            std::unique_lock<std::mutex> lock(mutex);
            worker_ready = true;
            ready.notify_one();
            ready.wait(lock, [&release] { return release; });
        }
        ref.reset();
    });

    {
        std::unique_lock<std::mutex> lock(mutex);
        ready.wait(lock, [&worker_ready] { return worker_ready; });
        release = true;
        ready.notify_one();
    }
    // The unlock above publishes the start signal. Neither reset is ordered
    // after the other, so one thread releases and the other deletes.
    main_ref.reset();
    worker.join();

    REQUIRE(threaded_node_destructions == 1);
}
#endif

TEST_CASE("Handoff changes mode", "[cooperative_intrusive_rcp]")
{
    g_inc_count = 0;
    g_dec_count = 0;

    RCP<Node> a = make_rcp<Node>();
    REQUIRE(a->use_count() == 1);
    REQUIRE(not a->is_external_owned());

    // Hand off to external-owned mode
    a->set_self_external(fake_runtime);

    REQUIRE(a->is_external_owned());
    REQUIRE(a->self_external() == fake_runtime);
    REQUIRE(a->use_count() == 0);
    REQUIRE(not a->is_uniquely_owned_by_cpp());
    REQUIRE(not a->is_uniquely_owned());

    // The inc hook should have been called once per pre-existing reference
    // (replayed count). With 1 C++ reference, that's 1 call.
    REQUIRE(g_inc_count == 1);

    // The fake dec-hook doesn't free the object (no real Py_DECREF), so the
    // Node leaks even after releasing the RCP handle below.  ASAN leak
    // detection is disabled for this test via CTest ENVIRONMENT.
    a = SymEngine::null;
}

TEST_CASE("External-owned inc/dec route through hooks",
          "[cooperative_intrusive_rcp]")
{
    g_inc_count = 0;
    g_dec_count = 0;

    RCP<Node> a = make_rcp<Node>();
    a->set_self_external(fake_runtime);
    REQUIRE(a->is_external_owned());

    // After handoff, hook was called once (replay)
    int inc_after_handoff = g_inc_count;
    int dec_after_handoff = g_dec_count;

    // Copying an RCP should call the inc hook
    {
        RCP<Node> b = a;
        REQUIRE(g_inc_count == inc_after_handoff + 1);
        REQUIRE(a->use_count() == 0); // Always 0 when external-owned
        REQUIRE(a->is_external_owned());

        // Destroying b should call the dec hook
    }
    REQUIRE(g_dec_count == dec_after_handoff + 1);
    REQUIRE(a->use_count() == 0);
    REQUIRE(a->is_external_owned());

    // The fake dec-hook doesn't free the object (no real Py_DECREF), so the
    // Node leaks even after releasing the RCP handle below.  ASAN leak
    // detection is disabled for this test via CTest ENVIRONMENT.
    a = SymEngine::null;
}

TEST_CASE("use_count() is never 1 when external-owned",
          "[cooperative_intrusive_rcp]")
{
    g_inc_count = 0;
    g_dec_count = 0;

    RCP<Node> a = make_rcp<Node>();
    a->set_self_external(fake_runtime);

    // Across several hook-driven inc/decs, the uniqueness predicate must
    // always be false — the safety contract for dictionary stealing.
    REQUIRE(a->use_count() == 0);
    REQUIRE(not a->is_uniquely_owned_by_cpp());
    REQUIRE(not a->is_uniquely_owned());

    {
        RCP<Node> b = a;
        REQUIRE(a->use_count() == 0);
        REQUIRE(not a->is_uniquely_owned_by_cpp());
        REQUIRE(not a->is_uniquely_owned());

        {
            RCP<Node> c = b;
            REQUIRE(a->use_count() == 0);
            REQUIRE(not a->is_uniquely_owned_by_cpp());
            REQUIRE(not a->is_uniquely_owned());
        }
        REQUIRE(a->use_count() == 0);
        REQUIRE(not a->is_uniquely_owned_by_cpp());
        REQUIRE(not a->is_uniquely_owned());
    }
    REQUIRE(a->use_count() == 0);
    REQUIRE(not a->is_uniquely_owned_by_cpp());
    REQUIRE(not a->is_uniquely_owned());

    // Cleanup: drop our RCP reference. dec_ref() invokes the fake dec-hook
    // (incrementing g_dec_count) and returns false, so the RCP destructor
    // does not delete the object.  The Node leaks — acceptable for a test
    // that validates the counter, not real Python lifetime management.
    a = SymEngine::null;
}

TEST_CASE("detach external and replay inline references",
          "[cooperative_intrusive_rcp]")
{
    g_inc_count = 0;
    g_dec_count = 0;

    RCP<Node> a = make_rcp<Node>();
    RCP<Node> b = a;
    REQUIRE(a->use_count() == 2);
    REQUIRE(not a->is_external_owned());

    a->set_self_external(fake_runtime);
    REQUIRE(a->is_external_owned());
    REQUIRE(a->self_external() == fake_runtime);
    REQUIRE(a->use_count() == 0);
    REQUIRE(g_inc_count == 2);

    void *detached = a->detach_external();
    REQUIRE(detached == fake_runtime);
    REQUIRE(not a->is_external_owned());
    REQUIRE(a->self_external() == nullptr);
    REQUIRE(a->use_count() == 0);

    a->inc_ref();
    a->inc_ref();
    REQUIRE(a->use_count() == 2);

    b = SymEngine::null;
    REQUIRE(a->use_count() == 1);
}

TEST_CASE("Add::from_dict retains an external Mul dictionary",
          "[cooperative_intrusive_rcp]")
{
    g_inc_count = 0;
    g_dec_count = 0;

    RCP<const SymEngine::Basic> x = SymEngine::symbol("x");
    RCP<const SymEngine::Basic> y = SymEngine::symbol("y");
    RCP<const SymEngine::Basic> term = SymEngine::mul(x, y);
    RCP<const SymEngine::Mul> term_mul
        = SymEngine::rcp_static_cast<const SymEngine::Mul>(term);
    const auto original_size = term_mul->get_dict().size();

    term->set_self_external(fake_runtime);
    REQUIRE(not term->is_uniquely_owned());

    SymEngine::umap_basic_num terms;
    terms.emplace(term, SymEngine::integer(2));
    RCP<const SymEngine::Basic> collapsed
        = SymEngine::Add::from_dict(SymEngine::zero, std::move(terms));

    // An external wrapper can still hand this object back to a foreign
    // runtime. The copy path must therefore leave the original Mul intact.
    REQUIRE(term_mul->get_dict().size() == original_size);
    REQUIRE(SymEngine::eq(*collapsed, *SymEngine::mul(SymEngine::integer(2),
                                                      SymEngine::mul(x, y))));

    // The fake runtime hooks deliberately do not free the external object.
    term = SymEngine::null;
}

TEST_CASE("Add::from_dict collapses a uniquely C++-owned Mul",
          "[cooperative_intrusive_rcp]")
{
    RCP<const SymEngine::Basic> x = SymEngine::symbol("x");
    RCP<const SymEngine::Basic> y = SymEngine::symbol("y");
    RCP<const SymEngine::Basic> term = SymEngine::mul(x, y);
    REQUIRE(term->is_uniquely_owned());

    SymEngine::umap_basic_num terms;
    terms.emplace(std::move(term), SymEngine::integer(2));
    REQUIRE(term.is_null());

    RCP<const SymEngine::Basic> collapsed
        = SymEngine::Add::from_dict(SymEngine::zero, std::move(terms));
    REQUIRE(SymEngine::eq(*collapsed, *SymEngine::mul(SymEngine::integer(2),
                                                      SymEngine::mul(x, y))));
}

// ---------------------------------------------------------------------------
// Immortality and the foreign hand-off.
//
// An immortal object keeps no reference count, so there is nothing to fold
// into a wrapper's count when a foreign runtime attaches one -- and nothing
// the runtime would ever get back, since C++ never releases the object. It is
// therefore given exactly one permanent reference and stays immortal, which is
// what keeps self_external() from returning a dangling pointer once the
// runtime's own users are done with the wrapper. Everything a binding can
// observe behaves as it does for an ordinary externalised object.
// ---------------------------------------------------------------------------

// A second fake runtime object, so the immortal cases below do not collide
// with the mortal ones above on the single external-owner slot.
alignas(4) static char fake_runtime4_storage;
static void *fake_runtime4 = &fake_runtime4_storage;

// Parked in statics for the same reason constants.cpp keeps its singletons in
// static storage: an immortal object outlives every handle to it, and what
// makes that a retention rather than a leak is that something still names it.
static Node *immortal_node = nullptr;
static Node *immortal_external_node = nullptr;

TEST_CASE("Immortal objects keep no count", "[cooperative_intrusive_rcp]")
{
    RCP<Node> a = make_rcp<Node>();
    a->v = 7;
    immortal_node = a.get();
    REQUIRE(a->use_count() == 1);

    a->mark_immortal();
    REQUIRE(a->is_immortal());
    REQUIRE(not a->is_external_owned());
    // "Effectively infinite": no count is kept, so no honest finite answer
    // exists, and every ownership gate must read it as "not mine to steal".
    REQUIRE(a->use_count() == ~0u);
    REQUIRE(not a->is_uniquely_owned_by_cpp());
    REQUIRE(not a->is_uniquely_owned());

    {
        RCP<Node> b = a;
        RCP<Node> c = b;
        REQUIRE(a->use_count() == ~0u);
    }
    a.reset();
    REQUIRE(immortal_node->v == 7); // survived its last handle, as promised
    REQUIRE(immortal_node->is_immortal());
}

TEST_CASE("Immortal hand-off grants one permanent reference",
          "[cooperative_intrusive_rcp]")
{
    g_inc_count = 0;
    g_dec_count = 0;

    RCP<Node> a = make_rcp<Node>();
    {
        RCP<Node> b = a; // a second C++ reference, deliberately
        a->mark_immortal();
    }

    a->set_self_external(fake_runtime4);

    // One incref, not one-per-outstanding-reference: there is no outstanding
    // reference count to replay, and C++ will never hand any back.
    REQUIRE(g_inc_count == 1);
    REQUIRE(g_dec_count == 0);

    REQUIRE(a->is_external_owned());
    REQUIRE(a->self_external() == fake_runtime4);
    REQUIRE(a->is_immortal());
    REQUIRE(a->use_count() == 0);
    REQUIRE(not a->is_uniquely_owned_by_cpp());
    REQUIRE(not a->is_uniquely_owned());

    // C++ reference traffic stays off the hooks *and* off the state word --
    // which is the point: an immortal object keeps its scaling even under a
    // language binding.
    {
        RCP<Node> b = a;
        RCP<Node> c = b;
    }
    REQUIRE(g_inc_count == 1);
    REQUIRE(g_dec_count == 0);

    // Shutdown: the binding takes its wrapper back. The object returns to
    // immortal-without-a-wrapper, not to a count of zero, so the RCP released
    // below is not an underflow.
    immortal_external_node = a.get();
    REQUIRE(a->detach_external() == fake_runtime4);
    REQUIRE(not a->is_external_owned());
    REQUIRE(a->is_immortal());
    REQUIRE(a->self_external() == nullptr);
    REQUIRE(a->detach_external() == nullptr);
    a.reset();
    REQUIRE(immortal_external_node->is_immortal());
    REQUIRE(g_dec_count == 0);
}

#if defined(WITH_SYMENGINE_THREAD_SAFE)
// ---------------------------------------------------------------------------
// mark_immortal() landing *inside* set_self_external().
//
// Externalization folds the object's C++ count into the wrapper's own, one
// foreign incref per outstanding reference, and publishes the wrapper with a
// compare-exchange. Marking recognised on the initial load takes a different
// protocol -- one permanent reference, because an immortal object has no count
// to fold -- but a mark that lands *after* that load and before the exchange
// used to be read back as an ordinary count. The band's counts are
// astronomically large, so "replay this count" meant ~2^62 foreign increfs:
// not a wrong answer but a hang, which is worse, because a hung test reports
// nothing at all.
//
// The regression is deterministic, not a stress test. The foreign incref hook
// blocks on its first call, which puts externalization exactly in the window:
// the count has been read and the wrapper is not yet published. The mark is
// then made from this thread, and only afterwards is the hook released.
// Promptness is asserted by the hook itself: past a small bound it aborts,
// so the pre-fix behaviour is a loud failure in milliseconds rather than a
// test run that never ends.
// ---------------------------------------------------------------------------

alignas(4) static char fake_runtime5_storage;
static void *fake_runtime5 = &fake_runtime5_storage;
static Node *immortal_raced_node = nullptr;

static std::mutex g_gate_mutex;
static std::condition_variable g_gate_cv;
static bool g_gate_entered = false;
static bool g_gate_open = false;
static std::atomic<long> g_gated_incs(0);
static std::atomic<long> g_gated_decs(0);

// Replaying a band count would be ~2^62 calls. Anything past a handful is
// already the bug, and stopping here is what turns a hang into a verdict.
static const long gated_incref_bound = 4096;

static void gated_inc_hook(void *) noexcept
{
    const long n = ++g_gated_incs;
    if (n > gated_incref_bound) {
        std::fprintf(stderr,
                     "set_self_external replayed %ld foreign increfs: it read "
                     "the immortal band as an ordinary count\n",
                     n);
        std::abort();
    }
    if (n == 1) {
        std::unique_lock<std::mutex> lock(g_gate_mutex);
        g_gate_entered = true;
        g_gate_cv.notify_all();
        while (!g_gate_open) {
            g_gate_cv.wait(lock);
        }
    }
}

static void gated_dec_hook(void *) noexcept
{
    ++g_gated_decs;
}

TEST_CASE("mark_immortal racing set_self_external",
          "[cooperative_intrusive_rcp]")
{
    SymEngine::cooperative_intrusive_init(gated_inc_hook, gated_dec_hook);

    RCP<Node> a = make_rcp<Node>();
    a->v = 11;
    RCP<Node> b = a;
    RCP<Node> c = b;
    REQUIRE(a->use_count() == 3);
    immortal_raced_node = a.get();

    Node *raw = a.get();
    std::thread externalizer(
        [raw]() { raw->set_self_external(fake_runtime5); });

    {
        std::unique_lock<std::mutex> lock(g_gate_mutex);
        while (!g_gate_entered) {
            g_gate_cv.wait(lock);
        }
    }
    // Externalization is now parked inside its first foreign incref: it has
    // read the count and has not published the wrapper. This is the window.
    a->mark_immortal();
    {
        std::lock_guard<std::mutex> lock(g_gate_mutex);
        g_gate_open = true;
    }
    g_gate_cv.notify_all();
    externalizer.join();

    // Bounded work: the outstanding references may be replayed (the loser did
    // not know yet), but nothing beyond them.
    REQUIRE(g_gated_incs.load() <= 8);
    // ... and the wrapper is left holding exactly the one permanent reference
    // the immortal hand-off promises, not three and not zero.
    REQUIRE(g_gated_incs.load() - g_gated_decs.load() == 1);

    REQUIRE(a->is_external_owned());
    REQUIRE(a->self_external() == fake_runtime5);
    REQUIRE(a->is_immortal());
    REQUIRE(a->use_count() == 0);

    // Later C++ traffic stays off the hooks, as for any immortal object.
    const long incs = g_gated_incs.load();
    const long decs = g_gated_decs.load();
    {
        RCP<Node> d = a;
        RCP<Node> e = d;
    }
    REQUIRE(g_gated_incs.load() == incs);
    REQUIRE(g_gated_decs.load() == decs);

    a.reset();
    b.reset();
    c.reset();
    REQUIRE(immortal_raced_node->v == 11);
    REQUIRE(immortal_raced_node->is_immortal());

    SymEngine::cooperative_intrusive_init(counting_inc_hook, counting_dec_hook);
}
#endif // WITH_SYMENGINE_THREAD_SAFE

#endif // WITH_SYMENGINE_COOPERATIVE_INTRUSIVE_RCP
