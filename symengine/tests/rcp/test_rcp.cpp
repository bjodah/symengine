#include "catch.hpp"

#include <array>
#include <atomic>
#include <thread>
#include <vector>

#include <symengine/symengine_rcp.h>
#include <symengine/constants.h>
#include <symengine/integer.h>

using SymEngine::EnableRCPFromThis;
using SymEngine::make_rcp;
using SymEngine::null;
using SymEngine::Ptr;
using SymEngine::RCP;

// This is the canonical use of EnableRCPFromThis:

class Mesh : public EnableRCPFromThis<Mesh>
{
public:
    int x, y;
};

TEST_CASE("Test make_rcp", "[rcp]")
{
    RCP<Mesh> m = make_rcp<Mesh>();
    Ptr<Mesh> p = m.ptr();
    REQUIRE(not(m == null));
    REQUIRE(p->use_count() == 1);
    RCP<Mesh> m2 = m;
    REQUIRE(p->use_count() == 2);
    RCP<Mesh> m3 = m2;
    REQUIRE(p->use_count() == 3);
}

void f(Mesh &m)
{
    REQUIRE(m.use_count() == 1);
    // rcp_from_this() gives up non const version of RCP<Mesh> because 'm' is
    // not const
    RCP<Mesh> m2 = m.rcp_from_this();
    REQUIRE(m.use_count() == 2);
    m2->x = 6;
}

void f_const(const Mesh &m)
{
    REQUIRE(m.use_count() == 1);
    // rcp_from_this() gives up const version of RCP<Mesh> because 'm' is const
    RCP<const Mesh> m2 = m.rcp_from_this();
    REQUIRE(m.use_count() == 2);
}

TEST_CASE("Test rcp_from_this", "[rcp]")
{
    RCP<Mesh> m = make_rcp<Mesh>();
    REQUIRE(m->use_count() == 1);
    m->x = 5;
    REQUIRE(m->x == 5);
    f(*m);
    REQUIRE(m->use_count() == 1);
    REQUIRE(m->x == 6);

    f_const(*m);
    REQUIRE(m->use_count() == 1);
}

TEST_CASE("Test rcp_from_this const", "[rcp]")
{
    RCP<const Mesh> m = make_rcp<const Mesh>();
    REQUIRE(m->use_count() == 1);
    f_const(*m);
    REQUIRE(m->use_count() == 1);
}

// This is not a canonical way how to use EnableRCPFromThis, since we use
// 'const Mesh2' for the internal weak pointer, so we can only get
// 'RCP<const Mesh2>' out of rcp_from_this(). But it is legitimate code, so we
// test it as well.

class Mesh2 : public EnableRCPFromThis<const Mesh2>
{
public:
    int x, y;
};

void f2_const(const Mesh2 &m)
{
    REQUIRE(m.use_count() == 1);
    // rcp_from_this() gives up const version of RCP<Mesh> because 'm' is const
    RCP<const Mesh2> m2 = m.rcp_from_this();
    REQUIRE(m.use_count() == 2);
}

void f2_hybrid(Mesh2 &m)
{
    REQUIRE(m.use_count() == 1);
    // rcp_from_this() gives up const version of RCP<Mesh> even though 'm' is
    // not const, because the internal pointer inside Mesh2 is const.
    RCP<const Mesh2> m2 = m.rcp_from_this();
    REQUIRE(m.use_count() == 2);
}

TEST_CASE("Test rcp_from_this const 2", "[rcp]")
{
    RCP<const Mesh2> m = make_rcp<const Mesh2>();
    REQUIRE(m->use_count() == 1);
    f2_const(*m);
    REQUIRE(m->use_count() == 1);

    RCP<Mesh2> m2 = make_rcp<Mesh2>();
    REQUIRE(m2->use_count() == 1);
    f2_const(*m2);
    REQUIRE(m2->use_count() == 1);
    f2_hybrid(*m2);
    REQUIRE(m2->use_count() == 1);
}

#if defined(WITH_SYMENGINE_COOPERATIVE_INTRUSIVE_RCP)                          \
    || defined(WITH_SYMENGINE_RCP)
TEST_CASE("Test RCP move construct", "[rcp]")
{
    RCP<Mesh> m = make_rcp<Mesh>();
    REQUIRE(m->use_count() == 1);
    RCP<Mesh> m2 = std::move(m);
    REQUIRE(m.is_null());
    REQUIRE(m2->use_count() == 1);
}

TEST_CASE("Test RCP move assign", "[rcp]")
{
    RCP<Mesh> m = make_rcp<Mesh>();
    REQUIRE(m->use_count() == 1);
    RCP<Mesh> m2;
    m2 = std::move(m);
    REQUIRE(m.is_null());
    REQUIRE(m2->use_count() == 1);
}
#endif

TEST_CASE("Test RCP reset", "[rcp]")
{
    RCP<Mesh> m = make_rcp<Mesh>();
    REQUIRE(m->use_count() == 1);
    m.reset();
    REQUIRE(m.is_null());

    // reset on null is safe
    m.reset();
    REQUIRE(m.is_null());
}

TEST_CASE("Test RCP copy does not affect source", "[rcp]")
{
    RCP<Mesh> m = make_rcp<Mesh>();
    REQUIRE(m->use_count() == 1);
    {
        RCP<Mesh> m2 = m;
        REQUIRE(m->use_count() == 2);
        REQUIRE(m2->use_count() == 2);
    }
    REQUIRE(m->use_count() == 1);
}

TEST_CASE("Test RCP field access", "[rcp]")
{
    RCP<Mesh> m = make_rcp<Mesh>();
    m->x = 10;
    m->y = 20;
    REQUIRE(m->x == 10);
    REQUIRE(m->y == 20);
}

TEST_CASE("Test make_rcp const", "[rcp]")
{
    RCP<const Mesh> m = make_rcp<const Mesh>();
    REQUIRE(m->use_count() == 1);
    REQUIRE(not(m == null));
}

TEST_CASE("Test rcp_from_this increments count", "[rcp]")
{
    RCP<Mesh> m = make_rcp<Mesh>();
    REQUIRE(m->use_count() == 1);
    {
        RCP<Mesh> m2 = m->rcp_from_this();
        REQUIRE(m->use_count() == 2);
        REQUIRE(m2->use_count() == 2);
    }
    REQUIRE(m->use_count() == 1);
}

TEST_CASE("Test RCP null comparison", "[rcp]")
{
    RCP<Mesh> m;
    REQUIRE(m == null);
    m = make_rcp<Mesh>();
    REQUIRE(not(m == null));
}

TEST_CASE("Test RCP copy assign to null RCP", "[rcp]")
{
    RCP<Mesh> m;
    RCP<Mesh> m2 = make_rcp<Mesh>();
    m = m2;
    REQUIRE(m->use_count() == 2);
    REQUIRE(m2->use_count() == 2);
}

TEST_CASE("Test RCP copy assign from null", "[rcp]")
{
    RCP<Mesh> m = make_rcp<Mesh>();
    REQUIRE(m->use_count() == 1);
    RCP<Mesh> m2;
    m = m2;
    REQUIRE(m.is_null());
}

#if defined(WITH_SYMENGINE_COOPERATIVE_INTRUSIVE_RCP)
TEST_CASE("Test cooperative_intrusive counter queries", "[rcp]")
{
    RCP<Mesh> m = make_rcp<Mesh>();
    REQUIRE(m->is_uniquely_owned_by_cpp());
    REQUIRE(m->is_uniquely_owned());
    REQUIRE(not m->is_external_owned());

    RCP<Mesh> m2 = m;
    REQUIRE(not m->is_uniquely_owned_by_cpp());
    REQUIRE(not m->is_uniquely_owned());
    REQUIRE(not m->is_external_owned());
}
#endif

TEST_CASE("Test backend-neutral unique ownership query", "[rcp]")
{
    RCP<Mesh> m = make_rcp<Mesh>();
#if defined(WITH_SYMENGINE_RCP) && defined(WITH_SYMENGINE_THREAD_SAFE)
    // Thread-safe legacy RCP intentionally preserves the historical no-steal
    // behavior until that optimization can be reviewed independently.
    REQUIRE(not m->is_uniquely_owned());
#else
    REQUIRE(m->is_uniquely_owned());
    RCP<Mesh> m2 = m;
    REQUIRE(not m->is_uniquely_owned());
#endif
}

// ---------------------------------------------------------------------------
// Immortality. The intrusive backends stop keeping the count; the Teuchos
// backend keeps its count in a separate node and has no equivalent, so
// mark_immortal() is documented there as a no-op and is_immortal() stays
// false. Everything the *library* relies on -- that a marked object survives
// its last RCP and is never a steal candidate -- is asserted per backend.
//
// The retained pointer is parked in a static rather than left on the stack,
// which is not bookkeeping but the point: an immortal object is one a static
// still names, exactly as constants.cpp names `zero` and `one`. That is what
// makes it *reachable* rather than leaked, and it is what a leak checker is
// entitled to be shown.
// ---------------------------------------------------------------------------
static Mesh *immortal_mesh = nullptr;

TEST_CASE("Test mark_immortal", "[rcp]")
{
    RCP<Mesh> m = make_rcp<Mesh>();
    m->x = 4242;
    Mesh *raw = m.get();
    immortal_mesh = raw;
    REQUIRE(not m->is_immortal());

    m->mark_immortal();

#if defined(WITH_SYMENGINE_COOPERATIVE_INTRUSIVE_RCP)                          \
    || defined(WITH_SYMENGINE_RCP)
    REQUIRE(m->is_immortal());
    // Never a steal candidate: no count is kept, so "exactly one reference"
    // is not something that can be established about it.
    REQUIRE(not m->is_uniquely_owned());

    {
        // Copies and destructions leave the count alone -- that is the whole
        // point, and it is why the object below survives.
        RCP<Mesh> m2 = m;
        RCP<Mesh> m3 = m2;
        m3 = m;
        m2.reset();
        REQUIRE(m->is_immortal());
    }

    // Dropping the last handle does not delete an immortal object.
    m.reset();
    REQUIRE(immortal_mesh->x == 4242);
    REQUIRE(immortal_mesh->is_immortal());
    // Deliberately retained: deleting it here would contradict the guarantee
    // the mark makes to anything still holding a raw pointer.
#else
    REQUIRE(not m->is_immortal());
    immortal_mesh = nullptr;
    (void)raw;
#endif
}

TEST_CASE("Test the library's own constants are immortal", "[rcp]")
{
#if defined(WITH_SYMENGINE_COOPERATIVE_INTRUSIVE_RCP)                          \
    || defined(WITH_SYMENGINE_RCP)
    REQUIRE(SymEngine::zero->is_immortal());
    REQUIRE(SymEngine::one->is_immortal());
    REQUIRE(SymEngine::minus_one->is_immortal());
    REQUIRE(SymEngine::two->is_immortal());
    REQUIRE(SymEngine::I->is_immortal());
    REQUIRE(SymEngine::pi->is_immortal());
    REQUIRE(SymEngine::Inf->is_immortal());
    REQUIRE(SymEngine::Nan->is_immortal());
    // ... and are therefore never stolen from by the Add/Mul move-out paths.
    REQUIRE(not SymEngine::one->is_uniquely_owned());
#endif
    // An ordinary object is not immortal in any backend.
    REQUIRE(not SymEngine::integer(7)->is_immortal());
}

// ---------------------------------------------------------------------------
// The immortal band.
//
// Immortality is not a single sentinel value but a range of counts, and the
// reason is ABI: reference traffic is inline, so a consumer compiled against a
// header that predates immortality runs *its* arithmetic on objects this
// library marked. A top-of-range sentinel makes that consumer's first
// increment wrap the count to zero and its next destruction delete a
// process-global singleton.
//
// The band is what makes that harmless, so these assertions are about the band
// itself rather than about any object: that ordinary counts are outside it,
// that the value mark_immortal() writes is far from both edges, and that
// walking a long way from that value in either direction stays inside. The
// last is also what makes the marking race in inc_ref()/dec_ref() benign: a
// loser can only move the count by one.
// ---------------------------------------------------------------------------
#if defined(WITH_SYMENGINE_RCP)
TEST_CASE("Immortal band, original intrusive backend", "[rcp]")
{
    using SymEngine::detail::immortal_refcount;
    using SymEngine::detail::immortal_refcount_hi;
    using SymEngine::detail::immortal_refcount_lo;
    using SymEngine::detail::refcount_is_immortal;

    // Counts a program can actually reach are mortal.
    for (unsigned int c = 0; c < 1024; ++c) {
        REQUIRE(not refcount_is_immortal(c));
    }
    REQUIRE(not refcount_is_immortal(immortal_refcount_lo - 1u));
    REQUIRE(not refcount_is_immortal(immortal_refcount_hi));
    REQUIRE(not refcount_is_immortal(~0u));

    // What mark_immortal() writes, and both edges.
    REQUIRE(refcount_is_immortal(immortal_refcount));
    REQUIRE(refcount_is_immortal(immortal_refcount_lo));
    REQUIRE(refcount_is_immortal(immortal_refcount_hi - 1u));

    // The margin: how far unbalanced traffic from a stale consumer would have
    // to walk before the mark stops being recognised. Spelled as a number
    // rather than a symbol, so that shrinking the band fails here.
    REQUIRE(immortal_refcount - immortal_refcount_lo == 1073741824u);
    REQUIRE(immortal_refcount_hi - immortal_refcount == 1073741824u);
    for (unsigned int k = 1; k <= 1u << 20; k <<= 1) {
        REQUIRE(refcount_is_immortal(immortal_refcount + k));
        REQUIRE(refcount_is_immortal(immortal_refcount - k));
    }
}
#endif

#if defined(WITH_SYMENGINE_COOPERATIVE_INTRUSIVE_RCP)
TEST_CASE("Immortal band, cooperative state word", "[rcp]")
{
    using SymEngine::detail::cooperative_immortal;
    using SymEngine::detail::cooperative_immortal_hi;
    using SymEngine::detail::cooperative_immortal_lo;
    using SymEngine::detail::cooperative_state_is_immortal;

    // C++-owned state words encode count c as 2c+1.
    for (uintptr_t c = 0; c < 1024; ++c) {
        REQUIRE(not cooperative_state_is_immortal(2 * c + 1));
    }
    REQUIRE(not cooperative_state_is_immortal(cooperative_immortal_lo - 2));
    REQUIRE(not cooperative_state_is_immortal(cooperative_immortal_hi));
    REQUIRE(cooperative_state_is_immortal(cooperative_immortal));
    REQUIRE(cooperative_state_is_immortal(cooperative_immortal_lo));
    REQUIRE(cooperative_state_is_immortal(cooperative_immortal_hi - 2));

    // ~0 was the sentinel before the band existed, and is deliberately *not*
    // in it: the band's upper edge has to be an ordinary count so that drift
    // out of the top resumes truthful counting instead of running into the
    // wrap. No state word can hold ~0 -- only mark_immortal() writes a value
    // above a plausible count, and it writes the band's centre.
    REQUIRE(not cooperative_state_is_immortal(~static_cast<uintptr_t>(0)));

    for (uintptr_t k = 2; k <= (uintptr_t)1 << 20; k <<= 1) {
        REQUIRE(cooperative_state_is_immortal(cooperative_immortal + k));
        REQUIRE(cooperative_state_is_immortal(cooperative_immortal - k));
    }
}
#endif

// ---------------------------------------------------------------------------
// mark_immortal() racing reference traffic.
//
// inc_ref() is a load, a test, and a read-modify-write; marking is a store.
// They are separate atomic operations, so a marker can land between the test
// and the RMW, and the RMW then moves the count off the mark. Each operation
// is individually atomic, so **TSan cannot see this** -- there is no data
// race, only a lost update. Nor did the tests above reach it: they mark first
// and then generate traffic, never the other way round.
//
// Two things are asserted here, because neither alone is evidence:
//
//  * The band test above says the *value* consequence is bounded: a loser can
//    move the count by one, and one from the band's centre is still deep
//    inside the band. That is deterministic and it is what makes the race
//    benign rather than merely unlikely. Before the band existed, one past the
//    sentinel was zero.
//  * This test drives the transition itself, with several threads copying and
//    destroying an object while another marks it. It is a probabilistic
//    reproducer and is honest about that: what would make it *fail* to catch a
//    regression is the marking store landing outside every worker's
//    load-to-RMW window in all rounds. It is here because the value argument
//    above is about arithmetic, and this is about the object staying alive.
// ---------------------------------------------------------------------------
#if defined(WITH_SYMENGINE_THREAD_SAFE)                                        \
    && (defined(WITH_SYMENGINE_COOPERATIVE_INTRUSIVE_RCP)                      \
        || defined(WITH_SYMENGINE_RCP))

// An immortal object is never destroyed. That is the property this test
// checks, so it gets its own class with a counting destructor rather than
// reusing Mesh: it makes detection exact instead of leaving it to whether a
// freed object's bytes happen to still look right.
static std::atomic<int> race_mesh_destructions(0);

class RaceMesh : public EnableRCPFromThis<RaceMesh>
{
public:
    int x = 0, y = 0;
    ~RaceMesh()
    {
        ++race_mesh_destructions;
    }
};

// Marked objects are deliberately never freed, so the rounds below park their
// object here for the same reason constants.cpp parks its singletons in
// statics: what makes retention a retention rather than a leak is that
// something still names it. This must have trivial static destruction: a
// vector would discard the roots before LeakSanitizer's exit-time scan.
static constexpr std::size_t race_mesh_rounds = 200;
static std::array<RaceMesh *, race_mesh_rounds> marked_meshes{};

TEST_CASE("mark_immortal races reference traffic", "[rcp]")
{
    const int workers = 4;
    const int ops = 4000;

    race_mesh_destructions.store(0);

    for (std::size_t round = 0; round < race_mesh_rounds; ++round) {
        RCP<RaceMesh> base = make_rcp<RaceMesh>();
        base->x = 0x5eed;
        base->y = static_cast<int>(round);
        marked_meshes[round] = base.get();

        std::atomic<bool> go(false);
        std::atomic<int> ready(0);
        std::vector<std::thread> ts;
        for (int w = 0; w < workers; ++w) {
            ts.emplace_back([&]() {
                ++ready;
                while (!go.load(std::memory_order_acquire)) {
                }
                for (int i = 0; i < ops; ++i) {
                    RCP<RaceMesh> copy = base;  // inc_ref
                    RCP<RaceMesh> copy2 = copy; // inc_ref
                }                               // dec_ref, dec_ref
            });
        }
        while (ready.load() != workers) {
        }
        go.store(true, std::memory_order_release);
        // Into the middle of the storm, not before it and not after it.
        for (int i = 0; i < 64; ++i) {
            std::this_thread::yield();
        }
        base->mark_immortal();
        for (auto &t : ts) {
            t.join();
        }

        REQUIRE(base->is_immortal());
        REQUIRE(base->use_count() == ~0u);
        base.reset();
        // The count is not a reliable witness on its own: an increment that
        // wrapped a top-of-range sentinel to zero used to be *undone* by the
        // next decrement wrapping back to it, so the final value looked right
        // while a delete had already happened in between. The destructor count
        // is the witness that cannot be undone.
        REQUIRE(race_mesh_destructions.load() == 0);
        REQUIRE(marked_meshes[round]->x == 0x5eed);
        REQUIRE(marked_meshes[round]->y == static_cast<int>(round));
        REQUIRE(marked_meshes[round]->is_immortal());
    }

    // Every object every round marked is still there.
    for (std::size_t i = 0; i < marked_meshes.size(); ++i) {
        REQUIRE(marked_meshes[i]->x == 0x5eed);
        REQUIRE(marked_meshes[i]->is_immortal());
    }
    REQUIRE(race_mesh_destructions.load() == 0);
}

#endif
