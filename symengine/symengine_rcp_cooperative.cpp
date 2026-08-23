/*
    Cooperative intrusive reference counter.

    This is an adapted port of nanobind/intrusive/counter.inl:
    Copyright (c) 2023 Wenzel Jakob.
    nanobind is distributed under the BSD 3-Clause License.

*/

#include <symengine/symengine_config.h>

#if defined(WITH_SYMENGINE_COOPERATIVE_INTRUSIVE_RCP)

#include <symengine/symengine_rcp.h>

#include <cstdio>
#include <cstdlib>

namespace
{

// The atomicity policy is defined once, in symengine_rcp.h, and is gated there
// on WITH_SYMENGINE_THREAD_SAFE. These aliases exist so that the slow paths
// below read the same way the fast paths in the header do -- and, more to the
// point, so that there is no second place where a build could decide to be
// atomic or not.
using SymEngine::detail::cooperative_immortal;
using SymEngine::detail::cooperative_immortal_wrapped;
using SymEngine::detail::cooperative_state_cmpxchg;
using SymEngine::detail::cooperative_state_cmpxchg_decrement;
using SymEngine::detail::cooperative_state_is_immortal;
using SymEngine::detail::cooperative_state_load;

#define SYMENGINE_ATOMIC_LOAD(ptr) cooperative_state_load(ptr)
#define SYMENGINE_ATOMIC_CMPXCHG(ptr, cmp, xchg)                               \
    cooperative_state_cmpxchg(ptr, cmp, xchg)
#define SYMENGINE_ATOMIC_DECREMENT(ptr, cmp, xchg)                             \
    cooperative_state_cmpxchg_decrement(ptr, cmp, xchg)

SymEngine::cooperative_incref_hook g_cooperative_incref = nullptr;
SymEngine::cooperative_decref_hook g_cooperative_decref = nullptr;

[[noreturn]] void cooperative_abort(const char *message, const void *ptr)
{
    std::fprintf(stderr, "%s (%p)\n", message, ptr);
    std::abort();
}

void ensure_hooks_ready(const void *ptr)
{
    if (!g_cooperative_incref || !g_cooperative_decref) [[unlikely]] {
        cooperative_abort("cooperative_intrusive backend used before "
                          "cooperative_intrusive_init()",
                          ptr);
    }
}

} // namespace

namespace SymEngine
{

void cooperative_intrusive_init(cooperative_incref_hook incref,
                                cooperative_decref_hook decref) noexcept
{
    g_cooperative_incref = incref;
    g_cooperative_decref = decref;
}

// The two exported entry points. Their bodies are one call each, because the
// logic lives in the header's inlined `*_fast()` pair; they exist so that the
// mangled symbols an already-linked consumer resolved against this library
// keep resolving, and so that anything not on the hot path can keep calling
// the plain names.
void symengine_cooperative_intrusive_counter::inc_ref() const noexcept
{
    inc_ref_fast();
}

bool symengine_cooperative_intrusive_counter::dec_ref() const noexcept
{
    return dec_ref_fast();
}

// The cold halves, entered from the fast paths with the state word they last
// observed. `v` is that word: in the contended case the failed exchange has
// already refreshed it, so resuming the loop here is exactly the loop the
// undivided function used to run.
void symengine_cooperative_intrusive_counter::inc_ref_slow(
    uintptr_t v) const noexcept
{
    while (true) {
        if (v & 1) {
            // An immortal object keeps no count; the fast path normally
            // catches this, but a contended retry can land here. Re-tested at
            // the top of every iteration, so a mark_immortal() that lands
            // after a failed exchange is seen on the next pass.
            if (cooperative_state_is_immortal(v))
                break;
            if (!SYMENGINE_ATOMIC_CMPXCHG(&m_state, &v, v + 2))
                continue;
        } else if (v & cooperative_immortal_wrapped) {
            break;
        } else {
            ensure_hooks_ready(this);
            g_cooperative_incref((void *)v);
        }
        break;
    }
}

bool symengine_cooperative_intrusive_counter::dec_ref_slow(
    uintptr_t v) const noexcept
{
    while (true) {
        if (v & 1) {
            if (cooperative_state_is_immortal(v))
                return false;

            if (v == 1) [[unlikely]] {
                cooperative_abort("symengine_cooperative_intrusive_counter::"
                                  "dec_ref underflow",
                                  this);
            }

            // Non-final decrements release; the final decrement also acquires
            // before its caller deletes (same discipline as dec_ref_fast).
            if (!SYMENGINE_ATOMIC_DECREMENT(&m_state, &v, v - 2))
                continue;

            if (v == 3)
                return true;
        } else if (v & cooperative_immortal_wrapped) {
            return false;
        } else {
            ensure_hooks_ready(this);
            g_cooperative_decref((void *)v);
        }

        return false;
    }
}

void symengine_cooperative_intrusive_counter::mark_immortal() const noexcept
{
    uintptr_t v = SYMENGINE_ATOMIC_LOAD(&m_state);
    while (true) {
        if (!(v & 1)) [[unlikely]] {
            if (v & cooperative_immortal_wrapped)
                return; // already immortal, with a wrapper attached
            cooperative_abort("symengine_cooperative_intrusive_counter::mark_"
                              "immortal on an external-owned object",
                              this);
        }
        // Terminal state: re-tested after every failed exchange, so a second
        // marker (or a reference operation that pushed the count one off the
        // band's centre) does not restart the transition.
        if (cooperative_state_is_immortal(v))
            return;
        if (SYMENGINE_ATOMIC_CMPXCHG(&m_state, &v, cooperative_immortal))
            return;
    }
}

bool symengine_cooperative_intrusive_counter::is_immortal() const noexcept
{
    uintptr_t v = SYMENGINE_ATOMIC_LOAD(&m_state);
    if (v & 1)
        return cooperative_state_is_immortal(v);
    return (v & cooperative_immortal_wrapped) != 0;
}

void symengine_cooperative_intrusive_counter::set_self_external(
    void *o) noexcept
{
    uintptr_t v = SYMENGINE_ATOMIC_LOAD(&m_state);

    if (!(v & 1)) [[unlikely]] {
        cooperative_abort("symengine_cooperative_intrusive_counter::set_self_"
                          "external already external-owned",
                          this);
    }

    // The two low bits of the state word are the union's tags: bit 0
    // distinguishes a count from a pointer, bit 1 marks an immortal object
    // that has a wrapper. So a wrapper pointer has to have both clear. This
    // is checked rather than assumed because getting it wrong is silent: a
    // pointer with bit 1 set would read back as immortal and every reference
    // operation on the object would be dropped. Every runtime this backend
    // targets allocates its objects with at least pointer alignment.
    if (((uintptr_t)o & 3u) != 0) [[unlikely]] {
        cooperative_abort("symengine_cooperative_intrusive_counter::set_self_"
                          "external needs a 4-byte-aligned foreign object",
                          o);
    }

    ensure_hooks_ready(this);

    // Immortality is re-tested after every failed exchange, not only on the
    // initial load. mark_immortal() can land at any point in the loop below,
    // and the band's counts are astronomically large: reading one as an
    // ordinary count would replay ~2^62 foreign increfs, i.e. hang. The
    // sentinel is a terminal state, so recognising it anywhere in the loop
    // means switching protocols, never resuming the count-folding one.
    uintptr_t transferred = 0;
    while (true) {
        if (!(v & 1)) [[unlikely]] {
            // Another externalization won; the count this one was folding is
            // now somebody else's wrapper pointer. Forbidden by the
            // single-external-owner contract, and never silently continued:
            // `v >> 1` on a pointer is not a count.
            cooperative_abort("symengine_cooperative_intrusive_counter::set_"
                              "self_external raced another external owner",
                              this);
        }
        if (cooperative_state_is_immortal(v)) {
            // An immortal object has no count to fold into the wrapper's,
            // because C++ will never hand a reference back: it grants the
            // wrapper exactly one permanent reference, which is what keeps
            // self_external() from ever returning a dangling pointer, and
            // stays immortal (bit 1) so that reference operations remain
            // no-ops. Everything the caller can observe -- self_external(),
            // is_external_owned(), detach_external() -- behaves as it does for
            // an ordinary externalised object.
            //
            // A reference already replayed into the wrapper becomes that
            // permanent one; the rest are released *after* publication, so the
            // wrapper's count never passes through zero.
            if (transferred == 0) {
                g_cooperative_incref(o);
                transferred = 1;
            }
            while (!SYMENGINE_ATOMIC_CMPXCHG(
                &m_state, &v, (uintptr_t)o | cooperative_immortal_wrapped)) {
                // Marking is not reversible and the band is terminal, so the
                // only way the state can have changed away from it is another
                // set_self_external() racing this one -- which the
                // single-external-owner contract already forbids.
                if (!cooperative_state_is_immortal(v)) {
                    cooperative_abort("symengine_cooperative_intrusive_counter"
                                      "::set_self_external raced another "
                                      "external owner",
                                      this);
                }
            }
            while (transferred > 1) {
                g_cooperative_decref(o);
                --transferred;
            }
            return;
        }

        for (uintptr_t count = v >> 1; transferred < count; ++transferred)
            g_cooperative_incref(o);
        if (SYMENGINE_ATOMIC_CMPXCHG(&m_state, &v, (uintptr_t)o))
            break;
    }

    for (uintptr_t count = v >> 1; transferred > count; --transferred)
        g_cooperative_decref(o);
}

void *symengine_cooperative_intrusive_counter::self_external() const noexcept
{
    uintptr_t v = SYMENGINE_ATOMIC_LOAD(&m_state);
    if (v & 1)
        return nullptr;
    return (void *)(v & ~cooperative_immortal_wrapped);
}

unsigned int symengine_cooperative_intrusive_counter::use_count() const noexcept
{
    uintptr_t v = SYMENGINE_ATOMIC_LOAD(&m_state);
    // An immortal object reports UINT_MAX: no count is kept, and "effectively
    // infinite" is the only honest answer to "how many references are there?"
    // for something that is never released. Spelled out rather than left to
    // `v >> 1`, which would leak the band's biased value.
    if ((v & 1) && cooperative_state_is_immortal(v))
        return ~0u;
    if (v & 1)
        return static_cast<unsigned int>(v >> 1);
    return 0;
}

bool symengine_cooperative_intrusive_counter::is_external_owned() const noexcept
{
    return (SYMENGINE_ATOMIC_LOAD(&m_state) & 1) == 0;
}

bool symengine_cooperative_intrusive_counter::is_uniquely_owned_by_cpp()
    const noexcept
{
    uintptr_t v = SYMENGINE_ATOMIC_LOAD(&m_state);
    return (v & 1) && (v >> 1) == 1;
}

void *symengine_cooperative_intrusive_counter::detach_external() const noexcept
{
    uintptr_t v = SYMENGINE_ATOMIC_LOAD(&m_state);

    // Shutdown normally provides exclusive access, but retain the same CAS
    // discipline as the other mode transitions so a stale load cannot reset a
    // state that changed before publication. The shutdown caller still owns
    // synchronization with foreign-runtime hook calls.
    while (true) {
        if (v & 1)
            return nullptr;
        // An immortal object goes back to being immortal without a wrapper,
        // not to a count of zero: it is still never deleted, and a count of
        // zero would make the next dec_ref() an underflow abort.
        const bool immortal = (v & cooperative_immortal_wrapped) != 0;
        const uintptr_t next = immortal ? cooperative_immortal : (uintptr_t)1;
        if (SYMENGINE_ATOMIC_CMPXCHG(&m_state, &v, next))
            return (void *)(v & ~cooperative_immortal_wrapped);
    }
}

} // namespace SymEngine

#endif
