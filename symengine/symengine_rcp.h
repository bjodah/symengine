#ifndef SYMENGINE_RCP_H
#define SYMENGINE_RCP_H

#include <iostream>
#include <cstddef>
#include <stdexcept>
#include <string>
#if __cplusplus <= 201703L
#include <ciso646>
#else
#include <version>
#endif

#include <symengine/symengine_config.h>
#include <symengine/symengine_assert.h>

#if defined(WITH_SYMENGINE_COOPERATIVE_INTRUSIVE_RCP)

#include <cstdint>

#elif defined(WITH_SYMENGINE_RCP)

#if defined(WITH_SYMENGINE_THREAD_SAFE)
#include <atomic>
#endif

#else

// Include all Teuchos headers here:
#include <symengine/utilities/teuchos/Teuchos_RCP.hpp>
#include <symengine/utilities/teuchos/Teuchos_TypeNameTraits.hpp>

#endif

// A hint, not a promise, and only about code layout: reference operations on
// ordinary objects greatly outnumber those on the handful of immortal
// singletons, and told nothing GCC lays the atomic read-modify-write out as
// the cold path -- costing every ordinary increment a taken branch and a jump
// back. The builtin applies the bias directly to the tested expression and is
// available on both compiler families used here. Undefined at the end of the
// file.
#if defined(__GNUC__) || defined(__clang__)
#define SYMENGINE_RCP_MORTAL(x) __builtin_expect(!!(x), 1)
#else
#define SYMENGINE_RCP_MORTAL(x) (x)
#endif

namespace SymEngine
{

// ---------------------------------------------------------------------------
// noexcept became part of a function-pointer *type* in C++17, so the two
// branches below declare different types -- and therefore different mangled
// names for cooperative_intrusive_init(), `PDoFvPvE` against `PFvPvE`. This
// library declares cxx_std_11, so both branches are reachable, and a foreign
// runtime's binding must be compiled against the same side of this test as
// the library it links: a mismatch is an undefined symbol at link time, not
// a silent disagreement.
#if __cplusplus >= 201703L
using cooperative_incref_hook = void (*)(void *) noexcept;
using cooperative_decref_hook = void (*)(void *) noexcept;
#else
using cooperative_incref_hook = void (*)(void *);
using cooperative_decref_hook = void (*)(void *);
#endif

void cooperative_intrusive_init(cooperative_incref_hook incref,
                                cooperative_decref_hook decref) noexcept;

// Counter wrapper for the cooperative_intrusive backend.
// ---------------------------------------------------------------------------
#if defined(WITH_SYMENGINE_COOPERATIVE_INTRUSIVE_RCP)

namespace detail
{
// ---------------------------------------------------------------------------
// The cooperative counter's atomicity policy, defined **once**.
//
// Both the inlined fast paths below and the out-of-line slow paths in
// symengine_rcp_cooperative.cpp go through these primitives, so the two
// halves of the counter cannot disagree about whether this build is
// thread-safe. That single definition is the whole safety argument: a
// ThreadSanitizer run cannot police the choice, because the library a
// sanitized consumer links is not itself instrumented, so the gate has to be
// right by construction rather than by test.
//
// WITH_SYMENGINE_THREAD_SAFE lives in the generated, installed
// symengine_config.h, which this header includes. A consumer therefore always
// compiles the inlined fast path with the same policy the library was built
// with; there is no way to select one here and the other there.
//
// The compare-and-exchange is load-bearing in *both* configurations and is
// not merely an atomicity device: bit 0 of the state word is the tag of a
// union (set = C++-owned reference count in the upper bits, clear = pointer
// to a foreign-runtime object), so every mode transition has to decide on
// exactly the word it tested. The non-thread-safe build keeps that discipline
// and drops only the atomicity, which turns the sequence into a plain
// load / test / add / store.
// ---------------------------------------------------------------------------
#if defined(WITH_SYMENGINE_THREAD_SAFE)

#if !defined(_MSC_VER)

inline uintptr_t cooperative_state_load(const uintptr_t *ptr) noexcept
{
    return __atomic_load_n(ptr, __ATOMIC_RELAXED);
}

inline bool cooperative_state_cmpxchg(uintptr_t *ptr, uintptr_t *cmp,
                                      uintptr_t xchg) noexcept
{
    return __atomic_compare_exchange_n(ptr, cmp, xchg, true, __ATOMIC_RELAXED,
                                       __ATOMIC_RELAXED);
}

// A decrement releases; the final encoded `3 -> 1` decrement also acquires
// before its caller deletes. Keeping both halves on the CAS is equivalent to a
// release CAS plus an acquire fence, and is directly understood by TSan.
// Increments stay relaxed: this is the classic asymmetric refcount discipline.
inline bool cooperative_state_cmpxchg_decrement(uintptr_t *ptr, uintptr_t *cmp,
                                                uintptr_t xchg) noexcept
{
    if (*cmp == 3) {
        return __atomic_compare_exchange_n(ptr, cmp, xchg, true,
                                           __ATOMIC_ACQ_REL, __ATOMIC_RELAXED);
    }
    return __atomic_compare_exchange_n(ptr, cmp, xchg, true, __ATOMIC_RELEASE,
                                       __ATOMIC_RELAXED);
}

#else // _MSC_VER

extern "C" void *_InterlockedCompareExchangePointer(void *volatile *Destination,
                                                    void *Exchange,
                                                    void *Comparand);
#pragma intrinsic(_InterlockedCompareExchangePointer)

inline uintptr_t cooperative_state_load(const uintptr_t *ptr) noexcept
{
    return *((volatile const uintptr_t *)ptr);
}

inline bool cooperative_state_cmpxchg(uintptr_t *ptr, uintptr_t *cmp,
                                      uintptr_t xchg) noexcept
{
    uintptr_t cmpv = *cmp;
    uintptr_t prev = (uintptr_t)_InterlockedCompareExchangePointer(
        (void *volatile *)ptr, (void *)xchg, (void *)cmpv);
    if (prev == cmpv) {
        return true;
    }
    *cmp = prev;
    return false;
}

// _InterlockedCompareExchangePointer is a full barrier on every MSVC target,
// so it supplies both decrement orderings without another fence.
inline bool cooperative_state_cmpxchg_decrement(uintptr_t *ptr, uintptr_t *cmp,
                                                uintptr_t xchg) noexcept
{
    return cooperative_state_cmpxchg(ptr, cmp, xchg);
}

#endif // _MSC_VER

#else // !WITH_SYMENGINE_THREAD_SAFE

// This build has declared that it does not share SymEngine objects between
// threads. The tag-bit dispatch above is preserved exactly; only the atomicity
// is dropped.
inline uintptr_t cooperative_state_load(const uintptr_t *ptr) noexcept
{
    return *ptr;
}

inline bool cooperative_state_cmpxchg(uintptr_t *ptr, uintptr_t *cmp,
                                      uintptr_t xchg) noexcept
{
    if (*ptr == *cmp) {
        *ptr = xchg;
        return true;
    }
    // Cannot happen without concurrency, but the caller's retry loops are
    // written against this contract and the compiler folds the branch away.
    *cmp = *ptr;
    return false;
}

// No threads, no ordering: the decrement variant exists so the callers read
// identically in every configuration.
inline bool cooperative_state_cmpxchg_decrement(uintptr_t *ptr, uintptr_t *cmp,
                                                uintptr_t xchg) noexcept
{
    return cooperative_state_cmpxchg(ptr, cmp, xchg);
}

#endif // WITH_SYMENGINE_THREAD_SAFE

// ---------------------------------------------------------------------------
// Immortality.
//
// An object that is never destroyed does not need its reference count kept:
// no sequence of increments and decrements can change the answer to "is it
// still alive?". Not keeping it is worth a great deal on a shared object,
// because a reference count is a *write*, and a write to a line several cores
// hold forces that line through the coherence protocol every single time. The
// library's process-global singletons (constants.cpp: `zero`, `one`,
// `minus_one`, ...) are exactly such objects, and are referenced by nearly
// every expression built -- `Mul` stores `one` in its coefficient slot, `Add`
// stores `zero`, `sub()` multiplies by `minus_one`. See
// benchmarks/rcp_scaling.cpp for what that costs.
//
// Immortality is a *band* of count values, not one sentinel. See the
// "immortal band" comment below cooperative_immortal_wrapped for why -- in
// short, because reference traffic is inline, and code compiled against an
// older header does not know the mark exists.
//
// Two state-word encodings:
//
//   cooperative_immortal         bit 0 set, so the object stays "C++-owned"
//                                for every existing predicate. The value is
//                                the centre of the immortal band and is not a
//                                plausible count. `use_count()` reports
//                                UINT_MAX and `is_uniquely_owned_by_cpp()`
//                                reports false, which is what a caller asking
//                                "may I steal this?" needs to hear about a
//                                singleton.
//
//   cooperative_immortal_wrapped
//                                bit 1 of an *external* (bit 0 clear) state
//                                word: an immortal object that a foreign
//                                runtime has attached a wrapper to. C++
//                                reference operations stay no-ops -- the
//                                foreign runtime holds one permanent
//                                reference and C++ never gives it back --
//                                so an immortal object keeps its scaling
//                                even under a language binding. This costs
//                                bit 1 of the wrapper pointer, i.e. foreign
//                                objects must be 4-byte aligned; every
//                                runtime this backend targets allocates
//                                its objects with at least pointer
//                                alignment, and set_self_external() checks.
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// The immortal band, and why immortality is not a single sentinel.
//
// Reference-count traffic is inline: `RCP<T>`'s copy and destruction reach the
// arithmetic below through a header, so the code that runs is whatever header
// the *consumer* was compiled against. The library's SONAME does not change
// when a purely additive change like this one lands, so ordinary shared-object
// replacement pairs an old consumer's inline arithmetic with a new library's
// marked objects. An old consumer does not know the mark exists; it increments
// and decrements the sentinel like any other count.
//
// A single top-of-range sentinel makes that fatal. `~0u + 1` is 0, `0 + 1` is
// 1, and the next destruction deletes a process-global singleton. So the sentinel
// is chosen to be a value old inline code handles *harmlessly*: a large biased
// count, tested by range rather than by equality.
//
//   count < lo            ordinary, mortal, counted by everybody
//   lo <= count < hi      immortal: this header elides the read-modify-write
//   count >= hi           ordinary again
//
// mark_immortal() writes the band's centre. An old consumer's traffic is
// RAII-balanced, so it hovers there and never approaches zero: correct, merely
// un-optimised, which is exactly what that consumer had before.
//
// The objection to answer is drift, since a mixed process can produce
// *unbalanced* operations -- a reference the new library creates (elided) and
// the old consumer destroys (counted) is a net decrement, and the mirror image
// is a net increment. Leaving the band takes `hi - centre` of them, which is
// 2^30 for the 32-bit count below and 2^61 for the cooperative state word.
// But the two edges are *absorbing*, which is the property that matters more
// than the size of the margin: the moment a count leaves the band this header
// resumes counting truthfully, so the very traffic that was unbalanced becomes
// balanced again and the count cannot walk further. Past the low edge it can
// only fall by the number of *simultaneously live* references whose increment
// was elided -- and 2^30 live references is 8 GB of pointers, the same order
// of impossibility the unmodified counter already relies on to not overflow.
//
// The cost is that an ordinary object holding >= `lo` references would be
// misread as immortal. That is the same trade the old sentinel made at `~0u`,
// moved down by four; 2^30 references cannot fit in an address space that also
// has to hold the objects making them.
//
// Values, for a 64-bit uintptr_t (the count is v >> 1, so the band's counts
// are 2^61, 2^62 and 3*2^61):
const uintptr_t cooperative_immortal_lo
    = (static_cast<uintptr_t>(1) << (sizeof(uintptr_t) * 8 - 2)) | 1;
const uintptr_t cooperative_immortal
    = (static_cast<uintptr_t>(1) << (sizeof(uintptr_t) * 8 - 1)) | 1;
const uintptr_t cooperative_immortal_hi
    = (static_cast<uintptr_t>(3) << (sizeof(uintptr_t) * 8 - 2)) | 1;
const uintptr_t cooperative_immortal_wrapped = 2;

//! \return whether the C++-owned state word `v` is in the immortal band.
//! The caller must already have established `(v & 1) != 0`; an external state
//! word is a pointer, and pointers are not counts.
inline bool cooperative_state_is_immortal(uintptr_t v) noexcept
{
    return static_cast<uintptr_t>(v - cooperative_immortal_lo)
           < static_cast<uintptr_t>(cooperative_immortal_hi
                                    - cooperative_immortal_lo);
}

} // namespace detail

//! Counter wrapper for the cooperative_intrusive backend.
//! In external-owned mode, use_count() returns 0 because the foreign runtime
//! holds the live references. Use is_uniquely_owned_by_cpp() for ownership-
//! gated optimizations.
class symengine_cooperative_intrusive_counter
{
public:
    symengine_cooperative_intrusive_counter() noexcept : m_state(1) {}

    // inc_ref()/dec_ref() are the ABI: they keep their exported, out-of-line
    // definitions in symengine_rcp_cooperative.cpp so that anything already
    // linked against this library keeps resolving them. They forward to the
    // *_fast() pair below, which is the single definition of the logic.
    void inc_ref() const noexcept;
    bool dec_ref() const noexcept;

    //! The hot path, inlined: one load, one tag test, one CAS.
    //!
    //! This is what `RCP<T>` copy and destruction reach (via
    //! `EnableRCPFromThis`), and inlining it here is the point of the split --
    //! out of line it was a cross-DSO PLT call on every reference-count
    //! operation, which a profile of the Rubi corpus put at 17% of run time.
    //! Everything that is not "the object is C++-owned and the count moves by
    //! one" goes out of line to `*_slow`, which owns the retry loop, the
    //! foreign-runtime hooks and the underflow abort.
    void inc_ref_fast() const noexcept
    {
        uintptr_t v = detail::cooperative_state_load(&m_state);
        if ((v & 1) != 0) {
            // Immortal: return without writing anything. That is the whole
            // point -- the line stays Shared in every core's cache instead of
            // being pulled Modified into one of them.
            if (!SYMENGINE_RCP_MORTAL(
                    !detail::cooperative_state_is_immortal(v))) {
                return;
            }
            if (detail::cooperative_state_cmpxchg(&m_state, &v, v + 2)) {
                return;
            }
        } else if ((v & detail::cooperative_immortal_wrapped) != 0) {
            return;
        }
        inc_ref_slow(v);
    }

    //! \return true when this call dropped the last C++ reference.
    //! Non-final drops release; the final drop also acquires before deletion.
    bool dec_ref_fast() const noexcept
    {
        uintptr_t v = detail::cooperative_state_load(&m_state);
        if ((v & 1) != 0) {
            if (!SYMENGINE_RCP_MORTAL(
                    !detail::cooperative_state_is_immortal(v))) {
                return false;
            }
            // `v == 1` is a count of zero, i.e. an underflow; it is handed to
            // the slow path so the abort lives in one place.
            if (v != 1
                && detail::cooperative_state_cmpxchg_decrement(&m_state, &v,
                                                               v - 2)) {
                return v == 3;
            }
        } else if ((v & detail::cooperative_immortal_wrapped) != 0) {
            return false;
        }
        return dec_ref_slow(v);
    }

    //! Declare that this object outlives the process: stop keeping its
    //! reference count altogether. Increments and decrements become no-ops
    //! and it is never deleted, so the storage is retained until exit.
    //!
    //! Only correct for genuine process-global singletons. The current count
    //! is *discarded*, so anything still holding a reference keeps a valid
    //! pointer forever -- which is the guarantee, not a leak of the ordinary
    //! kind, but it does mean an object marked by mistake is never freed.
    //!
    //! Must be called while the object is C++-owned; marking one a foreign
    //! runtime already owns aborts, since the count has by then been folded
    //! into the wrapper's own and cannot be discarded on its behalf.
    //!
    //! The caller must also hold a reference. Marking is atomic, but nothing
    //! can rescue a call that races the drop of the last one: the other thread
    //! has already been told to delete.
    void mark_immortal() const noexcept;
    bool is_immortal() const noexcept;

    //! Hand ownership to a foreign runtime. `o` must be at least 4-byte
    //! aligned: the two low bits of the state word are the union's tags (see
    //! detail::cooperative_immortal above), and a wrapper pointer that
    //! collided with either would be misread. Checked, not assumed.
    void set_self_external(void *o) noexcept;
    void *self_external() const noexcept;
    unsigned int use_count() const noexcept;
    bool is_external_owned() const noexcept;
    bool is_uniquely_owned_by_cpp() const noexcept;

    // Detach the foreign wrapper: atomically reset to C++-owned (refcount 0)
    // and return the previous runtime object pointer so the caller can
    // release it in the foreign runtime. This is used during runtime shutdown.
    // The CAS prevents a stale state snapshot from overwriting a concurrent
    // mode transition; callers must still synchronize shutdown with hook calls.
    // Returns nullptr if already C++-owned.
    void *detach_external() const noexcept;

private:
    // The paths the inlined fast ones decline: foreign-runtime ownership, a
    // contended compare-and-exchange, and the underflow abort. Out of line
    // because they need the process-wide hook pointers, and because they are
    // cold by construction.
    void inc_ref_slow(uintptr_t v) const noexcept;
    bool dec_ref_slow(uintptr_t v) const noexcept;

    mutable uintptr_t m_state;
};
static_assert(
    sizeof(symengine_cooperative_intrusive_counter) == sizeof(void *),
    "cooperative-intrusive counter must stay pointer-sized for the C ABI");
#endif

#if defined(WITH_SYMENGINE_RCP)
namespace detail
{
// The immortal band for the original intrusive backend, whose count is a bare
// `unsigned int`. Same construction and same reasoning as the cooperative
// state word above (see "The immortal band" there): a biased count tested by
// range, so that inline arithmetic compiled against a header that predates
// immortality keeps a marked object alive instead of wrapping its count to
// zero and deleting a process-global singleton.
//
//   [2^30, 3*2^30)   the band          mark_immortal() writes 2^31, its centre
//
// so 2^30 net unbalanced operations are needed to leave it in either
// direction, and both edges are absorbing: outside the band this header counts
// truthfully again, which is what makes the drift self-limiting rather than a
// slow walk to zero.
static_assert(
    sizeof(unsigned int) >= 4,
    "the immortal band assumes a reference count of at least 32 bits");
enum : unsigned int {
    immortal_refcount_lo = 1u << 30,
    immortal_refcount = 1u << 31,
    immortal_refcount_hi = 3u << 30
};

//! \return whether the reference count `c` is in the immortal band.
inline bool refcount_is_immortal(unsigned int c)
{
    return static_cast<unsigned int>(c - immortal_refcount_lo)
           < static_cast<unsigned int>(immortal_refcount_hi
                                       - immortal_refcount_lo);
}
} // namespace detail
#endif

// ---------------------------------------------------------------------------
// Shared between both intrusive backends (cooperative_intrusive and symengine):
// Ptr<T>, ENull, rcp(), rcp_*_cast, outArg, ptrFromRef, typeName,
// print_stack_on_segfault.
// ---------------------------------------------------------------------------
#if defined(WITH_SYMENGINE_COOPERATIVE_INTRUSIVE_RCP)                          \
    || defined(WITH_SYMENGINE_RCP)

/* Ptr */

// Ptr is always pointing to a valid object (can never be nullptr).

template <class T>
class Ptr
{
public:
    inline explicit Ptr(T *ptr) : ptr_(ptr)
    {
        SYMENGINE_ASSERT(ptr_ != nullptr)
    }
    inline Ptr(const Ptr<T> &ptr) : ptr_(ptr.ptr_) {}
    template <class T2>
    inline Ptr(const Ptr<T2> &ptr) : ptr_(ptr.get())
    {
    }
    Ptr<T> &operator=(const Ptr<T> &ptr)
    {
        ptr_ = ptr.get();
        return *this;
    }
#if defined(HAVE_DEFAULT_CONSTRUCTORS)
    inline Ptr(Ptr &&) = default;
    Ptr<T> &operator=(Ptr &&) = default;
#endif
    inline T *operator->() const
    {
        return ptr_;
    }
    inline T &operator*() const
    {
        return *ptr_;
    }
    inline T *get() const
    {
        return ptr_;
    }
    inline T *getRawPtr() const
    {
        return get();
    }
    inline const Ptr<T> ptr() const
    {
        return *this;
    }

private:
    T *ptr_;
};

template <typename T>
inline Ptr<T> outArg(T &arg)
{
    return Ptr<T>(&arg);
}

/** \brief Create a pointer to a object from an object reference.
 *
 * \relates Ptr
 */
template <typename T>
inline Ptr<T> ptrFromRef(T &arg)
{
    return Ptr<T>(&arg);
}

/* RCP */

enum ENull { null };

// RCP can be null. Functionally it should be equivalent to Teuchos::RCP.

// ---------------------------------------------------------------------------
// RCP<T> is shared by both intrusive backends. It reaches the count only
// through the object's own inc_ref()/dec_ref(), which EnableRCPFromThis<T>
// defines per backend -- so the pointer semantics below, and the branches
// that decide when to delete, are the same code for cooperative_intrusive as
// for symengine, and are written once.
// ---------------------------------------------------------------------------
template <class T>
class RCP
{
public:
    RCP(ENull null_arg = null) : ptr_(nullptr) {}
    explicit RCP(T *p) : ptr_(p)
    {
        SYMENGINE_ASSERT(ptr_ != nullptr)
        ptr_->inc_ref();
    }
    // Copy constructor
    RCP(const RCP<T> &rp) : ptr_(rp.ptr_)
    {
        if (not is_null())
            ptr_->inc_ref();
    }
    // Copy constructor
    template <class T2>
    RCP(const RCP<T2> &r_ptr) : ptr_(r_ptr.get())
    {
        if (not is_null())
            ptr_->inc_ref();
    }
    // Move constructor
    RCP(RCP<T> &&rp) SYMENGINE_NOEXCEPT : ptr_(rp.ptr_)
    {
        rp.ptr_ = nullptr;
    }
    // Move constructor
    template <class T2>
    RCP(RCP<T2> &&r_ptr)
    SYMENGINE_NOEXCEPT : ptr_(r_ptr.get())
    {
        r_ptr._set_null();
    }
    ~RCP() SYMENGINE_NOEXCEPT
    {
        if (ptr_ != nullptr and ptr_->dec_ref())
            delete ptr_;
    }
    T *operator->() const
    {
        SYMENGINE_ASSERT(ptr_ != nullptr)
        return ptr_;
    }
    T &operator*() const
    {
        SYMENGINE_ASSERT(ptr_ != nullptr)
        return *ptr_;
    }
    T *get() const
    {
        return ptr_;
    }
    Ptr<T> ptr() const
    {
        return Ptr<T>(get());
    }
    bool is_null() const
    {
        return ptr_ == nullptr;
    }
    template <class T2>
    bool operator==(const RCP<T2> &p2) const
    {
        return ptr_ == p2.ptr_;
    }
    template <class T2>
    bool operator!=(const RCP<T2> &p2) const
    {
        return ptr_ != p2.ptr_;
    }
    // Copy assignment
    RCP<T> &operator=(const RCP<T> &r_ptr)
    {
        T *r_ptr_ptr_ = r_ptr.ptr_;
        if (not r_ptr.is_null())
            r_ptr_ptr_->inc_ref();
        if (not is_null() and ptr_->dec_ref())
            delete ptr_;
        ptr_ = r_ptr_ptr_;
        return *this;
    }
    // Move assignment
    RCP<T> &operator=(RCP<T> &&r_ptr)
    {
        std::swap(ptr_, r_ptr.ptr_);
        return *this;
    }
    void reset()
    {
        if (not is_null() and ptr_->dec_ref())
            delete ptr_;
        ptr_ = nullptr;
    }
    // Don't use this function directly:
    void _set_null()
    {
        ptr_ = nullptr;
    }

private:
    T *ptr_;
};

template <class T>
inline RCP<T> rcp(T *p)
{
    return RCP<T>(p);
}

template <class T2, class T1>
inline RCP<T2> rcp_static_cast(const RCP<T1> &p1)
{
    // Make the compiler check if the conversion is legal
    T2 *check = static_cast<T2 *>(p1.get());
    return RCP<T2>(check);
}

template <class T2, class T1>
inline RCP<T2> rcp_dynamic_cast(const RCP<T1> &p1)
{
    if (not p1.is_null()) {
        T2 *p = nullptr;
        // Make the compiler check if the conversion is legal
        p = dynamic_cast<T2 *>(p1.get());
        if (p) {
            return RCP<T2>(p);
        }
    }
    throw std::runtime_error("rcp_dynamic_cast: cannot convert.");
}

template <class T2, class T1>
inline RCP<T2> rcp_const_cast(const RCP<T1> &p1)
{
    // Make the compiler check if the conversion is legal
    T2 *check = const_cast<T2 *>(p1.get());
    return RCP<T2>(check);
}

template <class T>
inline bool operator==(const RCP<T> &p, ENull)
{
    return p.get() == nullptr;
}

template <typename T>
std::string typeName(const T &t)
{
    return "RCP<>";
}

void print_stack_on_segfault();

#endif // shared intrusive block

// ---------------------------------------------------------------------------
// Teuchos backend aliases
// ---------------------------------------------------------------------------
#if !defined(WITH_SYMENGINE_COOPERATIVE_INTRUSIVE_RCP)                         \
    && !defined(WITH_SYMENGINE_RCP)

using Teuchos::null;
using Teuchos::outArg;
using Teuchos::print_stack_on_segfault;
using Teuchos::Ptr;
using Teuchos::ptrFromRef;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcp_static_cast;
using Teuchos::typeName;

#endif

template <class T>
class EnableRCPFromThis
{
    // Public interface
public:
    //! Get RCP<T> pointer to self (it will cast the pointer to T)
    inline RCP<T> rcp_from_this()
    {
#if defined(WITH_SYMENGINE_COOPERATIVE_INTRUSIVE_RCP)                          \
    || defined(WITH_SYMENGINE_RCP)
        return rcp(static_cast<T *>(this));
#else
        return rcp_static_cast<T>(weak_self_ptr_.create_strong());
#endif
    }

    //! Get RCP<const T> pointer to self (it will cast the pointer to const T)
    inline RCP<const T> rcp_from_this() const
    {
#if defined(WITH_SYMENGINE_COOPERATIVE_INTRUSIVE_RCP)                          \
    || defined(WITH_SYMENGINE_RCP)
        return rcp(static_cast<const T *>(this));
#else
        return rcp_static_cast<const T>(weak_self_ptr_.create_strong());
#endif
    }

    //! Get RCP<T2> pointer to self (it will cast the pointer to T2)
    template <class T2>
    inline RCP<const T2> rcp_from_this_cast() const
    {
#if defined(WITH_SYMENGINE_COOPERATIVE_INTRUSIVE_RCP)                          \
    || defined(WITH_SYMENGINE_RCP)
        return rcp(static_cast<const T2 *>(this));
#else
        return rcp_static_cast<const T2>(weak_self_ptr_.create_strong());
#endif
    }

    //! In the cooperative_intrusive backend, use_count() returns 0 for
    //! external-owned objects, so it is NOT a total reference count in that
    //! mode. In every backend the value is advisory while references can be
    //! changed concurrently; use is_uniquely_owned() for ownership-gated
    //! optimizations.
    unsigned int use_count() const
    {
#if defined(WITH_SYMENGINE_COOPERATIVE_INTRUSIVE_RCP)
        return refcount_.use_count();
#elif defined(WITH_SYMENGINE_RCP)
        // An immortal object reports UINT_MAX in every backend: no count is
        // kept, so "effectively infinite" is the only honest answer, and the
        // biased value the band actually holds is an encoding detail no caller
        // should see.
        {
            const unsigned int c = refcount_;
            return detail::refcount_is_immortal(c) ? ~0u : c;
        }
#else
        return weak_self_ptr_.strong_count();
#endif
    }

    //! Whether this object can be safely treated as uniquely C++-owned for
    //! move-out optimizations. Externally-owned cooperative objects are never
    //! unique: a foreign runtime can still reach their wrapper. Neither is an
    //! immortal one -- its count is not kept, so "one reference" is never
    //! something that can be established about it.
    bool is_uniquely_owned() const noexcept
    {
#if defined(WITH_SYMENGINE_COOPERATIVE_INTRUSIVE_RCP)
        return refcount_.is_uniquely_owned_by_cpp();
#elif defined(WITH_SYMENGINE_RCP)
#if defined(WITH_SYMENGINE_THREAD_SAFE)
        // Keep the historical no-steal behavior in thread-safe builds. A
        // relaxed load could be considered in a future, separately reviewed
        // optimization, but it is intentionally not enabled here.
        return false;
#else
        return refcount_ == 1;
#endif
#else
        return weak_self_ptr_.strong_count() == 1;
#endif
    }

    //! Declare that this object outlives the process, so that its reference
    //! count need not be kept at all: increments and decrements become no-ops
    //! and it is never deleted. See detail::cooperative_immortal in this
    //! header, and benchmarks/rcp_scaling.cpp for the motivation -- a shared
    //! object's reference count is a write to a line every core wants, and
    //! that is what stops a thread-parallel workload from scaling.
    //!
    //! Only correct for genuine process-global singletons: the current count
    //! is discarded and the storage is retained until exit. The Teuchos
    //! backend keeps its count in a separate node rather than in the object,
    //! and has no equivalent, so this is a no-op there.
    void mark_immortal() const noexcept
    {
#if defined(WITH_SYMENGINE_COOPERATIVE_INTRUSIVE_RCP)
        refcount_.mark_immortal();
#elif defined(WITH_SYMENGINE_RCP)
        // The band's centre, not a top-of-range sentinel: see
        // detail::refcount_is_immortal. Plain store rather than an exchange --
        // the band is terminal, and inc_ref()/dec_ref() repair the ±1 a
        // concurrent loser can leave behind.
        refcount_ = detail::immortal_refcount;
#endif
    }

    //! \return whether mark_immortal() has been called on this object.
    bool is_immortal() const noexcept
    {
#if defined(WITH_SYMENGINE_COOPERATIVE_INTRUSIVE_RCP)
        return refcount_.is_immortal();
#elif defined(WITH_SYMENGINE_RCP)
        return detail::refcount_is_immortal(refcount_);
#else
        return false;
#endif
    }

#if defined(WITH_SYMENGINE_COOPERATIVE_INTRUSIVE_RCP)

    // Everything below is private interface for cooperative_intrusive
private:
    mutable symengine_cooperative_intrusive_counter refcount_;

public:
    EnableRCPFromThis() = default;

    // The *_fast() pair rather than inc_ref()/dec_ref(): this is the hot path
    // (every RCP copy and destruction lands here) and the fast variants are
    // inlined, where the ABI-stable ones are an out-of-line call.
    void inc_ref() const noexcept
    {
        refcount_.inc_ref_fast();
    }
    bool dec_ref() const noexcept
    {
        return refcount_.dec_ref_fast();
    }
    void set_self_external(void *o) const noexcept
    {
        refcount_.set_self_external(o);
    }
    void *self_external() const noexcept
    {
        return refcount_.self_external();
    }

    bool is_external_owned() const noexcept
    {
        return refcount_.is_external_owned();
    }
    bool is_uniquely_owned_by_cpp() const noexcept
    {
        return refcount_.is_uniquely_owned_by_cpp();
    }

    // Detach the foreign wrapper: atomically reset to C++-owned (refcount 0)
    // and return the previous runtime object pointer so the caller can release
    // it in the foreign runtime. Used by runtime-specific shutdown cleanup.
    void *detach_external() const noexcept
    {
        return refcount_.detach_external();
    }

private:
#elif defined(WITH_SYMENGINE_RCP)

    // Everything below is private interface for symengine RCP
private:
//! Public variables if defined with SYMENGINE_RCP
// The reference counter is defined either as "unsigned int" (faster, but
// not thread safe) or as std::atomic<unsigned int> (slower, but thread
// safe). Semantically they are almost equivalent, except that the
// pre-decrement operator `operator--()` returns a copy for std::atomic
// instead of a reference to itself.
// The refcount_ is defined as mutable, because it does not change the
// state of the instance, but changes when more copies
// of the same instance are made.
#if defined(WITH_SYMENGINE_THREAD_SAFE)
    mutable std::atomic<unsigned int> refcount_; // reference counter
#else
    mutable unsigned int refcount_; // reference counter
#endif // WITH_SYMENGINE_THREAD_SAFE

public:
    EnableRCPFromThis() : refcount_(0) {}

    // The reference operations RCP<T> uses. They are one branch wider than
    // the bare `refcount_++` / `--refcount_ == 0` they replace, and that
    // branch is what buys immortality: on the objects that matter it turns a
    // contended read-modify-write into a load.
    //
    // The second test, on the value the read-modify-write returns, is what
    // makes marking linearizable. The first test and the RMW are two separate
    // atomic operations, so mark_immortal() can land between them; told
    // nothing, the increment would then move the count off the mark. Reading
    // the pre-RMW value says whether that happened, and the repair puts the
    // count back, so no number of racing losers can accumulate a drift. Both
    // arms elide the RMW entirely once the count is in the band, which is the
    // property benchmarks/rcp_scaling.cpp measures.
    //
    // A compare-exchange loop would be the other way to write this. It is not
    // used because it would replace `lock xadd` with a retry loop on the
    // *mortal* path -- every ordinary object in the process -- to linearize a
    // transition that happens a few dozen times per process, at static init.
    void inc_ref() const noexcept
    {
        if (SYMENGINE_RCP_MORTAL(!detail::refcount_is_immortal(refcount_))) {
            // Post-increment returns the previous value for both
            // `unsigned int` and `std::atomic<unsigned int>`.
            if (!SYMENGINE_RCP_MORTAL(
                    !detail::refcount_is_immortal(refcount_++))) {
                refcount_--;
            }
        }
    }
    //! \return true when this call dropped the last reference.
    bool dec_ref() const noexcept
    {
        if (!SYMENGINE_RCP_MORTAL(!detail::refcount_is_immortal(refcount_))) {
            return false;
        }
        const unsigned int prev = refcount_--;
        if (!SYMENGINE_RCP_MORTAL(!detail::refcount_is_immortal(prev))) {
            // Marking won the race. Undo, and report no last reference: an
            // immortal object has none to drop.
            refcount_++;
            return false;
        }
        return prev == 1;
    }

private:
#else

    // Everything below is private interface for Teuchos
private:
    mutable RCP<T> weak_self_ptr_;

    void set_weak_self_ptr(const RCP<T> &w)
    {
        weak_self_ptr_ = w;
    }

    void set_weak_self_ptr(const RCP<const T> &w) const
    {
        weak_self_ptr_ = rcp_const_cast<T>(w);
    }
#endif // WITH_SYMENGINE_COOPERATIVE_INTRUSIVE_RCP / WITH_SYMENGINE_RCP

#if defined(WITH_SYMENGINE_COOPERATIVE_INTRUSIVE_RCP)                          \
    || defined(WITH_SYMENGINE_RCP)
    template <class T_>
    friend class RCP;
#endif

    template <typename T_, typename... Args>
    friend inline RCP<T_> make_rcp(Args &&...args);
};

template <typename T, typename... Args>
inline RCP<T> make_rcp(Args &&...args)
{
#if defined(WITH_SYMENGINE_COOPERATIVE_INTRUSIVE_RCP)                          \
    || defined(WITH_SYMENGINE_RCP)
    return rcp(new T(std::forward<Args>(args)...));
#else
    RCP<T> p = rcp(new T(std::forward<Args>(args)...));
    p->set_weak_self_ptr(p.create_weak());
    return p;
#endif
}

} // namespace SymEngine

#undef SYMENGINE_RCP_MORTAL

#endif
