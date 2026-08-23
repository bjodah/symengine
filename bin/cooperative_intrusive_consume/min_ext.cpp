// Minimal, hand-written (no litgen) nanobind extension that drives SymEngine's
// `cooperative_intrusive` RCP backend through nanobind's *intrusive pointer*
// machinery.
//
// Purpose: a self-contained smoke test proving that a foreign runtime (CPython,
// via nanobind) can cooperatively share ownership of SymEngine objects when
// SymEngine is built with -DSYMENGINE_RCP_BACKEND=cooperative_intrusive.  The
// backend itself has no dependency on nanobind; the cooperation is established
// purely through the generic incref/decref hooks installed at module init.
//
// It also covers the one state no in-tree test can reach with a *real*
// runtime: an immortal object that a foreign runtime has attached a wrapper to
// (`detail::cooperative_immortal_wrapped`).  SymEngine's process-global
// constants are immortal -- see constants.cpp and UPSTREAM-NOTES-false-
// sharing.md -- so a binding meets that state the first time it hands
// `SymEngine::one` to Python.  The C++ unit tests exercise it against a fake
// hook pair that only counts calls; here it runs against CPython's own
// reference counts, which is the thing the design is actually about.
//
// This is intentionally tiny -- it exposes Symbol, Add, the constants and
// enough introspection to assert the above.  The full-featured bindings live
// in the sister `nbsymengine` project.

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <symengine/add.h>
#include <symengine/basic.h>
#include <symengine/constants.h>
#include <symengine/integer.h>
#include <symengine/number.h>
#include <symengine/symbol.h>
#include <symengine/symengine_rcp.h>

#include <atomic>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace nb = nanobind;
using namespace SymEngine;

// --------------------------------------------------------------------------
// Minimal type caster: RCP<T> (and RCP<const T>) <-> Python wrapper.
// --------------------------------------------------------------------------
namespace nanobind::detail
{
template <typename T>
struct type_caster<SymEngine::RCP<T>> {
    using U = std::remove_const_t<T>;
    using Caster = make_caster<T>;
    NB_TYPE_CASTER(SymEngine::RCP<T>, Caster::Name)

    bool from_python(handle src, uint8_t flags, cleanup_list *cl) noexcept
    {
        Caster caster;
        if (!caster.from_python(src, flags, cl))
            return false;
        // Building an RCP<T> from the T* increments the intrusive count --
        // unless the object is immortal, where the whole point is that it
        // does not.
        value = SymEngine::RCP<T>(caster.operator U *());
        return true;
    }

    static handle from_cpp(const SymEngine::RCP<T> &v, rv_policy policy,
                           cleanup_list *cl) noexcept
    {
        if (!v.get())
            return none().release();
        // Identity reuse: if this object already owns a Python wrapper, hand
        // back the very same wrapper (this is what makes `x is x` hold).
        if (void *o = v->self_external())
            return handle(reinterpret_cast<PyObject *>(o)).inc_ref();
        handle h = Caster::from_cpp(v.get(), policy, cl);
        // An immortal object outlives the process and is owned by nobody, so
        // its wrapper must never be the thing that deletes it.  Asking for
        // `rv_policy::reference` at the call site does not achieve that:
        // nb_type_put_common() *overrides* the policy to `take_ownership` for
        // every type registered with nb::intrusive_ptr, on the assumption that
        // the object's own reference count decides when the wrapper dies.  For
        // an ordinary cooperative object that assumption holds.  For an
        // immortal one it does not -- nothing is counting -- so the ownership
        // has to be taken back explicitly, right here, where no call site can
        // forget it.
        //
        // Without this, tp_dealloc runs `delete` on a process-global
        // singleton.  It cannot fire while the hand-off's permanent Python
        // reference is outstanding, which is exactly why it is easy to miss:
        // it fires the first time a binding gives that reference back at
        // shutdown (see detach_constant_wrappers below), and then
        // `constants.cpp`'s statics point at freed memory.
        if (h.ptr() != nullptr && v->is_immortal())
            nanobind::inst_set_state(h, /*ready=*/true, /*destruct=*/false);
        return h;
    }
};
} // namespace nanobind::detail

// --------------------------------------------------------------------------
// Install the cooperative hooks: route C++ inc/dec of externally-owned objects
// to CPython's reference counting.
// --------------------------------------------------------------------------
static void initialize_intrusive_hooks()
{
    static std::atomic<bool> s_python_is_dead{false};
    s_python_is_dead.store(false, std::memory_order_relaxed);
    Py_AtExit([]() noexcept {
        s_python_is_dead.store(true, std::memory_order_relaxed);
    });
    SymEngine::cooperative_intrusive_init(
        [](void *o) noexcept {
            if (s_python_is_dead.load(std::memory_order_relaxed))
                return;
            nb::gil_scoped_acquire g;
            Py_INCREF(reinterpret_cast<PyObject *>(o));
        },
        [](void *o) noexcept {
            if (s_python_is_dead.load(std::memory_order_relaxed))
                return;
            nb::gil_scoped_acquire g;
            Py_DECREF(reinterpret_cast<PyObject *>(o));
        });
}

// --------------------------------------------------------------------------
// The library's process-global constants.
//
// constants.cpp calls mark_immortal() on every one of these: their reference
// counts are not kept at all, because `Add` stores `zero` in its coefficient
// slot and `Mul` stores `one`, so a live count on them is a write to a cache
// line every core wants.  Three consequences a binding has to respect, all of
// them visible in the tests below:
//
//   * The wrapper must not own the object -- handled once, in the caster
//     above, rather than by annotating every call site.
//   * The hand-off grants exactly *one* permanent Python reference, not one
//     per outstanding C++ reference, and C++ never hands it back.  So no
//     amount of later C++ reference traffic touches the wrapper, and giving
//     the wrapper back at shutdown is one Py_DECREF rather than the
//     refcount-folding arithmetic an ordinary singleton needs.
//   * mark_immortal() must never be called on an object that is already
//     externally owned; it aborts, because the count has by then been folded
//     into the wrapper's own and cannot be discarded on its behalf.  Nothing
//     here does: the library marks these during static initialization, long
//     before any binding can see them.
//
// Holding this vector for the life of the process is likewise free: its RCP
// destructors run during C++ static destruction and are no-ops on immortal
// objects, so there is nothing to order against `constants.cpp`'s own
// teardown.
// --------------------------------------------------------------------------
static const std::vector<std::pair<const char *, RCP<const Basic>>> &constants()
{
    static const std::vector<std::pair<const char *, RCP<const Basic>>> v = {
        {"zero", zero},
        {"one", one},
        {"minus_one", minus_one},
        {"two", two},
        {"I", I},
        {"pi", pi},
        {"E", E},
        {"EulerGamma", EulerGamma},
        {"Catalan", Catalan},
        {"GoldenRatio", GoldenRatio},
    };
    return v;
}

static RCP<const Basic> constant_by_name(const std::string &name)
{
    for (const auto &c : constants()) {
        if (name == c.first)
            return c.second;
    }
    throw nb::key_error(name.c_str());
}

// Give the constants' wrappers back, so the interpreter can free them before
// nanobind's leak checker looks.
//
// For an ordinary singleton this is delicate arithmetic: set_self_external()
// folded every outstanding C++ reference into the wrapper's Python refcount,
// and untangling it means reading Py_REFCNT and compensating with inc_ref()
// calls so that neither side drops the last reference at the wrong moment.
// For an immortal one there is nothing to untangle.  The grant was exactly one
// reference, so giving it back is exactly one Py_DECREF; detach_external()
// returns the object to immortal-*without*-a-wrapper rather than to a count of
// zero, so no later release can underflow; and the wrapper does not own the
// object, so its deallocation frees nothing.  Idempotent: a second call finds
// the constants C++-owned again and detaches nothing.
static int detach_constant_wrappers() noexcept
{
    int detached = 0;
    for (const auto &c : constants()) {
        if (void *o = c.second->detach_external()) {
            Py_DECREF(reinterpret_cast<PyObject *>(o));
            ++detached;
        }
    }
    return detached;
}

// --------------------------------------------------------------------------
// Introspection, so the Python side can assert the cooperation rather than
// merely fail to crash.
// --------------------------------------------------------------------------
static PyObject *wrapper_of(const RCP<const Basic> &b)
{
    return reinterpret_cast<PyObject *>(b->self_external());
}

NB_MODULE(min_ext, m)
{
    initialize_intrusive_hooks();

    // The intrusive_ptr policy tells nanobind to hand ownership of a freshly
    // created wrapper to the C++ object via set_self_external(), which flips
    // the cooperative counter into "externally owned" mode.
    nb::class_<Basic>(
        m, "Basic",
        nb::intrusive_ptr<Basic>(
            [](Basic *o, PyObject *po) noexcept { o->set_self_external(po); }))
        .def("__str__", [](const Basic &b) { return b.__str__(); });

    // Concrete subtypes registered so the caster returns proper Python types
    // (e.g. isinstance(symbol("x"), Symbol) holds).
    nb::class_<Symbol, Basic>(m, "Symbol");
    nb::class_<Add, Basic>(m, "Add");
    nb::class_<Number, Basic>(m, "Number");
    nb::class_<Integer, Number>(m, "Integer");
    nb::class_<Constant, Basic>(m, "Constant");

    // Factory functions, lowercase to mirror SymEngine's own C++ API
    // (SymEngine::symbol / SymEngine::add) and to avoid clashing with the
    // Symbol / Add type objects above.
    m.def("symbol", [](const std::string &name) -> RCP<const Basic> {
        return SymEngine::symbol(name);
    });
    m.def("add",
          [](const RCP<const Basic> &a, const RCP<const Basic> &b)
              -> RCP<const Basic> { return SymEngine::add(a, b); });

    // The immortal process-global constants.  No rv_policy annotation: the
    // caster takes ownership of an immortal object back on its own, which is
    // the property being tested -- an annotation here would be ignored.
    m.def("constant", &constant_by_name, nb::arg("name"));
    m.def("constant_names", []() {
        nb::list names;
        for (const auto &c : constants())
            names.append(nb::str(c.first));
        return names;
    });

    // --- introspection -----------------------------------------------------

    //! (is_immortal, is_external_owned, use_count) for a constant, read
    //! *without* handing it to Python -- the only way to observe a constant
    //! that no wrapper has been attached to yet, since asking for it in any
    //! other way is what attaches one.
    m.def("state_by_name", [](const std::string &name) {
        RCP<const Basic> c = constant_by_name(name);
        return nb::make_tuple(c->is_immortal(), c->is_external_owned(),
                              c->use_count());
    });

    m.def("is_immortal",
          [](const RCP<const Basic> &b) { return b->is_immortal(); });
    m.def("is_external_owned",
          [](const RCP<const Basic> &b) { return b->is_external_owned(); });
    m.def("use_count",
          [](const RCP<const Basic> &b) { return b->use_count(); });
    m.def("is_uniquely_owned",
          [](const RCP<const Basic> &b) { return b->is_uniquely_owned(); });

    //! Whether the object's Python wrapper is `obj` itself -- i.e. whether the
    //! identity-reuse path in the caster above will fire for it.
    m.def("wrapper_is_self", [](const RCP<const Basic> &b, nb::handle obj) {
        return wrapper_of(b) == obj.ptr();
    });

    //! Py_REFCNT of the object's Python wrapper, or -1 if it has none.
    m.def("wrapper_refcount", [](const RCP<const Basic> &b) -> Py_ssize_t {
        PyObject *o = wrapper_of(b);
        return o ? Py_REFCNT(o) : -1;
    });

    //! Whether nanobind would run the C++ destructor when this wrapper dies.
    //! False for an immortal object: nothing may delete a process-global
    //! singleton, least of all a garbage collector.
    m.def("wrapper_will_destruct",
          [](nb::handle obj) { return nb::inst_state(obj).second; });

    //! Hold `n` extra C++ references to `b`, and report Py_REFCNT of its
    //! wrapper *while they are live*.  This is the measurement the whole
    //! immortality change exists for: on an ordinary externalised object every
    //! C++ reference is a Py_INCREF on the wrapper, so the answer grows by
    //! `n`; on an immortal one the reference operations do not touch the state
    //! word at all, so the answer does not move.
    m.def(
        "refcount_under_cpp_refs",
        [](const RCP<const Basic> &b, int n) -> Py_ssize_t {
            std::vector<RCP<const Basic>> hold;
            hold.reserve(static_cast<std::size_t>(n));
            for (int i = 0; i < n; ++i)
                hold.push_back(b);
            PyObject *o = wrapper_of(b);
            return o ? Py_REFCNT(o) : -1;
        },
        nb::arg("b"), nb::arg("n"));

    // --- teardown ----------------------------------------------------------
    m.def("detach_constant_wrappers", &detach_constant_wrappers);

    // A real binding registers this rather than relying on its caller; the
    // smoke test also calls it explicitly, which is safe because it is
    // idempotent.
    nb::module_::import_("atexit").attr("register")(
        m.attr("detach_constant_wrappers"));
}
