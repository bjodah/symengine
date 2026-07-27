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
// This is intentionally tiny -- it exposes Symbol, Add and __str__ only.  The
// full-featured bindings live in the sister `nbsymengine` project.

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <symengine/add.h>
#include <symengine/basic.h>
#include <symengine/number.h>
#include <symengine/symbol.h>
#include <symengine/symengine_rcp.h>

#include <atomic>
#include <string>
#include <type_traits>

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
        // Building an RCP<T> from the T* increments the intrusive count.
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
        return Caster::from_cpp(v.get(), policy, cl);
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

NB_MODULE(min_ext, m)
{
    initialize_intrusive_hooks();

    // The intrusive_ptr policy tells nanobind to hand ownership of a freshly
    // created wrapper to the C++ object via set_self_external(), which flips the
    // cooperative counter into "externally owned" mode.
    nb::class_<Basic>(m, "Basic",
                      nb::intrusive_ptr<Basic>([](Basic *o, PyObject *po) noexcept {
                          o->set_self_external(po);
                      }))
        .def("__str__", [](const Basic &b) { return b.__str__(); });

    // Concrete subtypes registered so the caster returns proper Python types
    // (e.g. isinstance(symbol("x"), Symbol) holds).
    nb::class_<Symbol, Basic>(m, "Symbol");
    nb::class_<Add, Basic>(m, "Add");

    // Factory functions, lowercase to mirror SymEngine's own C++ API
    // (SymEngine::symbol / SymEngine::add) and to avoid clashing with the
    // Symbol / Add type objects above.
    m.def("symbol",
          [](const std::string &name) -> RCP<const Basic> {
              return SymEngine::symbol(name);
          });
    m.def("add",
          [](const RCP<const Basic> &a, const RCP<const Basic> &b) -> RCP<const Basic> {
              return SymEngine::add(a, b);
          });
}
