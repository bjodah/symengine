#include <symengine/symengine_rcp.h>

#if SYMENGINE_RCP_KIND == SYMENGINE_RCP_KIND_TEUCHOS
#include <symengine/utilities/teuchos/Teuchos_RCP.hpp>
#elif SYMENGINE_RCP_KIND == SYMENGINE_RCP_KIND_COOPERATIVE
static void (*object_inc_ref_py)(PyObject *) noexcept = nullptr;
static void (*object_dec_ref_py)(PyObject *) noexcept = nullptr;
#endif

namespace SymEngine
{

#if SYMENGINE_RCP_KIND == SYMENGINE_RCP_KIND_EXCLUSIVE

void print_stack_on_segfault()
{
#ifdef WITH_SYMENGINE_TEUCHOS
    Teuchos::print_stack_on_segfault();
#endif
}
#elif SYMENGINE_RCP_KIND == SYMENGINE_RCP_KIND_TEUCHOS
// nothing todo
#elif SYMENGINE_RCP_KIND == SYMENGINE_RCP_KIND_COOPERATIVE
/* The code in this block is adapted from code in the tests directory of Wenzel
   Jakobs BSD-licensed nanobind. [TODO] Copy LICENSE file of nanobind into tree of symengine. */
void Ptr::inc_ref() const noexcept {
    uintptr_t value = this->count_or_pointer.load(std::memory_order_relaxed);

    while (true) {
        if (value & 1) {
            if (!this->count_or_pointer.compare_exchange_weak(value,
                                               value + 2,
                                               std::memory_order_relaxed,
                                               std::memory_order_relaxed))
                continue;
        } else {
            object_inc_ref_py((PyObject *) value);
        }

        break;
    }
}
void Ptr::dec_ref() const noexcept {
    uintptr_t value = this->count_or_pointer.load(std::memory_order_relaxed);

    while (true) {
        if (value & 1) {
            if (value == 1) {
                fprintf(stderr, "Ptr::dec_ref(%p): reference count underflow!", this);
                assert(false);
            } else if (value == 3) {
                delete this;
            } else {
                if (!this->count_or_pointer.compare_exchange_weak(value,
                                                                  value - 2,
                                                                  std::memory_order_relaxed,
                                                                  std::memory_order_relaxed))
                {
                    continue;
                }
            }
        } else {
            object_dec_ref_py((PyObject *) value);
        }
        break;
    }
}

void Ptr::set_self_py(PyObject *o) noexcept {
    uintptr_t value = this->count_or_pointer.load(std::memory_order_relaxed);
    if (value & 1) {
        value >>= 1;
        for (uintptr_t i = 0; i < value; ++i) {
            object_inc_ref_py(o);
        }
        this->count_or_pointer.store((uintptr_t) o);
    } else {
        fprintf(stderr, "Ptr::set_self_py(%p): a Python object was already present!", this);
        assert(false);
    }
}

PyObject *Ptr::self_py() const noexcept {
    uintptr_t value = this->count_or_pointer.load(std::memory_order_relaxed);
    if (value & 1) {
        return nullptr;
    } else {
        return (PyObject *) value;
    }
}

void object_init_py(void (*object_inc_ref_py_)(PyObject *) noexcept,
                    void (*object_dec_ref_py_)(PyObject *) noexcept) {
    object_inc_ref_py = object_inc_ref_py_;
    object_dec_ref_py = object_dec_ref_py_;
}

#else
#error "Unkown value for SYMENGINE_RCP_KIND"
#endif

} // namespace SymEngine
