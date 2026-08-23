"""Smoke test for the cooperative_intrusive backend driven from Python.

Exercises shared ownership between CPython (nanobind) and SymEngine's
cooperative_intrusive RCP in both of the states a real binding meets:

  * ordinary objects, whose C++ references are folded into their Python
    wrapper's reference count, so every `RCP` copy is a `Py_INCREF`; and

  * the library's *immortal* process-global constants (`zero`, `one`,
    `minus_one`, `pi`, ...), whose references are not counted at all.  The
    hand-off grants their wrapper exactly one permanent Python reference and
    C++ never hands it back, so no later reference traffic touches it.

The second is the state the immortality change exists for -- see
UPSTREAM-NOTES-false-sharing.md -- and this is the only test in the tree that
puts it in front of a real foreign runtime rather than a counting stub.
"""
import gc
import sys

import min_ext

#: use_count() of an object that keeps no count.  "Effectively infinite" is
#: the only honest answer for something that is never released.
UINT_MAX = 2**32 - 1

#: How many simultaneous C++ references to take when measuring whether
#: reference traffic reaches the Python wrapper.
CPP_REFS = 500


def check_ordinary_objects() -> None:
    """Baseline: the cooperation that predates immortality still works."""
    x = min_ext.symbol("x")
    y = min_ext.symbol("y")

    assert isinstance(x, min_ext.Symbol), type(x)
    assert str(x) == "x", str(x)

    s = min_ext.add(x, y)
    assert isinstance(s, min_ext.Add), type(s)
    # SymEngine normalises Add ordering deterministically.
    assert str(s) in ("x + y", "y + x"), str(s)

    # An ordinary object is mortal, is externally owned once Python has seen
    # it, and -- being externally owned -- is never "mine to steal".
    assert not min_ext.is_immortal(x)
    assert min_ext.is_external_owned(x)
    assert not min_ext.is_uniquely_owned(x)

    # Its wrapper owns it: when the last reference goes, nanobind runs the
    # destructor.  This is what must *not* be true of a constant.
    assert min_ext.wrapper_will_destruct(x)

    # And every live C++ reference really is a Python reference on the
    # wrapper.  This is the cost immortality removes.
    base = min_ext.refcount_under_cpp_refs(x, 0)
    under = min_ext.refcount_under_cpp_refs(x, CPP_REFS)
    assert under == base + CPP_REFS, (base, under)

    del x, y, s
    gc.collect()


def check_constants_start_immortal_and_uncounted() -> None:
    """A constant nobody has asked for yet: immortal, unwrapped, uncounted."""
    # `GoldenRatio` is deliberately untouched by every other check here, so
    # this observes the pristine state rather than one this test created.
    immortal, external, count = min_ext.state_by_name("GoldenRatio")
    assert immortal, "constants.cpp should have marked it immortal"
    assert not external, "no binding has attached a wrapper yet"
    assert count == UINT_MAX, count

    # And it is not one lucky constant: constants.cpp marks every entry of
    # DEFINE_CONSTANTS, so a mark that stopped covering one of them is a
    # scaling regression nothing else here would notice.  state_by_name()
    # reads them without handing them to Python, so this stays a read.
    names = min_ext.constant_names()
    assert len(names) >= 10, names
    for name in names:
        immortal, external, count = min_ext.state_by_name(name)
        assert immortal, name
        if not external:
            assert count == UINT_MAX, (name, count)


def check_wrapping_a_constant_keeps_it_immortal() -> None:
    """The hand-off must not undo immortality, and must not take ownership."""
    one = min_ext.constant("one")
    assert isinstance(one, min_ext.Integer), type(one)
    assert str(one) == "1", str(one)

    immortal, external, count = min_ext.state_by_name("one")
    assert immortal, "attaching a wrapper must not clear the mark"
    assert external, "the wrapper is the external owner"
    # `use_count()` reports 0 rather than UINT_MAX here, because external
    # ownership is the more specific answer: the count is the foreign
    # runtime's now.  `is_immortal()` is the predicate that stays true, and
    # the ownership gate stays shut either way.
    assert count == 0, count
    assert not min_ext.is_uniquely_owned(one)

    # Identity reuse: asking again returns the very same Python object.
    assert min_ext.constant("one") is one
    assert min_ext.wrapper_is_self(one, one)

    # Nothing may delete a process-global singleton -- least of all a garbage
    # collector.  nanobind overrides the return-value policy to
    # `take_ownership` for intrusive types, so the caster has to take that
    # ownership back; this is the assertion that it did.
    assert not min_ext.wrapper_will_destruct(one), (
        "the wrapper owns the constant: dropping it would delete SymEngine::one"
    )


def check_cpp_references_never_reach_the_wrapper() -> None:
    """The point of the whole exercise, measured against CPython's counts."""
    one = min_ext.constant("one")

    # 500 simultaneous C++ references to `one`.  On an ordinary externalised
    # object each one is a Py_INCREF on the wrapper (check_ordinary_objects
    # measures exactly that); on an immortal one the reference operations do
    # not touch the state word at all.
    base = min_ext.refcount_under_cpp_refs(one, 0)
    under = min_ext.refcount_under_cpp_refs(one, CPP_REFS)
    assert under == base, (base, under)

    # The same thing through a realistic workload rather than a synthetic
    # loop: `SymEngine::add` puts `zero` in every Add's coefficient slot, so
    # these expressions hold 500 live C++ references to it.
    zero = min_ext.constant("zero")
    before = min_ext.wrapper_refcount(zero)
    assert before > 0, before
    exprs = [
        min_ext.add(min_ext.symbol("x%d" % i), min_ext.symbol("y%d" % i))
        for i in range(CPP_REFS)
    ]
    assert str(exprs[0]) in ("x0 + y0", "y0 + x0"), str(exprs[0])
    assert min_ext.wrapper_refcount(zero) == before, min_ext.wrapper_refcount(zero)

    del exprs
    gc.collect()
    assert min_ext.wrapper_refcount(zero) == before, min_ext.wrapper_refcount(zero)


def check_detach_round_trip() -> None:
    """Shutdown: give the wrappers back, and survive having done so."""
    one = min_ext.constant("one")

    detached = min_ext.detach_constant_wrappers()
    assert detached >= 1, detached

    # Back to immortal-*without*-a-wrapper, not to a count of zero -- so no
    # later release can underflow.
    immortal, external, count = min_ext.state_by_name("one")
    assert immortal and not external, (immortal, external)
    assert count == UINT_MAX, count
    assert min_ext.wrapper_refcount(one) == -1

    # Idempotent: nothing left to detach.
    assert min_ext.detach_constant_wrappers() == 0

    # The object outlived its wrapper, which is the guarantee.  Dropping the
    # last Python reference here is precisely the moment a wrapper that had
    # kept ownership would `delete` SymEngine::one; the next line would then
    # be a use-after-free rather than a "1".
    del one
    gc.collect()
    assert str(min_ext.constant("one")) == "1"

    # That last call re-attached a wrapper; hand it back so the interpreter
    # can finalize with nothing of ours left on the C++ statics.  (The module
    # also registers this with atexit, so a caller that forgets is fine too.)
    min_ext.detach_constant_wrappers()


def main() -> int:
    for check in (
        check_ordinary_objects,
        check_constants_start_immortal_and_uncounted,
        check_wrapping_a_constant_keeps_it_immortal,
        check_cpp_references_never_reach_the_wrapper,
        check_detach_round_trip,
    ):
        check()
        print("  ok:", check.__name__)

    print("cooperative_intrusive python smoke: OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
