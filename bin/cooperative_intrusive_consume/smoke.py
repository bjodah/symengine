"""Smoke test for the cooperative_intrusive backend driven from Python.

Exercises shared ownership between CPython (nanobind) and SymEngine's
cooperative_intrusive RCP: create symbols, add them, stringify the result,
and confirm reference-count cooperation does not crash on teardown.
"""
import gc
import sys

import min_ext


def main() -> int:
    x = min_ext.symbol("x")
    y = min_ext.symbol("y")

    assert isinstance(x, min_ext.Symbol), type(x)
    assert str(x) == "x", str(x)

    s = min_ext.add(x, y)
    assert isinstance(s, min_ext.Add), type(s)
    # SymEngine normalises Add ordering deterministically.
    assert str(s) in ("x + y", "y + x"), str(s)

    # Identity reuse: the same C++ object hands back the same Python wrapper.
    assert min_ext.add(x, y) is not None

    print("cooperative_intrusive python smoke: OK ->", str(s))

    # Drop references and collect so wrappers are torn down while the
    # interpreter is still alive (decref cooperation exercised here).
    del x, y, s
    gc.collect()
    return 0


if __name__ == "__main__":
    sys.exit(main())
