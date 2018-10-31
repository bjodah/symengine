#include <pybind11/pybind11.h>



/*

I'm debugging an issue in one my codes where ubsan in clang is flagging undefined behaviour:
```
pykinetgine/tests/test_crn.py::test_CRN /drone/syme-af64b72-rel/lib/cmake/symengine/../../../include/symengine/symengine_rcp.h:147:13: runtime error: member call on address 0x0000014bdf10 which does not point to an object of type 'const SymEngine::Basic'
0x0000014bdf10: note: object is of type 'SymEngine::Mul'
 00 00 00 00  e8 1d 8a ea a7 7f 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  0f 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              vptr for 'SymEngine::Mul'
    #0 0x7fa7ec1fb385 in SymEngine::RCP<SymEngine::Basic const>::~RCP() /drone/syme-af64b72-rel/lib/cmake/symengine/../../../include/symengine/symengine_rcp.h:147:13
    #1 0x7fa7ec411b0e in kinetgine::ReacA::~ReacA() /drone/src/github.com/bjodah/kinetgine/./kinetgine/crn.hpp:36:8
```
My class ``ReacA`` is really simple and has a member

 */
