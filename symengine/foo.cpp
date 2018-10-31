// export UBSAN_OPTIONS=print_stacktrace=1,halt_on_error=1,external_symbolizer_path=/usr/lib/llvm-6.0/bin/llvm-symbolizer
// CMAKE_PREFIX_PATH=$(pwd)/build-dbg clang++-6.0 -fsanitize=undefined $(cmake --find-package -DNAME=SymEngine -DCOMPILER_ID=CLANG -DLANGUAGE=CXX -DMODE=COMPILE) symengine/foo.cpp $(cmake --find-package -DNAME=SymEngine -DCOMPILER_ID=CLANG -DLANGUAGE=CXX -DMODE=LINK)

#include <vector>
#include <symengine/parser.h>
#include <iostream>

namespace se = SymEngine;
using se_basic = se::RCP<const se::Basic>;

struct Bar {
    std::vector<int> baz;
    Bar(std::vector<int> baz_) : baz(baz_) {}
};

struct Foo : public Bar {
    se_basic expr;
    Foo(std::vector<int> baz, se_basic arg) : Bar(baz), expr(arg) {}
};

int main(int argc, char ** argv){
    std::vector<Foo> foos;
    for (int i=1; i<argc; ++i) {
        foos.push_back(Foo({i}, se::parse(argv[i])));
        std::cout << foos.back().expr->__str__() << std::endl;
    }
}
