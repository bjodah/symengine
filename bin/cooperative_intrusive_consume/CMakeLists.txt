# Standalone consumer: builds a minimal nanobind extension against an *installed*
# SymEngine that was configured with -DSYMENGINE_RCP_BACKEND=cooperative_intrusive.
# This proves both that the install is consumable via find_package(SymEngine) and
# that the cooperative_intrusive backend cooperates with nanobind's intrusive ptr.
cmake_minimum_required(VERSION 3.18)
project(cooperative_intrusive_consume LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(SymEngine CONFIG REQUIRED)

# Locate nanobind via its installed CMake package (pip install nanobind).
find_package(Python 3.8 REQUIRED COMPONENTS Interpreter Development.Module)
execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(min_ext min_ext.cpp)
target_link_libraries(min_ext PRIVATE ${SYMENGINE_LIBRARIES})
