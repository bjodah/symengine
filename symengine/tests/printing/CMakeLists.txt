project(test_printing)

add_executable(test_printing test_printing.cpp)
target_link_libraries(test_printing symengine catch)
add_test(test_printing ${PROJECT_BINARY_DIR}/test_printing)

add_executable(test_ccode test_ccode.cpp)
target_link_libraries(test_ccode symengine catch)
add_test(test_ccode ${PROJECT_BINARY_DIR}/test_ccode)

add_executable(test_srepr test_srepr.cpp)
target_link_libraries(test_srepr symengine catch)
add_test(test_srepr ${PROJECT_BINARY_DIR}/test_srepr)
