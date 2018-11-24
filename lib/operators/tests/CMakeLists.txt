
include_directories(
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/lib/solvers
)

add_executable(operator_test mgf_operator_test.c)
target_link_libraries(operator_test mogenops mogenmop moa)
