
include_directories(
    ${CMAKE_SOURCE_DIR}
    ../
    ${CMAKE_SOURCE_DIR}/lib/operators
)

set(mops_TEST
    mgn_zdt_test.c)

add_executable(zdt_test ${mops_TEST})
target_link_libraries(zdt_test mgnmoefa)
