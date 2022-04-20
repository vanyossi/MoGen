
include_directories(
    ${CMAKE_SOURCE_DIR}
)

set(mops_TEST
    mgn_zdt_test.c)

add_executable(zdt_test ${mops_TEST})

target_link_libraries(zdt_test bench_mops)
