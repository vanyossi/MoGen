
include_directories(
    ${CMAKE_SOURCE_DIR}
)

set(mops_TEST
    zdt_test.c)

add_executable(test_zdt ${mops_TEST})

target_link_libraries(test_zdt bench_mops)
