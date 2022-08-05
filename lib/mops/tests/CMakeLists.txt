
include_directories(
    ${CMAKE_SOURCE_DIR}
)

set(mops_TEST
    mgn_zdt_test.c)

add_executable(zdt_test ${mops_TEST})
