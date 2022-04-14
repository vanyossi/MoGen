
include_directories(
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/lib/solvers
)

add_executable(houseoptim houseoptim.c)
target_link_libraries(houseoptim moa mogenmop mogenops bench_mops)
