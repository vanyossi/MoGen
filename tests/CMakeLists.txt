
include_directories(
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/lib
    ${CMAKE_SOURCE_DIR}/lib/math
    ${CMAKE_SOURCE_DIR}/lib/solvers
)

add_executable(moptest mop_test.c)
target_link_libraries(moptest mogenmop)

add_executable(mopmono_test mop_mono.c)
target_link_libraries(mopmono_test mopmono)

add_executable(nsga2_test nsga2_test.c)
target_link_libraries(nsga2_test moa mogenops bench_mops)

add_executable(moead_test moead_test.c)
target_link_libraries(moead_test moa mogenops bench_mops)

add_executable(moeadtest mgn_moeadtest.c)
target_link_libraries(moeadtest mgnmoefa)
