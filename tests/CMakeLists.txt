
include_directories(
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/lib
    ${CMAKE_SOURCE_DIR}/lib/math
    ${CMAKE_SOURCE_DIR}/lib/solvers
    ${CMAKE_SOURCE_DIR}/lib/solvers/decomposition
    ${CMAKE_SOURCE_DIR}/lib/mops
    ${CMAKE_SOURCE_DIR}/lib/operators
)

# mgn tests
add_executable(moeadtest mgn_moeadtest.c)
target_link_libraries(moeadtest mgnmoefa)

add_executable(edtest edtest.c)
target_link_libraries(edtest mgnmoefa)

add_executable(moead-rbftest mgn_moead-rbf_test.c)
target_link_libraries(moead-rbftest mgnmoefa)
