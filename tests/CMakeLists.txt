
include_directories(
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/lib
    ${CMAKE_SOURCE_DIR}/lib/math
    ${CMAKE_SOURCE_DIR}/lib/solvers
    ${CMAKE_SOURCE_DIR}/lib/3rdparty
    ${CMAKE_SOURCE_DIR}/lib/solvers/decomposition
    ${CMAKE_SOURCE_DIR}/lib/mops
    ${CMAKE_SOURCE_DIR}/lib/operators
)

message(message "pareto fron dir ${FP_DIR}")
if(FP_DIR)
    add_compile_definitions(FP_DIR="${FP_DIR}")
endif()

# mgn tests
add_executable(moeadtest mgn_moeadtest.c)
target_link_libraries(moeadtest mgnmoefa)

add_executable(moead_de_test moead_de_test.c)
target_link_libraries(moead_de_test mgnmoefa fpli_hv fp_indicator)

add_executable(edtest edtest.c)
target_link_libraries(edtest mgnmoefa fpli_hv fp_indicator)

add_executable(moead-rbftest mgn_moead-rbf_test.c)
target_link_libraries(moead-rbftest mgnmoefa fpli_hv fp_indicator)

add_executable(moead_fcrbf_test mgn_moead_fcrbf_test.c)
target_link_libraries(moead_fcrbf_test mgnmoefa fpli_hv fp_indicator)


add_executable(moead_fcrbf_inv_test
        ${fcrbf_inv_SOURCES}
        ${CMAKE_SOURCE_DIR}/lib/solvers/mgn_moead_fcrbf_inv.c
        mgn_moead_fcrbf_inv_test.c)
target_link_libraries(moead_fcrbf_inv_test mgnmoefa fpli_hv fp_indicator)

add_executable(moead_fcrbf_nn_test
    ${CMAKE_SOURCE_DIR}/lib/solvers/mgn_moead_fcrbf_nn.c
    mgn_moead_fcrbf_nn_test.c)
target_link_libraries(moead_fcrbf_nn_test mgnmoefa fpli_hv fp_indicator)

add_executable(moead_de_test_e
    moead_de_test_exp.c)
target_link_libraries(moead_de_test_e mgnmoefa fpli_hv fp_indicator)
