include_directories(
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/lib/
    ${CMAKE_SOURCE_DIR}/lib/math
    ${CMAKE_SOURCE_DIR}/lib/solvers
    ${CMAKE_SOURCE_DIR}/lib/operators
)

add_executable(mgntest_gsl_additional gsl_additional_test.c )
target_link_libraries(mgntest_gsl_additional mgnmoefa)
