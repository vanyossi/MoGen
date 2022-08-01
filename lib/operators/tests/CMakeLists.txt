
include_directories(
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/lib/
    ${CMAKE_SOURCE_DIR}/lib/math
    ${CMAKE_SOURCE_DIR}/lib/solvers
    ${CMAKE_SOURCE_DIR}/lib/operators
)

add_executable(mgnoptest mgn_optest.c )
target_link_libraries(mgnoptest mgnmoefa)

add_executable(gsl_kmeanstest gsl_kmeanstest.c ../gsl_kmeans.c)
target_link_libraries(gsl_kmeanstest mgnmoefa)

add_executable(mgnrbftest mgn_rbf_test.c)
target_link_libraries(mgnrbftest)
