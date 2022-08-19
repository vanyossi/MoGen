
include_directories(
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/lib/
    ${CMAKE_SOURCE_DIR}/lib/math
    ${CMAKE_SOURCE_DIR}/lib/solvers
    ${CMAKE_SOURCE_DIR}/lib/operators
)

add_executable(mgnoptest mgn_optest.c )
target_link_libraries(mgnoptest mgnmoefa)

add_executable(gsl_kmeans_test gsl_kmeans_test.c ../gsl_kmeans_new.c)
target_link_libraries(gsl_kmeans_test mgnmoefa)

add_executable(mgn_rbftest mgn_rbf_test.c)
target_link_libraries(mgn_rbftest mgnmoefa)
configure_file(training_set.txt ${CMAKE_CURRENT_BINARY_DIR}/training_set.txt COPYONLY)
configure_file(testing_set.txt ${CMAKE_CURRENT_BINARY_DIR}/testing_set.txt COPYONLY)
