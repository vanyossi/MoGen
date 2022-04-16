
include_directories(
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/lib/
    ${CMAKE_SOURCE_DIR}/lib/math
    ${CMAKE_SOURCE_DIR}/lib/solvers
    ${CMAKE_SOURCE_DIR}/lib/operators
)

add_executable(operator_test mgf_operator_test.c)
target_link_libraries(operator_test mogenops mogenmop moa)

add_executable(mgnoptest mgn_optest.c )
target_link_libraries(mgnoptest ${GSL_LIBRARIES} ${BLAS_LIBRARIES})
if(APPLE)
    target_link_libraries(mgnoptest ${vecLib_LINKER_LIBS})
endif()
