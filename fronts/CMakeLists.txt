

include_directories(
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/lib
    ${CMAKE_SOURCE_DIR}/lib/math
    ${CMAKE_SOURCE_DIR}/lib/operators
    ./decomposition
)

add_executable(cec09UF cec09UF.c)
target_link_libraries(cec09UF mgnmoefa)

#if(APPLE)
#    target_link_libraries(mogenmop ${vecLib_LINKER_LIBS})
#endif()
