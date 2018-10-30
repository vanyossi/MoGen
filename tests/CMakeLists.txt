
include_directories(
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/lib/solvers/mono
)

add_executable(moptest mop_test.c)
target_link_libraries(moptest mogenmop)

add_executable(mopmono_test mop_mono.c)
target_link_libraries(mopmono_test mopmono)
