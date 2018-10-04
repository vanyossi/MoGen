
include_directories(
    ${CMAKE_SOURCE_DIR}
)

add_executable(moptest mop_test.c)
target_link_libraries(moptest mogenmop)
