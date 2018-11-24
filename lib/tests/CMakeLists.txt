
include_directories(
    ${CMAKE_SOURCE_DIR}
)

add_executable(poptest pop_test.c)
target_link_libraries(poptest mogenmop moa mogenops)
