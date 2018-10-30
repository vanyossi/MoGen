
include_directories(
    ${CMAKE_SOURCE_DIR}
)

add_executable(moeadzdt1 moead_zdt1.c)
target_link_libraries(moeadzdt1 mogenmop mogenops)
