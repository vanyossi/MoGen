
include_directories(
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_INSTALL_PREFIX}/include
)

add_executable(poptest pop_test.c)
target_link_libraries(poptest mogenmop moa mogenops)

add_executable(populationtest population_test.c)
target_link_libraries(populationtest mgnmoefa)

add_executable(initializertest initializertest.c ../mgn_initializer.c)
target_link_libraries(initializertest mgnmoefa)

