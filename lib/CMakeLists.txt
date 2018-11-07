add_subdirectory(math)
add_subdirectory(operators)
add_subdirectory(tests)

include_directories(
    ${CMAKE_SOURCE_DIR}
)

set(mogen_mop_SOURCES
    mogen_mop.c
    mgf_population.c
    mgf_individual.c
    multi_array.c
    operators/crossover.c
    rand.c
    mogen_moa.h
    mogen_moa.c
    mop_report.c
    time_tools.c
    )

add_library(mogenmop SHARED ${mogen_mop_SOURCES})
target_link_libraries(mogenmop moea_math)

set_target_properties(mogenmop PROPERTIES
    VERSION ${MOGEN_VERSION_MAJOR} SOVERSION ${MOGEN_VERSION_MAJOR}
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH_USE_LINK_PATH TRUE
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/${PROJECT_NAME};@loader_path/lib/mogen"
    )

add_subdirectory(solvers/mono)
