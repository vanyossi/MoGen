
include_directories(
    ${CMAKE_SOURCE_DIR}
)

set(mogen_mop_SOURCES
    mogen_mop.c
    mogen_population.c
    multi_array.c
    crossover.c
    )

add_library(mogenmop SHARED ${mogen_mop_SOURCES})
set_target_properties(mogenmop PROPERTIES
    VERSION ${MOGEN_VERSION_MAJOR} SOVERSION ${MOGEN_VERSION_MAJOR}
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH_USE_LINK_PATH TRUE
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/${PROJECT_NAME};@loader_path/lib/mogen"
    )
