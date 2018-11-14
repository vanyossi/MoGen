
add_subdirectory(mono)

include_directories(
    ${CMAKE_SOURCE_DIR}
)

set(MGF_MOA_SOURCES
    mgf_nsga2.c
    )

add_library(moa ${MGF_MOA_SOURCES})
target_link_libraries(moa mogenmop)

set_target_properties(moa PROPERTIES
    VERSION ${MOGEN_VERSION_MAJOR} SOVERSION ${MOGEN_VERSION_MAJOR}
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH_USE_LINK_PATH TRUE
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/${PROJECT_NAME};@loader_path/lib/mogen"
    )
