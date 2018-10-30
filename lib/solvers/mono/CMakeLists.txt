
include_directories(
    ${CMAKE_SOURCE_DIR}
)

set(mono_SOURCES
    secant.c
    )

add_library(mopmono SHARED ${mono_SOURCES})
target_link_libraries(mopmono mogenmop)

set_target_properties(mopmono PROPERTIES
    VERSION ${MOGEN_VERSION_MAJOR} SOVERSION ${MOGEN_VERSION_MAJOR}
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH_USE_LINK_PATH TRUE
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/${PROJECT_NAME};@loader_path/lib/mogen"
    )
