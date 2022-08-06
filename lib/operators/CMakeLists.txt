
add_subdirectory(tests)

include_directories(
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/lib/
    ${CMAKE_SOURCE_DIR}/lib/math
)

set(MOGEN_OP_SOURCES
    mgn_vector_distance.c
    )

add_library(mogenops STATIC ${MOGEN_OP_SOURCES})

#set_target_properties(mogenops PROPERTIES
#    VERSION ${MOGEN_VERSION_MAJOR} SOVERSION ${MOGEN_VERSION_MAJOR}
#    BUILD_WITH_INSTALL_RPATH TRUE
#    INSTALL_RPATH_USE_LINK_PATH TRUE
#    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/${PROJECT_NAME};@loader_path/lib/mogen"
#    )
