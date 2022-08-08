
add_subdirectory(tests)

include_directories(
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/lib/operators
)

#set(MOP_SOURCES
#    zdt.c)
#
#add_library(bench_mops SHARED ${MOP_SOURCES})
#target_link_libraries(bench_mops mogenmop)
#
#set_target_properties(bench_mops PROPERTIES
#    VERSION ${MOGEN_VERSION_MAJOR} SOVERSION ${MOGEN_VERSION_MAJOR}
#    BUILD_WITH_INSTALL_RPATH TRUE
#    INSTALL_RPATH_USE_LINK_PATH TRUE
#    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/${PROJECT_NAME};@loader_path/lib/mogen"
#    )
