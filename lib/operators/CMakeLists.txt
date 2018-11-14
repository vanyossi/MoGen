
include_directories(
    ${CMAKE_SOURCE_DIR}
)

set(MOGEN_OP_SOURCES
    crossover.c
    mgf_rank.c
    mgf_rank.h
    mgf_operators.c
    mgf_operators.h
    mgf_mutation.c
    mgf_mutation.h
    mgf_dominance.c
    )

add_library(mogenops SHARED ${MOGEN_OP_SOURCES})
target_link_libraries(mogenops mogenmop)

set_target_properties(mogenops PROPERTIES
    VERSION ${MOGEN_VERSION_MAJOR} SOVERSION ${MOGEN_VERSION_MAJOR}
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH_USE_LINK_PATH TRUE
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/${PROJECT_NAME};@loader_path/lib/mogen"
    )
