
include_directories(
    ${PROJECT_SOURCE_DIR}/math
    ${GSL_INCLUDE_DIRS}
)

set(moead_math_SOURCES
    array.c
    gsl_extension.c
    linalg.c
    moea_math.c
)

add_library(moea_math SHARED ${moead_math_SOURCES})
target_link_libraries(moea_math ${GSL_LIBRARY} m)
set_target_properties(moea_math PROPERTIES LINKER_LANGUAGE C
    VERSION ${MOGEN_VERSION_MAJOR} SOVERSION ${MOGEN_VERSION_MAJOR})

install(TARGETS moea_math ${INSTALL_TARGETS_DEFAULT_ARGS})
