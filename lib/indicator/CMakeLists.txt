
add_subdirectory(hv)

include_directories(
    ../math
)

set(IND_SOURCES
    indicadores.c
    )

add_library(fp_indicator STATIC ${IND_SOURCES})
