

set(HV_VERSION_MAJOR "1")
set(HV_VERSION_MINOR "3")

set(HV_SOURCES
    hv.c
    avl.c
    )

add_library(fpli_hv SHARED ${HV_SOURCES})

target_compile_definitions(fpli_hv PUBLIC
    -DVARIANT=4
    -DARCH="${CMAKE_SYSTEM_PROCESSOR}"
    )

set_target_properties(fpli_hv PROPERTIES
        VERSION ${HV_VERSION_MAJOR} SOVERSION ${HV_VERSION_MINOR}
        BUILD_WITH_INSTALL_RPATH TRUE
        INSTALL_RPATH_USE_LINK_PATH TRUE
#        INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/${PROJECT_NAME};@loader_path/lib/mgn"
        )
