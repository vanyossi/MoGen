
add_subdirectory(mono)

include_directories(
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/lib
    ${CMAKE_SOURCE_DIR}/lib/math
)

set(MGF_MOA_SOURCES
    mgf_nsga2.c
    mgf_moead.c
    decomposition/mgf_weights.c
    decomposition/mgf_scalarization.c
    mgf_pareto.c
)

add_library(moa ${MGF_MOA_SOURCES})
target_link_libraries(moa mogenmop)

set_target_properties(moa PROPERTIES
    VERSION ${MOGEN_VERSION_MAJOR} SOVERSION ${MOGEN_VERSION_MAJOR}
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH_USE_LINK_PATH TRUE
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/${PROJECT_NAME};@loader_path/lib/mogen"
    )

add_executable(testweight decomposition/tests/testweight.c)
target_link_libraries(testweight ${GSL_LIBRARIES} ${BLAS_LIBRARIES})
if(APPLE)
    target_link_libraries(testweight ${vecLib_LINKER_LIBS})
endif()

add_executable(scalartest decomposition/tests/mgn_scalar_test.c)
target_link_libraries(scalartest ${GSL_LIBRARIES} ${BLAS_LIBRARIES})
if(APPLE)
    target_link_libraries(scalartest ${vecLib_LINKER_LIBS})
endif()

#set(MGN_MOA_SOURCES
#    mgn_moead.h
#)
#add_library(mgnmoa ${MGN_MOA_SOURCES})
#set_target_properties(mgnmoa PROPERTIES LINKER_LANGUAGE C)
