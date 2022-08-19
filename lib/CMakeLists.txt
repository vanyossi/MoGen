add_subdirectory(math)
add_subdirectory(operators)
add_subdirectory(mops)
add_subdirectory(solvers)

add_subdirectory(tests)

include_directories(
    mops
    solvers
    solvers/decomposition
    operators
    math
    3rdparty
    ${CMAKE_CURRENT_BINARY_DIR}
)

set(mogen_mop_SOURCES
    time_tools.c)

add_library(mogenmop SHARED ${mogen_mop_SOURCES})
target_link_libraries(mogenmop moea_math)

if(APPLE)
    target_link_libraries(mogenmop ${vecLib_LINKER_LIBS})
endif()

set_target_properties(mogenmop PROPERTIES
    VERSION ${MOGEN_VERSION_MAJOR} SOVERSION ${MOGEN_VERSION_MAJOR}
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH_USE_LINK_PATH TRUE
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/${PROJECT_NAME};@loader_path/lib/mogen"
    )

install(TARGETS mogenmop ${INSTALL_TARGETS_DEFAULT_ARGS})

if(GNUPLOT_FOUND)
    configure_file(gnuplotbin.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/gnuplotbin.h)
endif()

set(MGN_LIB_SOURCES
#    math/gsl_vector_additional.c
#    solvers/decomposition/mgn_scalarization.c
    solvers/decomposition/mgn_weights.c
    solvers/mgn_moead.c
    solvers/mgn_de.c
    solvers/mgn_moead-rbf.c
    mgn_io.c
    mgn_io_datalist.c
    mgn_moa.c
    mgn_initializer.c
    mgn_random.c
    mgn_pop_proto.c
    mgn_poplist.c
    mgn_pop_matrix.c
    mgn_pop_helper.c
    individual.c
    population.c
#    mgn_types.c
    mgn_mop.c
    mgn_gnuplot.c
    mgn_pareto.c
#    mgn_pop_init.c
#    individual.c
#    population.c
    )


add_library(mgnmoefa SHARED ${MGN_LIB_SOURCES})

# ./moeadtest  26,10s user 0,03s system 99% cpu 26,215 total
target_link_libraries(mgnmoefa ${GSL_LIBRARY} ${BLAS_LIBRARIES} mogenops moea_math)
# ./moeadtest  58,79s user 0,52s system 81% cpu 1:12,42 total
# add if APPLE
# ./moeadtest  35,92s user 0,18s system 94% cpu 38,120 total
#target_link_libraries(mgnmoefa ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES})

if(APPLE)
#    message(status "veclib ${vecLib_LINKER_LIBS} -F${BLAS_LIBRARIES}/Versions/Current/Frameworks/vecLib.framework")
#    target_link_libraries(mgnmoefa "-F/System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework")
    target_link_libraries(mgnmoefa "-F${BLAS_LIBRARIES}/Versions/Current/Frameworks/vecLib.framework")
endif()

set_target_properties(mgnmoefa PROPERTIES
    VERSION ${MOGEN_VERSION_MAJOR} SOVERSION ${MOGEN_VERSION_MAJOR}
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH_USE_LINK_PATH TRUE
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/${PROJECT_NAME};@loader_path/lib/mgn"
    )
