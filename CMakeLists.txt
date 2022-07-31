
message(STATUS "Using CMake version: ${CMAKE_VERSION}")
cmake_minimum_required(VERSION 3.12)
project(mogen C)

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_C_STANDARD 99)

set(MOGEN_VERSION_STRING "0.2.0")
set(MOGEN_VERSION_MAJOR 0)
set(MOGEN_VERSION_MINOR 2)
set(MOGEN_VERSION_FIX 0)
set(MOGEN_YEAR 2018)

set(MOGEN_VERSION ${MOGEN_VERSION_STRING})

message(STATUS "mogen version - ${MOGEN_VERSION}")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
set(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_RELEASE "-O3")
set(CMAKE_C_FLAGS_RELWITHDEBINFO  "-O3")

if("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Weverything -Wno-unused-function")
endif()
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wno-padded -Wno-unused-function")


if(APPLE)
    set(CMAKE_SKIP_BUILD_PATH FALSE)
    set(CMAKE_MACOSX_RPATH TRUE)
#    set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/${PROJECT_NAME}")
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

    list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
    if("${isSystemDir}" STREQUAL "-1")
        set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/${PROJECT_NAME}")
    endif("${isSystemDir}" STREQUAL "-1")
endif()

## Modules
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

## Find GSL
find_package(GSL REQUIRED)
find_package(BLAS)
find_package(Gnuplot)

if(GNUPLOT_FOUND)
    if(APPLE)
        string(REPLACE "gnuplot.app" "Gnuplot.app" GNUPLOT_EXECUTABLE ${GNUPLOT_EXECUTABLE} )
        string(REPLACE "MacOS/gnuplot.sh" "Resources/bin/gnuplot-run.sh -p" GNUPLOT_EXECUTABLE ${GNUPLOT_EXECUTABLE} )
    endif()
    message("Gnuplot found in ${GNUPLOT_EXECUTABLE}")
endif()

if(APPLE)
    find_package(vecLib)
    include_directories(
        ${vecLib_INCLUDE_DIR}
    )
endif()

# if(NOT GSL_ROOT_DIR)
#     message(WARNING "GSL_ROOT_DIR not set, searching GSL on system")
# endif()
# FIND_PACKAGE(GSL REQUIRED)
# message(${GSL_INCLUDE_DIRS})


## install logic
set(INSTALL_TARGETS_DEFAULT_ARGS
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/${PROJECT_NAME}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/include/${PROJECT_NAME}
    )
## Include directories and subidrectories logic
include_directories(
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/lib
    ${CMAKE_INSTALL_PREFIX}/include
)


add_subdirectory(lib)
add_subdirectory(tests)
add_subdirectory(examples)
