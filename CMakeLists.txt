cmake_minimum_required( VERSION 2.6 )

project( gsSandbox )

# load cmake options
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# find libraries
find_package(gismo REQUIRED)
include_directories( ${GISMO_INCLUDE_DIRS})
link_directories   ( ${GISMO_LIBRARY_DIR} )
include( ${GISMO_USE_FILE} )

message(STATUS "G+SMo shared library found in ${GISMO_LIBRARY_DIR}")


# set build directories
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)


add_subdirectory( src )
