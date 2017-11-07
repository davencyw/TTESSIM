cmake_minimum_required(VERSION 2.8.3)
project(ttessim)

##___________________________
##
##FLAGS
##
#

SET (RELEASE_FLAGS " -O3 -ftree-vectorize -fopenmp -lboost_program_options -lboost_filesystem -lboost_system -mtune=native -fopt-info-vec-optimized -fstrict-aliasing")

SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread -mavx2 -mavx -msse -march=native")

##___________________________
##
##FILES & FOLDERS
##
#

include_directories(include)
#add_subdirectory(xxx)


##DIR
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/dll)
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)


file (GLOB SIMSRC
	src/main.cc
)

file (GLOB TESTINGSRC
	src/main_testing.cc
)

file (GLOB PDESOLVERSRC
	src/pdesolver.cc
)

file (GLOB TSTORAGEUNITSRC
	src/tstorageunit.cc
)


##___________________________
##
##LIBRARIES & EXECUTABLES
##
#

add_library(pdesolver ${PDESOLVERSRC})
add_library(pdesolver_testing ${PDESOLVERSRC})
target_compile_definitions(pdesolver_testing PRIVATE TESTING)
add_library(tstorageunit ${TSTORAGEUNITSRC})

add_executable (ttessim ${SIMSRC})
target_link_libraries(ttessim pdesolver tstorageunit)

add_executable (ttessim_testing  ${TESTINGSRC})
target_compile_definitions(ttessim_testing PRIVATE TESTING)
target_link_libraries(ttessim_testing pdesolver_testing tstorageunit)
