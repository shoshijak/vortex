cmake_minimum_required(VERSION 3.5)
project(Vortex-HPX CXX)

set(CXX g++-5) # Use g++ version 5 (required to use openMP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -g -fopenmp -march=native -mtune=native -D ORDER=12 -MP -MD")

# option to run the program with OpenMP rather than HPX
option(VORTEX_WITH_OPENMP "Use OpenMP for multithreading" OFF)

# option to print (a lot of) messages during execution
option(VORTEX_WITH_PRINT "Print messages for tree-build" OFF)

if(VORTEX_WITH_OPENMP)
	message("-- building VORTEX with OPENMP")
	add_definitions(-DRUN_WITH_OMP)
endif()

if(VORTEX_WITH_PRINT)
	message("-- Will print messages during tree-build")
	add_definitions(-DPRINT)
endif()

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY bin/)

include_directories(vrtx)

if(NOT VORTEX_WITH_OPENMP)
find_package(HPX)
endif()

# if with HPX: build all of this with
# cmake -DHWLOC_ROOT=/usr/local/include/hwloc
#	-DHPX_WITH_MALLOC=jemalloc
#   -DHPX_DIR=~/Documents/CSCS/hpx-build/lib/cmake/HPX


set(SOURCE_FILES
				src/kernels.cpp
				src/main.cpp
				src/hpx-main.cpp
				src/run-tests.cpp
				src/tree.cpp
)

set(HEADER_FILES
				vrtx/kernels.h
				vrtx/simulation.h
				vrtx/timer.h
				vrtx/tree_prepare.h
				vrtx/tree.h
)

#####################################
# regular non-hpx program
#####################################

if(VORTEX_WITH_OPENMP)
add_executable(vortex src/omp-main.cpp
    src/kernels.cpp
    src/run-tests.cpp
    src/tree.cpp)
endif()

#####################################
# hpx program												#
#####################################

if(NOT VORTEX_WITH_OPENMP)
message("building VORTEX with HPX")
add_hpx_executable(hpx-vortex
  ESSENTIAL
  SOURCES src/hpx-main.cpp
    src/kernels.cpp
    src/run-tests.cpp
    src/tree.cpp
  COMPONENT_DEPENDENCIES iostreams
)
endif()
