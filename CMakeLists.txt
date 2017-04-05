cmake_minimum_required(VERSION 3.5)
project(Vortex-HPX CXX)

set(CXX g++-5) # Use g++ version 5 (required to use openMP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -g -fopenmp -march=native -mtune=native -D ORDER=12 -MP -MD")

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY bin/)

include_directories(vrtx)

find_package(HPX)

# build all of this with
# cmake -DHWLOC_ROOT=/usr/local/include/hwloc
#	-DHPX_WITH_MALLOC=jemalloc
#       -DHPX_DIR=~/Documents/CSCS/hpx-build/lib/cmake/HPX


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
# regular non-hpx program						#
#####################################

#set (m4_cmd "/usr/bin/m4 ")
#set (m4_arg "-D ORDER=12 ")
#set (m4_inp "src/kernels.m4")
#set (m4_out "src/kernels.cpp")
#message("Unrolling with M4... ")
#execute_process(
#	COMMAND ${m4_cmd} ${m4_cmd}
#	INPUT_FILE ${m4_inp}
#	OUTPUT_FILE ${m4_out}
#)
#message("Generated ${m4_out} from ${m4_inp}")
#message("Creating executable... ")
#add_executable(vortex src/main.cpp src/kernels.cpp src/tree.cpp)

#####################################
# hpx program												#
#####################################

add_hpx_executable(hpx-vortex
  ESSENTIAL
  SOURCES src/hpx-main.cpp
        src/kernels.cpp
					src/run-tests.cpp
					src/tree.cpp
  COMPONENT_DEPENDENCIES iostreams
)
