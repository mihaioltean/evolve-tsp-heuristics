cmake_minimum_required(VERSION 3.3)
project(mpi)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )

set(SOURCE_FILES evolve_tsp_mpi_threads.cpp)
add_executable(mpi ${SOURCE_FILES})


# Require MPI for this project:
find_package(MPI REQUIRED)

set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

#SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -pg -fopenmp -DPROFILE ")
#SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -fopenmp -DNDEBUG")


include_directories(MPI_INCLUDE_PATH)
target_link_libraries(mpi ${MPI_LIBRARIES})