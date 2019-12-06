cmake_minimum_required(VERSION 3.10)
set(CMAKE_C_COMPILER "/opt/intel/compilers_and_libraries_2019.5.281/linux/bin/intel64/icc")
set(CMAKE_CXX_COMPILER "/opt/intel/compilers_and_libraries_2019.5.281/linux/bin/intel64/icpc")

project(test_project)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -g -lmkl_rt")

include_directories(./include/)
include_directories(./include/umfpack)
add_definitions(-DDEBUG)
add_subdirectory(src)

# add_executable(test_exec files/Laplace_1D.cpp)
# add_executable(test_exec files/CD1D_Galerkin.cpp)
# add_executable(test_exec files/Advection_1D_DG.cpp)
add_executable(test_exec files/Advection_1D_DG_Lagrange_Weak.cpp)


target_link_libraries(test_exec _sources /home/suganth/FEM/fem-assignment/lib/libumfpack_only.a)

set(CMAKE_EXE_LINKER_FLAGS "-ldl -lpthread")
