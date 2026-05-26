cmake_minimum_required(VERSION 3.23)

project(2D_FEM)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(
  -Wall
  -Wextra
  -Wpedantic
)

file(GLOB_RECURSE FEM_SOURCES
    src/*.cpp
)
list(APPEND FEM_SOURCES
    problems/poisson.cpp
)

add_library(fem ${FEM_SOURCES})

target_include_directories(fem PUBLIC include)

add_executable(poisson
    problems/bpx_poisson.cpp
)

target_link_libraries(poisson PRIVATE fem)

add_executable(multigrid_poisson
    problems/multigrid_poisson.cpp
)

target_link_libraries(multigrid_poisson PRIVATE fem)
