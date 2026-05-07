cmake_minimum_requirement(VERSION 3.23)

project(2D_FEM)

set(CMAKE_CXX_STANDART 20)
set(CMAKE_CXX_STANDART_REQUIRED ON)

add_compile_options(
  -Wall
  -Wextra
  -Wpedantic
)

file(GLOB_RECURSE FEM_SOURCES
    src/*.cpp
)

add_library(fem ${FEM_SOURCES})

target_include_directories(fem PUBLIC include)

add_executable(poisson
    apps/poisson.cpp
)

