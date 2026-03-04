cmake_minimum_required(VERSION 3.10)
project(ContinuousTimeIsing VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(ConTimeIsing ConTimeIsing.cpp)

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(ConTimeIsing PRIVATE -O3 -march=native -ffast-math)
else()
  target_compile_options(ConTimeIsing PRIVATE -O3)
endif()
