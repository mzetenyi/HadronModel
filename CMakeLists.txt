cmake_minimum_required(VERSION 3.0.0)
project(HadronModel VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED On)
set(CMAKE_CXX_EXTENSIONS Off)

add_subdirectory(lib)
add_subdirectory(exec)