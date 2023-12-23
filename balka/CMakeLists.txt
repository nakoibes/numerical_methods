cmake_minimum_required(VERSION 3.26)
project(balka)

set(CMAKE_CXX_STANDARD 17)
include_directories(eigen)
add_executable(balka main.cpp)
