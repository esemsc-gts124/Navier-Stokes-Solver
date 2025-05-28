cmake_minimum_required(VERSION 3.10)

project(NS_Solver
        VERSION 1.0
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)


add_executable(NSSolver ftcs_wave_equation.cpp matrix.cpp)

target_include_directories(NSSolver
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR})

