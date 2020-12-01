cmake_minimum_required(VERSION 3.17)
project(Master_Project_C_Part)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(Master_Project_C_Part
        cool.h
        CoolManager.h
        main.cpp)
