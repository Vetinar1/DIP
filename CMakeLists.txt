cmake_minimum_required(VERSION 3.15)
project(Master_Project_C_Part)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(Master_Project_C_Part
        CoolCool.h
        main.cpp CoolPoint.h CoolSimplex.h CoolConst.h CoolMLI.h CoolCool.cpp CoolSimplex.cpp)
