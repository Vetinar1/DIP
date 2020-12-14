cmake_minimum_required(VERSION 3.17)
project(Master_Project_C_Part)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(Master_Project_C_Part
        CoolCool.h
        CoolManager.h
        main.cpp CoolPoint.h CoolSimplex.h CoolConst.h CoolMLI.h CoolCInterface.cpp CoolCInterface.h CoolCool.cpp CoolManager.cpp CoolPoint.cpp CoolSimplex.cpp)
