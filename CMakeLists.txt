cmake_minimum_required(VERSION 3.12)
project(losat)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(losat
        gateLevelCkt.cpp
        main.cpp)

add_executable(fsim
        my_fsim.cpp)

add_executable(testckt
        gateLevelCkt.cpp
        main_test_ckt.cpp)
