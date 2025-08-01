cmake_minimum_required(VERSION 3.29)
project(Cpp_BetaDataStructure VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(Cpp_BetaDataStructure STATIC
        include/BigInt.h
        src/Test.cpp
)
