cmake_minimum_required(VERSION 3.19)
project(THADSNBit)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-g")

include_directories(BitSprayer)
include_directories(Epi)


add_executable(testarooni
        BitSprayer/bitspray.cpp
        BitSprayer/bitspray.h
        Epi/main.cpp
        Epi/setu.cpp
        Epi/setu.h
        Epi/stat.cpp
        Epi/filesystem.hpp
        Epi/stat.h )
