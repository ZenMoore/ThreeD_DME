cmake_minimum_required(VERSION 3.14)
project(ThreeD_DME)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)
include_directories(ZstTree)
include_directories(ZstTree/foundation)

add_executable(ThreeD_DME
        ZstTree/foundation/BaseDefine.h
        ZstTree/foundation/BiStates.h
        ZstTree/foundation/Coordinate.h
        ZstTree/foundation/DOUBLE.h
        ZstTree/foundation/MyInterval.h
        ZstTree/foundation/Trr.h
        ZstTree/calculate.cpp
        ZstTree/check.cpp
        ZstTree/header.h
        ZstTree/print.cpp
        ZstTree/read_file.cpp
        ZstTree/utils.cpp
        ZstTree/utils.h
        ZstTree/zst.cpp
        ZstTree/ZstTree.h
        main.cpp
        ThreeD_DME.h)
