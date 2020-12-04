cmake_minimum_required(VERSION 3.16)
project(HGS_CARP)

set(CMAKE_CXX_STANDARD 14)

include_directories(Program)

add_executable(HGS_CARP
        Program/Client.cpp
        Program/Client.h
        Program/commandline.cpp
        Program/commandline.h
        Program/Genetic.cpp
        Program/Genetic.h
        Program/Individual.cpp
        Program/Individual.h
        Program/LocalSearch.cpp
        Program/LocalSearch.h
        Program/main.cpp
        Program/Node.cpp
        Program/Node.h
        Program/Params.cpp
        Program/Params.h
        Program/Population.cpp
        Program/Population.h
        Program/Route.cpp
        Program/Route.h
        Program/SeqData.cpp
        Program/SeqData.h
        Program/Vehicle.cpp
        Program/Vehicle.h)
