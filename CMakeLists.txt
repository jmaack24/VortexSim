
cmake_minimum_required(VERSION 3.15)

project(VortexSim VERSION 1.0
                  DESCRIPTION "2D Point Vortex Simulator"
                  LANGUAGES CXX)

# set(CMAKE_CXX_STANDARD 11)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_POSITION_INDEPENDENT_CODE True)

# add_executable(vsim

add_subdirectory(src)
