cmake_minimum_required(VERSION 3.5)
project(slic)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
        test_slic.cpp
        slic.h slic.cpp)

find_package(OpenCV REQUIRED)

include_directories(${OpenCV_INCLUDE_DIRS})

add_executable(slic ${SOURCE_FILES})

target_link_libraries(slic ${OpenCV_LIBRARIES})