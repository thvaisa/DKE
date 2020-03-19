cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

#3d Density kernel estimation program
project(DKE)


file(DOWNLOAD https://raw.githubusercontent.com/lvandeve/lodepng/master/lodepng.cpp ${CMAKE_CURRENT_SOURCE_DIR}/lodepng.cpp)
file(DOWNLOAD https://raw.githubusercontent.com/lvandeve/lodepng/master/lodepng.h ${CMAKE_CURRENT_SOURCE_DIR}/lodepng.h)



find_package(PkgConfig REQUIRED)
pkg_check_modules(JSONCPP jsoncpp)


find_package(Eigen3 REQUIRED)


set (CMAKE_CXX_FLAGS "-O3 -mtune=native -std=c++11")
#set (CMAKE_CXX_FLAGS "-mtune=native -std=c++11 -pg")

include_directories(${EIGEN3_INCLUDE_DIR})

link_libraries(${JSONCPP_LIBRARIES})


set_source_files_properties(lodepng.cpp PROPERTIES COMPILE_FLAGS "-ansi -pedantic")

add_executable (DKE main.cpp lodepng.cpp)
target_link_libraries (DKE ${EIGEN3_LIBRARIES} ${JSONCPP_LIBRARIES})

