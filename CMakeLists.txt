cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

#3d Density kernel estimation program
project(DKE)

find_package(PCL 1.2 REQUIRED)
find_package(Eigen)

find_package(PkgConfig REQUIRED)
pkg_check_modules(JSONCPP jsoncpp)

include_directories(${jsoncpp_INCLUDE_DIRS})
link_libraries(${jsoncpp_LIBRARIES})


set (CMAKE_CXX_FLAGS "-O3 -mtune=native")

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})


set_source_files_properties(lodepng.cpp PROPERTIES COMPILE_FLAGS "-ansi -pedantic")

add_executable (DKE main.cpp lodepng.cpp)
target_link_libraries (DKE ${PCL_LIBRARIES} ${JSONCPP_LIBRARIES})
