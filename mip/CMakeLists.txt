# The name of our project is "mip"

cmake_minimum_required (VERSION 2.6)
project (mip)

SET(CPACK_PACKAGE_VERSION_MAJOR "0")
SET(CPACK_PACKAGE_VERSION_MINOR "0")
SET(CPACK_PACKAGE_VERSION_PATCH "1")
SET(CPACK_PACKAGE_NAME "vrp")
SET(CPACK_SYSTEM_NAME "i386")

set (CMAKE_VERBOSE_MAKEFILE OFF)

find_package (OpenMP) 
if(OPENMP_FOUND)
 set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set ( HEADERS
/usr/include
/usr/local/include
/usr/include/c++/4.4.3
)

set ( SRC
main
)

include_directories (${HEADERS} $ENV{GUROBI_HOME}/include)
link_directories ($ENV{GUROBI_HOME}/lib)
add_executable (mip ${SRC})
target_link_libraries (mip m pthread gurobi_c++ gurobi46 gurobi_g++4.1 gurobi_g++4.2)
