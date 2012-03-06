# The name of our project is "mip"

cmake_minimum_required (VERSION 2.6)
project (fl)

SET(CPACK_PACKAGE_VERSION_MAJOR "0")
SET(CPACK_PACKAGE_VERSION_MINOR "0")
SET(CPACK_PACKAGE_VERSION_PATCH "1")
SET(CPACK_PACKAGE_NAME "vrp")
SET(CPACK_SYSTEM_NAME "i386")

set (CMAKE_VERBOSE_MAKEFILE OFF)

set ( HEADERS
/usr/include
/usr/local/include
/usr/include/c++/4.4.3
)

set ( SRC
main
)

FIND_PACKAGE( OpenCV REQUIRED )

add_executable (fl ${SRC})
target_link_libraries (fl ${OpenCV_LIBS} m pthread)
