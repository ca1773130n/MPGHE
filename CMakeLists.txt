cmake_minimum_required(VERSION 2.8)
project(mpghe)
  
IF(MSVC_IDE)
  set(OpenCV_STATIC OFF)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
  add_definitions(-DUSING_CMAKE=1)
ELSE(MSVC_IDE)
  set(CFLAGS_WARN "-Wall -Wextra -Wno-unused-parameter -Wno-strict-aliasing")
  set(CMAKE_CXX_FLAGS "--std=c++11 -fPIC -O3 ${CFLAGS_WARN} ${CMAKE_CXX_FLAGS}")

  find_package(X11 REQUIRED)
  link_libraries(${X11_LIBRARIES})
  include_directories(${X11_INCLUDE_DIR})
ENDIF(MSVC_IDE)

find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/src/lib)

set(OBJECTS
	src/lib/mpghe.cpp
	src/lib/platform.cpp
	src/lib/profile.cpp
	src/lib/debug.cpp
)

set(DEMO_OBJECTS
	src/example/test.cpp
)

IF(MSVC_IDE)
set(OBJECTS
	getopt_win.c
	${OBJECTS}
)
ENDIF(MSVC_IDE)

list(APPEND "-std=c++11 -ftree-vectorize")

add_library(mpghe ${OBJECTS})

IF(BUILD_DEMO)
add_executable(testMPGHE ${DEMO_OBJECTS})
target_link_libraries(testMPGHE mpghe ${OpenCV_LIBS})
ENDIF(BUILD_DEMO)
