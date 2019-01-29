cmake_minimum_required (VERSION 3.4)

project(sample_sorter)

set(CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

#Adding cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/modules/)

#######################
## Internal Libraries
#######################

include_directories(include)
file(GLOB PROJECT_SOURCES src/*.cpp)

#######################
## External Libraries
#######################

# FFMPEG
find_package(FFMPEG REQUIRED)
include_directories(${FFMPEG_INCLUDE_DIR})
set(LIBS ${LIBS} ${FFMPEG_LIBRARIES})

# FFTW3
find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDES})
set(LIBS ${LIBS} ${FFTW_LIBRARIES})

####################
## Library Creation
####################

add_library(${PROJECT_NAME} SHARED ${PROJECT_SOURCES})
target_link_libraries(${PROJECT_NAME} ${LIBS})
set(LIBS ${LIBS} ${PROJECT_NAME})

#####################################
## Unit Tests
#####################################

# Add the tests
enable_testing()

# from list of files we'll create tests test_name.cpp -> test_name
file(GLOB TEST_SOURCES test/*.cpp)
foreach(_test_file ${TEST_SOURCES})
    get_filename_component(_test_name ${_test_file} NAME_WE)
    add_executable(${_test_name} ${_test_file})
    target_link_libraries(${_test_name} ${LIBS})
    add_test(${_test_name} ${_test_name})
endforeach()
