cmake_minimum_required (VERSION 3.4)
project (SampleSorter)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no c++11 support.")
endif()

#Adding cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/modules/)

#######################
## Internal Libraries
#######################

set(PROJECT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")
set(PROJECT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(PROJECT_SOURCES
  ${PROJECT_SOURCE_DIR}/AbletonSampleFile.cpp
  ${PROJECT_SOURCE_DIR}/SampleFile.cpp

  ${PROJECT_SOURCE_DIR}/AudioSample.cpp
  ${PROJECT_SOURCE_DIR}/TimeDomainProcessing.cpp
  ${PROJECT_SOURCE_DIR}/SpectralProcessing.cpp
  ${PROJECT_SOURCE_DIR}/Octave.cpp
  ${PROJECT_SOURCE_DIR}/Tempo.cpp
  ${PROJECT_SOURCE_DIR}/EqualLoudness.cpp
  ${PROJECT_SOURCE_DIR}/MusicTheory.cpp

  ${PROJECT_SOURCE_DIR}/PythonExternals.cpp

  ${PROJECT_SOURCE_DIR}/Plotting.cpp
)

include_directories("${PROJECT_INCLUDE_DIR}")

#######################
## External Libraries
#######################

# Boost
find_package(Boost COMPONENTS iostreams filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})
set(LIBS ${LIBS} ${Boost_LIBRARIES})
if(NOT Boost_FOUND)
  message(ERROR " Boost not found!")
endif(NOT Boost_FOUND)

# TinyXML2
find_package(TINYXML2 REQUIRED)
include_directories(${TINYXML2_INCLUDE_DIR})
set(LIBS ${LIBS} ${TINYXML2_LIBRARIES})
if(NOT TINYXML2_FOUND)
  message(ERROR " TinyXML2 not found!")
endif(NOT TINYXML2_FOUND)

# LibSndFile
find_package(LIBSNDFILE REQUIRED)
include_directories(${LIBSNDFILE_INCLUDE_DIR})
set(LIBS ${LIBS} ${LIBSNDFILE_LIBRARY})
if(NOT LIBSNDFILE_FOUND)
  message(ERROR " LIBSNDFILE not found!")
endif(NOT LIBSNDFILE_FOUND)

# FFTW3
find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDES})
set(LIBS ${LIBS} ${FFTW_LIBRARIES})
if(NOT FFTW_FOUND)
  message(ERROR " FFTW not found!")
endif(NOT FFTW_FOUND)


####################
## Library Creation
####################

add_library(${PROJECT_NAME} SHARED ${PROJECT_SOURCES})
target_link_libraries(${PROJECT_NAME} ${LIBS})
set(LIBS ${LIBS} ${PROJECT_NAME})

##################
## Google Test
#################

# Google Test
# For testing
include(ExternalProject)
ExternalProject_Add(
  gtest
  GIT_REPOSITORY https://github.com/google/googletest.git
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gtest
  TIMEOUT 10
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")

enable_testing()

ExternalProject_Get_Property(gtest source_dir binary_dir)
include_directories("${source_dir}/googletest/include")
set(LIBS ${LIBS} "${binary_dir}/googlemock/gtest/libgtest.a")

#####################################
## Unit Tests with CTest and GTest
#####################################

macro(MAKE_TEST test_name)
  add_executable(${test_name}
         test/${test_name}.cpp
         )
  target_link_libraries(${test_name} ${LIBS})
  add_test(${test_name} ${test_name})
endmacro()

MAKE_TEST(basic_test)
