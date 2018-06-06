cmake_minimum_required (VERSION 3.4)
project (SampleSorter)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

#Adding cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/modules/)

#######################
## Internal Libraries
#######################

set(PROJECT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")
set(PROJECT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(PROJECT_SOURCES
  ${PROJECT_SOURCE_DIR}/AudioFile.cpp

  ${PROJECT_SOURCE_DIR}/AbletonSampleFile.cpp
  ${PROJECT_SOURCE_DIR}/SampleFile.cpp

  ${PROJECT_SOURCE_DIR}/AudioSample.cpp
  ${PROJECT_SOURCE_DIR}/TimeDomainProcessing.cpp
  ${PROJECT_SOURCE_DIR}/SpectralProcessing.cpp
  ${PROJECT_SOURCE_DIR}/Octave.cpp
  ${PROJECT_SOURCE_DIR}/Tempo.cpp
  ${PROJECT_SOURCE_DIR}/EqualLoudness.cpp

  ${PROJECT_SOURCE_DIR}/Units.cpp

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

# FFMPEG
find_package(FFMPEG REQUIRED)
include_directories(${FFMPEG_INCLUDE_DIR})
set(LIBS 
  ${LIBS} 
  ${FFMPEG_LIBRARIES} 
  #${FFMPEG_avformat_LIBRARY}
  #${FFMPEG_avcodec_LIBRARY}
)
if(NOT FFMPEG_FOUND)
  message(ERROR "FFMPEG not found!")
endif(NOT FFMPEG_FOUND)

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

#####################################
## Unit Tests
#####################################

macro(MAKE_TEST test_name)
  add_executable(${test_name}
         test/${test_name}.cpp
         )
  target_link_libraries(${test_name} ${LIBS})
  add_test(${test_name} ${test_name})
endmacro()

# MAKE_TEST(make_examples)
