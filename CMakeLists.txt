cmake_minimum_required(VERSION 2.8)

project(SurrogateMOEAD CXX C)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build)
#set(INSTALL_LIB_DIR ${PROJECT_SOURCE_DIR}/build)

set(CMAKE_BUILD_TYPE RELEASE)
add_definitions(-DEIGEN_NO_DEBUG)
#set(CMAKE_BUILD_TYPE DEBUG)

list(APPEND CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
set_property(GLOBAL APPEND PROPERTY
             COMPILER_FLAGS " -Wall -Wextra -O2 -Wno-enum-compare -march=native -mavx")


######### DEPENDENCIES #########

# gtest
add_subdirectory(lib/gtest)
include_directories(lib/gtest)

# hv
add_subdirectory(lib/hv)
include_directories(lib/hv)

# eigen
include_directories(lib/eigen)
#add_subdirectory(lib/eigen)

# cppformat
include_directories(lib/cppformat)
add_subdirectory(lib/cppformat)

# inireader
include_directories(lib/inireader)
add_subdirectory(lib/inireader)



######### EXECUTABLES #########

add_subdirectory(src)

# ELMOEA/D
add_executable(moead_surrogate src/moead_surrogate.cpp)
target_link_libraries(moead_surrogate surrogate hv inireader)

# MOEA/D
add_executable(moead_alg src/moead_alg.cpp)
target_link_libraries(moead_alg surrogate hv)

# test
add_executable(test_all src/test/test_all.cpp)
target_link_libraries(test_all gtest gtest_main hv surrogate)
