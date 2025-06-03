cmake_minimum_required(VERSION 3.16)
project(SLM)
set(CMAKE_CXX_STANDARD 20)

#add_compile_options(
#  -Wno-unused-variable
#  -Wno-deprecated-declarations
#  -Wno-pragmas
#)

# 加入 Biovoltron 子專案
add_subdirectory(Biovoltron)

## Enable debug symbols in Debug mode
#set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)

## Set flags for Debug build
#set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -march=native")

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/inc
    ${CMAKE_CURRENT_SOURCE_DIR}/Biovoltron/include
)

# Source files
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

# Executable output path
file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# Build target
add_executable(main ${SOURCES})

# Link required libraries
target_link_libraries(main PRIVATE biovoltron)