cmake_minimum_required(VERSION 3.15)

# verbosity
set(CMAKE_VERBOSE_MAKEFILE ON)

# set the project name
project(excogito C)

# to compile libraries such as optimize with openMP. Avoid WARNINGS about #pragma in optimize.c
if(UNIX AND NOT APPLE)
 SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
endif() 

if(UNIX AND APPLE)
 SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Xpreprocessor -fopenmp")
endif() 

if(APPLE)
 SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Xpreprocessor -fopenmp")
endif()


# add link to include
include_directories(${CMAKE_SOURCE_DIR}/include /usr/local/include/) 

# add subdirectory lib
add_subdirectory(${CMAKE_SOURCE_DIR}/lib)

# add the executable
add_executable(excogito excogito.c)

# linking library
target_link_libraries(excogito PUBLIC lib)
 
if(UNIX AND NOT APPLE)
 target_link_libraries(excogito PRIVATE m)
 SET(GCC_COVERAGE_COMPILE_FLAGS "-Wall -fopenmp") 
endif()

if(UNIX AND APPLE)
 target_link_libraries(excogito PRIVATE m argp omp -L/usr/local/lib)
 SET(GCC_COVERAGE_COMPILE_FLAGS "-Wall -Xpreprocessor -fopenmp")
endif()

if(APPLE)
 target_link_libraries(excogito PRIVATE m argp omp -L/usr/local/lib)
 SET(GCC_COVERAGE_COMPILE_FLAGS "-Wall -Xpreprocessor -fopenmp")
endif()

SET(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

# putting all pieces together
target_include_directories(excogito PUBLIC
                         "${PROJECT_BINARY_DIR}"
                         "${PROJECT_SOURCE_DIR}/lib"
                          )