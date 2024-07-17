cmake_minimum_required(VERSION 2.8)

project(efanna2e)
include_directories(${PROJECT_SOURCE_DIR}/include)

# mimalloc
find_package(mimalloc 1.8 REQUIRED NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "Found mimalloc installed at: ${MIMALLOC_LIBRARY_DIR} (${MIMALLOC_VERSION_DIR})")

#OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    message(FATAL_ERROR "no OpenMP supprot")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")

add_definitions (-std=c++11 -O3 -lboost -pthread -march=native -Wall -DINFO)

add_subdirectory(src)
add_subdirectory(tests)
