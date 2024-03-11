cmake_minimum_required(VERSION 3.10)

project(CanonicalExpansion)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(GMP_PREFIX "" CACHE PATH "path to GMP")

find_path(GMP_INCLUDE_DIR gmp.h gmpxx.h
    PATHS ${GMP_PREFIX}/include /usr/include /usr/local/include)

find_library(GMP_LIBRARY NAMES gmp libgmp
    PATHS ${GMP_PREFIX}/lib /usr/lib /usr/local/lib)

if(GMP_INCLUDE_DIR AND GMP_LIBRARY)
    set(GMP_FOUND TRUE)
    message(STATUS "Found GMP: ${GMP_LIBRARY}")
    set(GMP_INCLUDE_DIRS ${GMP_INCLUDE_DIR})
    set(GMP_LIBRARIES ${GMP_LIBRARY})
else()
    if(GMP_FIND_REQUIRED)
        message(FATAL_ERROR "Could not find GMP")
    endif()
endif()

mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARY)

if(NOT GMP_FOUND)
    message(FATAL_ERROR "GMP not found")
endif()

set(SOURCES
    src/main.cpp
    src/algorithms.cpp
)

add_executable(canonicalExpansionApp ${SOURCES})

target_include_directories(canonicalExpansionApp PRIVATE ${GMP_INCLUDE_DIRS})
target_link_libraries(canonicalExpansionApp PRIVATE ${GMP_LIBRARIES})