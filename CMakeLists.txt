cmake_minimum_required(VERSION 3.10)

project(CanonicalExpansion)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}")

include(FindGMP)
find_package(GMP REQUIRED)

set(SOURCES
    src/main.cpp
    src/algorithms.cpp
)

add_executable(canonicalExpansionApp ${SOURCES})

# target_link_libraries(canonicalExpansionApp gmp gmpxx)
target_include_directories(canonicalExpansionApp PRIVATE ${GMP_INCLUDE_DIRS})
target_link_libraries(canonicalExpansionApp PRIVATE ${GMP_LIBRARIES})