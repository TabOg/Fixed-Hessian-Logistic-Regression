cmake_minimum_required(VERSION 3.5)
project(FHLR)

set(SEAL_DIR ${CMAKE_SOURCE_DIR}/SEAL/cmake)
add_definitions(-O3 -march=native -mtune=native -pedantic -Wall -Wextra)
find_package(SEAL 3.5 REQUIRED)

add_executable(main databasetools.cpp "Fixed Hessian.cpp" GD.cpp iDASH2017.cpp logregtools.cpp main.cpp NAG.cpp plaintextLR.cpp)
target_link_libraries(main PUBLIC SEAL::seal)
