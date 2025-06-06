cmake_minimum_required(VERSION 3.31)
project(RSA)

set(CMAKE_CXX_STANDARD 20)

set(SOURCE_FILES main.cpp)

add_executable(RSA ${SOURCE_FILES})

target_link_libraries(RSA PRIVATE gmp gmpxx)