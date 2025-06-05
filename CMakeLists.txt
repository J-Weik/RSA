cmake_minimum_required(VERSION 3.31)
project(RSA)

set(CMAKE_CXX_STANDARD 20)

add_executable(RSA main.cpp)

target_link_libraries(RSA PRIVATE gmp gmpxx)

