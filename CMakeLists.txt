cmake_minimum_required(VERSION 3.30)
project(RSA)
set(CMAKE_CXX_STANDARD 20)

set(SOURCE_FILES main.cpp)
add_executable(RSA ${SOURCE_FILES})

if (UNIX AND NOT APPLE)
    ## Build for Linux
    target_link_libraries(RSA gmp gmpxx)

elseif (WIN32)
    ## Build for Windows
    include_directories(${CMAKE_SOURCE_DIR}/extern/gmp/include)
    link_directories(${CMAKE_SOURCE_DIR}/extern/gmp/lib)

    # Enable static linking of C++ and GCC runtime libraries
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++ -static")

    target_link_libraries(RSA PRIVATE
            ${CMAKE_SOURCE_DIR}/extern/gmp/lib/libgmpxx.a
            ${CMAKE_SOURCE_DIR}/extern/gmp/lib/libgmp.a
    )

endif()