cmake_minimum_required(VERSION 3.15)
project(quadratic_sieve)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMake/Modules)
find_package(GMP REQUIRED)

INCLUDE_DIRECTORIES(
        "${GMP_INCLUDES}"
)


if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2  -flto -Wall -pedantic")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address,undefined -Wall -pedantic -D_GLIBCXX_DEBUG")
endif()


add_executable(quadratic_sieve
        main.cpp
        quadratic_sieve.hpp
        quadratic_sieve.cpp
        )

target_link_libraries(quadratic_sieve ${GMP_LIBRARIES})