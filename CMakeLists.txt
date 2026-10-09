cmake_minimum_required(VERSION 3.10)
project(LiLi_Test LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

add_library(lili_core src/info_matrix.cpp)
target_include_directories(lili_core
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${EIGEN3_INCLUDE_DIRS}
)
target_link_libraries(lili_core PUBLIC Eigen3::Eigen)

add_executable(lili_synthetic src/synthetic.cpp)
target_link_libraries(lili_synthetic PRIVATE lili_core)

add_executable(lili_cli src/cli.cpp)
target_link_libraries(lili_cli PRIVATE lili_core)

enable_testing()
add_executable(lili_tests tests/test_info_matrix.cpp)
target_link_libraries(lili_tests PRIVATE lili_core)
add_test(NAME lili_info_matrix_tests COMMAND lili_tests)
