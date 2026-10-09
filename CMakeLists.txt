cmake_minimum_required(VERSION 3.10)
project(LiLi_Test LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

add_library(lili_core
    src/info_matrix.cpp
    src/se3.cpp
    src/point_cloud.cpp
    src/scenes.cpp
    src/eval_quality.cpp
    src/icp.cpp
    src/zhang_detector.cpp
)
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

add_executable(lili_scenes src/scene_report.cpp)
target_link_libraries(lili_scenes PRIVATE lili_core)

add_executable(lili_eval src/icp_eval.cpp)
target_link_libraries(lili_eval PRIVATE lili_core)

enable_testing()
add_executable(lili_tests tests/test_info_matrix.cpp)
target_link_libraries(lili_tests PRIVATE lili_core)
add_test(NAME lili_info_matrix_tests COMMAND lili_tests)

add_executable(lili_se3_tests tests/test_se3.cpp)
target_link_libraries(lili_se3_tests PRIVATE lili_core)
add_test(NAME lili_se3_tests COMMAND lili_se3_tests)

add_executable(lili_quality_tests tests/test_quality_metric.cpp)
target_link_libraries(lili_quality_tests PRIVATE lili_core)
add_test(NAME lili_quality_tests COMMAND lili_quality_tests)

# M1 acceptance: the analytic degeneracy basis must keep the extended scan on
# the surface while a control direction off that basis must not.
add_test(NAME lili_scene_acceptance COMMAND lili_scenes)

add_executable(lili_icp_tests tests/test_icp.cpp)
target_link_libraries(lili_icp_tests PRIVATE lili_core)
add_test(NAME lili_icp_tests COMMAND lili_icp_tests)

# M2 acceptance: the shared ICP must place every scene inside its known
# degeneracy subspace, and the Zhang baseline must be reported against it.
add_test(NAME lili_m2_baseline COMMAND lili_eval --q-samples 8)
