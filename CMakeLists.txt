cmake_minimum_required(VERSION 3.20)
project(ASC_ode)

set (CMAKE_CXX_STANDARD 20)

include_directories(src)
include_directories(src nanoblas/src)

add_subdirectory (src)
#add_subdirectory (nanoblas)

add_subdirectory (mechsystem)

add_executable(test_ode_new demos/test_ode_new.cpp)
target_include_directories(test_ode_new PUBLIC ${PROJECT_SOURCE_DIR}/nanoblas/src)


add_executable (demo_autodiff demos/demo_autodiff.cpp)
target_include_directories(demo_autodiff PUBLIC ${PROJECT_SOURCE_DIR}/nanoblas/src)

add_executable(test_ode_circuit demos/test_ode_circuit.cpp)
target_include_directories(test_ode_circuit PUBLIC ${PROJECT_SOURCE_DIR}/nanoblas/src)

add_executable(test_explicit_rk demos/test_explicit_rk.cpp)
target_include_directories(test_explicit_rk PUBLIC ${PROJECT_SOURCE_DIR}/nanoblas/src)
