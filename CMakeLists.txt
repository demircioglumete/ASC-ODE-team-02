cmake_minimum_required(VERSION 3.20)
project(ASC_ode)

set (CMAKE_CXX_STANDARD 20)


include_directories(src nanoblas/src)

add_subdirectory (src)
add_subdirectory (nanoblas)


add_executable(test_ode_new demos/test_ode_new.cpp)
target_link_libraries(test_ode_new PUBLIC nanoblas)

add_executable (demo_autodiff demos/demo_autodiff.cpp)
target_link_libraries (demo_autodiff PUBLIC nanoblas)

add_executable(test_ode_circuit demos/test_ode_circuit.cpp)
target_link_libraries(test_ode_circuit PUBLIC nanoblas)
