# CMakeList.txt : CMake project for CMake_notes, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Add source to this project's executable.
add_executable (CMake_notes "main.cpp")

target_compile_options(CMake_notes PUBLIC -mfma -std=c++11)
target_compile_definitions(CMake_notes PUBLIC MAIN)
target_link_libraries (CMake_notes Eigen3::Eigen OpenMP::OpenMP_CXX)

# TODO: Add tests and install targets if needed.