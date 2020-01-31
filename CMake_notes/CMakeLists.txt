# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.8)

project ("CMake_notes")

find_package(Eigen3 REQUIRED)

add_definitions(-DMAIN)
add_compile_options(-std=c++11)

# Include sub-projects.
add_subdirectory ("CMake_notes")
