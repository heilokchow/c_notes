# CMakeList.txt : CMake project for CMake_notes, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Add source to this project's executable.
add_executable (CMake_notes "main.cpp" "template.h")

# TODO: Add tests and install targets if needed.