# We compile everything using gfortran
cmake_minimum_required(VERSION 3.16)

project(mcnpt)
enable_language(Fortran)

# gfortran flags for production runs
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-Wall -Wextra -Wimplicit-interface -fPIC -Werror -fmax-errors=1 -Ofast -march=native -funroll-loops")
endif()

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

# Look for and add the OpenMP flags
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_Fortran_FLAGS "${OpenMP_Fortran_FLAGS} ${dialect}")
endif()

#
# Compilation step
#
# This line will gather all the files and compile them
file(GLOB_RECURSE sources include/*.f90 src/*.f90)
# The executable here is `mcnpt`
add_executable(mcnpt ${sources})