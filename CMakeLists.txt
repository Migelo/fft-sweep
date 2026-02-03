cmake_minimum_required(VERSION 3.16)
project(sweep Fortran)

# Find MPI
find_package(MPI REQUIRED)

# Find FFTW3
find_library(FFTW3_LIB fftw3 HINTS /usr/lib /usr/local/lib)
find_library(FFTW3_MPI_LIB fftw3_mpi HINTS /usr/lib /usr/local/lib)
find_path(FFTW3_INCLUDE_DIR fftw3.f03 HINTS /usr/include /usr/local/include)

if(NOT FFTW3_LIB OR NOT FFTW3_MPI_LIB)
    message(FATAL_ERROR "FFTW3 libraries not found")
endif()

# Executable
add_executable(sweep sweep.f90)

# Include directories
target_include_directories(sweep PRIVATE 
    ${MPI_Fortran_INCLUDE_PATH}
    ${FFTW3_INCLUDE_DIR}
)

# Link libraries
target_link_libraries(sweep 
    ${MPI_Fortran_LIBRARIES}
    ${FFTW3_MPI_LIB}
    ${FFTW3_LIB}
    m
)

# Compiler flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2")
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
endif()
