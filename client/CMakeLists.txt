cmake_minimum_required(VERSION 3.10)

project(
  eonclient
  VERSION 1.0
  DESCRIPTION "EON client"
  LANGUAGES C CXX Fortran)

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n"
  )
endif()

set(CMAKE_COLOR_MAKEFILE ON)

# Compiler Flags
set(CMAKE_Fortran_FLAGS "-fbounds-check -ffpe-trap=invalid,zero,overflow"
)# -qopenmp
set(CMAKE_Fortran_FLAGS_RELEASE
    "-O3 -march=native -ffast-math -funroll-loops -std=legacy -fopenmp")
set(CMAKE_Fortran_FLAGS_DEBUG
    "-O0 -Warray-bounds  -std=legacy -Wall -Wno-unused ")
set(DOUBLE "-r8")

add_compile_options(# "-Wall" "-Wpedantic" "-Wextra" "-fexceptions"
                    "$<$<CONFIG:Debug>:-O0;-g3;-ggdb>")

# Messages
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})
message(STATUS "Fortran compiler: " ${CMAKE_Fortran_COMPILER})
message(STATUS "C compiler: " ${CMAKE_C_COMPILER})
message(STATUS "CXX compiler: " ${CMAKE_CXX_COMPILER})
message(STATUS "Python: " ${Python_EXECUTABLE})

if(${CMAKE_BUILD_TYPE} MATCHES "Release")
  message(STATUS "Fortran compiler flags:" ${CMAKE_Fortran_FLAGS_RELEASE})
  message(STATUS "CXX compiler flags:" ${CMAKE_CXX_FLAGS_RELEASE})
elseif(${CMAKE_BUILD_TYPE} MATCHES "Debug")
  message(STATUS "Fortran compiler flags:" ${CMAKE_Fortran_FLAGS_DEBUG})
  message(STATUS "CXX compiler flags:" ${CMAKE_CXX_FLAGS_DEBUG})
endif()

# Build
add_executable(eonclient ClientEON.cpp CommandLine.cpp)

# Versioning
add_custom_target(
  version
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/version.sh >
          ${CMAKE_CURRENT_SOURCE_DIR}/version.h
  COMMENT "run version update command in ${CMAKE_CURRENT_SOURCE_DIR}")

# Build order
add_dependencies(eonclient potentials eonlib version gprdimer)

# Libraries
target_link_libraries(eonclient gprdimer potentials eonlib)

# Sub-directories
add_subdirectory(gprdimer)
add_subdirectory(potentials)

# Dependencies
add_library(
  eonlib SHARED
  MinModeSaddleSearch.cpp
  Dimer.cpp
  EpiCenters.cpp
  Hessian.cpp
  ConjugateGradients.cpp
  HelperFunctions.cpp
  Matter.cpp
  Parameters.cpp
  INIFile.cpp
  Potential.cpp
  Quickmin.cpp
  ProcessSearchJob.cpp
  PointJob.cpp
  MinimizationJob.cpp
  HessianJob.cpp
  ParallelReplicaJob.cpp
  SafeHyperJob.cpp
  TADJob.cpp
  ReplicaExchangeJob.cpp
  Dynamics.cpp
  BondBoost.cpp
  FiniteDifferenceJob.cpp
  NudgedElasticBandJob.cpp
  TestJob.cpp
  BasinHoppingJob.cpp
  SaddleSearchJob.cpp
  ImprovedDimer.cpp
  NudgedElasticBand.cpp
  Lanczos.cpp
  Bundling.cpp
  Job.cpp
  CommandLine.cpp
  DynamicsJob.cpp
  Log.cpp
  LBFGS.cpp
  LowestEigenmode.cpp
  Optimizer.cpp
  Prefactor.cpp
  DynamicsSaddleSearch.cpp
  PrefactorJob.cpp
  FIRE.cpp
  GlobalOptimizationJob.cpp
  GlobalOptimization.cpp
  StructureComparisonJob.cpp
  MonteCarloJob.cpp
  MonteCarlo.cpp
  SteepestDescent.cpp
  BasinHoppingSaddleSearch.cpp
  BiasedGradientSquaredDescent.cpp
  AtomicGPDimer.cpp)

# Options
if(MPI)
  target_compile_definitions(eonclient PUBLIC -DEONMPI)
  # Python is only used with MPI
  find_package(Python COMPONENTS Interpreter Development)
  include_directories(${Python_INCLUDE_DIRS})
  # Need the MPI libraries
  find_package(MPI REQUIRED)
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
  set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
  target_link_libraries(eonclient ${Python_LIBRARIES} MPI::MPI_CXX)
endif()

if(UNIX)
  target_compile_definitions(eonclient PUBLIC -DLINUX)
endif()

if(OSX)
  target_compile_definitions(eonclient PUBLIC -DOSX)
endif()

if(WIN32)
  target_compile_definitions(eonlib PUBLIC -DWIN32) # No VASP
  target_compile_definitions(eonclient PUBLIC -DWIN32)
endif()

if(NO_FORTRAN)
  target_compile_definitions(eonlib PUBLIC -DNO_FORTRAN)
endif()

# Tests
option(PACKAGE_TESTS "Build the tests" OFF)
if(PACKAGE_TESTS)
  find_package(GTest REQUIRED)
  include(GoogleTest)
  add_subdirectory(gtests)
endif()