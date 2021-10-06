cmake_minimum_required(VERSION 2.6)

project("RuSseL1D" LANGUAGES Fortran) # remove the "LANGUAGES" keyword, if Cmake version is older than 2.6

#If selection is made as command-line argument, both the following lines must be commented out.
#set(CMAKE_BUILD_TYPE "Release" CACHE STRING "release_compiler_configuration" FORCE)
set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "debug_compiler_configuration" FORCE)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  message("Building project with GNU compiler")

  set(gcc_debug_flags "-O0 -fcheck=all -Wextra -m64 -pedantic-errors -frepack-arrays -fdump-core -fbounds-check -fimplicit-none -fbacktrace -frange-check -Waliasing -Wampersand -Wsurprising -Wunderflow -W -ffree-line-length-none")
  set(CMAKE_Fortran_FLAGS_DEBUG ${gcc_debug_flags} CACHE STRING "gcc_debug_flags" FORCE)

  set(gcc_release_flags "-O3 -ffree-line-length-none")
  set(CMAKE_Fortran_FLAGS_RELEASE ${gcc_release_flags} CACHE STRING "gcc_release_flags" FORCE)
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  message("Building project with Intel compiler")

  set(intel_debug_flags "-O0 -g -traceback -fpe0 -fp-stack-check -heap-arrays -ftrapuv -check pointers -check bounds -warn all")
  set(CMAKE_Fortran_FLAGS_DEBUG ${intel_debug_flags} CACHE STRING "intel_debug_flags" FORCE)

  set(intel_release_flags "-O3 -prec-sqrt -prec-div -align -static -ip -ipo -heap-arrays")
  set(CMAKE_Fortran_FLAGS_RELEASE ${intel_release_flags} CACHE STRING "intel_release_flags" FORCE)
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "MSVC")
  message("Building project with Visual Studio")
endif()

set(CMAKE_BINARY_DIR ./build)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR})

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../run)

add_subdirectory(src)

execute_process(
  COMMAND bash -c "chmod +x ./test_integrity/test_integrity.sh"
  COMMAND bash -c "chmod +x ./clean.sh"
)

add_custom_target(cleansimple rm -f ./run/o.*)
add_custom_target(cleaner rm -f ./run/o.* ./run/fort.* ./run/RSL1D)
add_custom_target(cleantest rm -f ./run/o.* ./run/fort.* ./run/RSL1D ./run/LOG.*)
add_custom_target(cleanout rm -f ./run/o.* ./run/fort.* LOG.*)
add_custom_target(rsltest ./test_integrity/test_integrity.sh)
