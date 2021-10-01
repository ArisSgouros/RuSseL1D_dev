cmake_minimum_required (VERSION 2.6)

#set(compiler "intel")
set(compiler "gnu")

set(CMAKE_BUILD_TYPE "Release")
#set(CMAKE_BUILD_TYPE "Debug")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

if (compiler MATCHES "gnu")
  set(CMAKE_Fortran_COMPILER /usr/bin/gfortran)
  message("Selected ${compiler} compiler.")
  if (CMAKE_BUILD_TYPE MATCHES "Release")
    set(CMAKE_Fortran_FLAGS "-O3 -ffree-line-length-none")
    message("Release configuration successful.")
  elseif(CMAKE_BUILD_TYPE MATCHES "Debug")
    set(CMAKE_Fortran_FLAGS "-O0  -g -fcheck=all -Wextra -m64 -pedantic-errors -frepack-arrays -fdump-core -fbounds-check -fimplicit-none -fbacktrace -frange-check -Waliasing -Wampersand -Wsurprising -Wunderflow -W -ffree-line-length-none")
    message("Debug configuration successful.")
  else()
    message("ERROR: select a compiler configuration to proceed.")
  endif()
elseif(compiler MATCHES "intel")
  set(CMAKE_Fortran_COMPILER ifort)
  message("Selected ${compiler} compiler!")
  if (CMAKE_BUILD_TYPE MATCHES "Release")
    set(CMAKE_Fortran_FLAGS "-O2 -prec-sqrt -prec-div -align -static -ip -ipo -heap-arrays")
    message("Release configuration successful.")
  elseif(CMAKE_BUILD_TYPE MATCHES "Debug")
    set(CMAKE_Fortran_FLAGS "-O0 -g -traceback -fpe0 -fp-stack-check -heap-arrays -ftrapuv -check pointers -check bounds -warn all")
    message("Debug configuration successful.")
  else()
    message("ERROR: select a compiler configuration to proceed.")
  endif()
else()
  message("ERROR: select a compiler type to proceed.")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../run)

enable_language(Fortran)
project (RuSseL1D)

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
