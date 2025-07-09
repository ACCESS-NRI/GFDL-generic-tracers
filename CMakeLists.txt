# https://cmake.org/cmake/help/v3.6/module/FindPkgConfig.html
# 3.11: add_library() doesn't require dummy source
# 3.12: target_link_libraries() supports OBJECT targets
# 3.20: cmake_path()
cmake_minimum_required(VERSION 3.20)
project(GFDLGTracers
  VERSION 2024.08.1
  DESCRIPTION "GFDLGTracers contains a collection of tracers and associated code for use with the MOM and GOLD ocean models."
  LANGUAGES Fortran
)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CheckFortranCompilerFlag)

if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo|MinSizeRel)$")
  message(STATUS "Setting build type to 'Relwithdebinfo' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "Relwithdebinfo"
      CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|GNU)$")
  message(
    WARNING "Compiler not officially supported: ${CMAKE_Fortran_COMPILER_ID}")
endif()

cmake_path(APPEND CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR} "include")

option(BUILD_SHARED_LIBS "Build shared/dynamic libraries" OFF)

################################################################################
# Fortran
################################################################################

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")

  # Copied from MOM5/bin/mkmf.template.nci.gfortran
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcray-pointer -fdefault-real-8 -ffree-line-length-none -fno-range-check -Waliasing -Wampersand -Warray-bounds -Wcharacter-truncation -Wconversion -Wline-truncation -Wintrinsics-std -Wsurprising -Wno-tabs -Wunderflow -Wunused-parameter -Wintrinsic-shadow -Wno-align-commons")
  check_fortran_compiler_flag("-fallow-invalid-boz" _boz_flag)
  check_fortran_compiler_flag("-fallow-argument-mismatch" _argmis_flag)
  if(_boz_flag)
    set(CMAKE_Fortran_FLAGS               "${CMAKE_Fortran_FLAGS} -fallow-invalid-boz" )
  endif()
  if(_argmis_flag)
    set(CMAKE_Fortran_FLAGS               "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch" )
  endif()
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -W -fbounds-check")

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")

  # Copied from MOM5/bin/mkmf.template.nci
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-alias -safe-cray-ptr -fpe0 -ftz -assume byterecl -i4 -r8 -traceback -nowarn -check noarg_temp_created -assume nobuffered_io -convert big_endian -grecord-gcc-switches -fp-model precise -fp-model source -align all")
  set(CMAKE_Fortran_FLAGS_RELEASE "-g3 -O2 -xCORE-AVX2 -debug all -check none")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g3 -O2 -xCORE-AVX2 -debug all -check none")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g3 -O0 -debug all -check -check noarg_temp_created -check nopointer -warn -warn noerrors -ftrapuv")

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")

  # Copied from MOM5/bin/mkmf.template.nci
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-alias -safe-cray-ptr -fpe0 -ftz -assume byterecl -i4 -r8 -traceback -nowarn -check noarg_temp_created -assume nobuffered_io -convert big_endian -grecord-gcc-switches -fp-model precise -fp-model source -align all")
  set(CMAKE_Fortran_FLAGS_RELEASE "-g3 -O2 -xCORE-AVX2 -debug all -check none")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g3 -O2 -xCORE-AVX2 -debug all -check none")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g3 -O0 -debug all -check -check noarg_temp_created -check nopointer -warn -warn noerrors -ftrapuv")

else()
  message(WARNING "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()

################################################################################

find_package(MPI REQUIRED COMPONENTS Fortran)
find_package(FMS REQUIRED COMPONENTS R8)
find_package(mocsy REQUIRED)

add_library(gtracers)

target_sources(gtracers PRIVATE
  generic_tracers/FMS_coupler_util.F90
  generic_tracers/FMS_ocmip2_co2calc.F90
  generic_tracers/generic_abiotic.F90
  generic_tracers/generic_age.F90
  generic_tracers/generic_argon.F90
  generic_tracers/generic_BLING.F90
  generic_tracers/generic_blres.F90
  generic_tracers/generic_CFC.F90
  generic_tracers/generic_COBALT.F90
  generic_tracers/generic_ERGOM.F90
  generic_tracers/generic_miniBLING.F90
  generic_tracers/generic_SF6.F90
  generic_tracers/generic_TOPAZ.F90
  generic_tracers/generic_tracer.F90
  generic_tracers/generic_tracer_utils.F90
  generic_tracers/generic_WOMBATlite.F90
  generic_tracers/generic_WOMBATmid.F90
)

set_target_properties(gtracers PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS})

# "Creating CMake Libraries - That others can find and use."
# https://www.youtube.com/watch?v=08f5Dav72aE
target_include_directories(gtracers
  PUBLIC
  # Generic Tracers has no header files but creates modules
  # "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/generic_tracers>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

target_link_libraries(gtracers PUBLIC
			FMS::fms_r8
			mocsy::mocsy)

install(TARGETS gtracers
        EXPORT GFDLGTracersTargets)

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ TYPE INCLUDE)

install(EXPORT GFDLGTracersTargets
  FILE GFDLGTracersTargets.cmake
  NAMESPACE GFDLGTracers::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GFDLGTracers
)

configure_package_config_file(
  cmake/GFDLGTracersConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/GFDLGTracersConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GFDLGTracers
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/GFDLGTracersConfig.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GFDLGTracers
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/GFDLGTracers.pc.in
               ${CMAKE_CURRENT_BINARY_DIR}/GFDLGTracers.pc @ONLY)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/GFDLGTracers.pc
              DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
              COMPONENT utilities)
