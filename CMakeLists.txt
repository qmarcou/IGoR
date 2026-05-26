### CMakeLists.txt ---

cmake_minimum_required(VERSION 3.31.0)

######################################################################

project(igor)

## ###################################################################
## cmake modules folder
## ###################################################################

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

## ###################################################################
## Version setup
## ###################################################################

set(${PROJECT_NAME}_VERSION_MAJOR 1)
set(${PROJECT_NAME}_VERSION_MINOR 4)
set(${PROJECT_NAME}_VERSION_PATCH 0)
set(${PROJECT_NAME}_VERSION
  ${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH})

## ###################################################################
## Output directory setup
## ###################################################################

include(GNUInstallDirs)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

## ###################################################################
## Default build type (RelWithDebInfo)
## ###################################################################

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'RelWithDebInfo' as none was specified.")
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

## #################################################################
## Install prefix
## #################################################################

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  if (DEFINED ENV{CONDA_PREFIX})
    set(CMAKE_INSTALL_LIBDIR "lib" CACHE PATH "${PROJECT_NAME} libdir" FORCE)
    set(CMAKE_INSTALL_PREFIX "$ENV{CONDA_PREFIX}" CACHE PATH "${PROJECT_NAME} install prefix" FORCE)
  else()
    set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install" CACHE PATH "${PROJECT_NAME} install prefix" FORCE)
  endif()
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

if(DEFINED ENV{PIXI_PROJECT_ROOT})
  set(IGOR_DATA_DIR ${CMAKE_SOURCE_DIR})
else()
  set(IGOR_DATA_DIR "${CMAKE_INSTALL_PREFIX}/share/igor")
endif()

# Install the demo directory to the data installation location
install(DIRECTORY ${CMAKE_SOURCE_DIR}/demo
        DESTINATION share/igor
        USE_SOURCE_PERMISSIONS)

## #############################################################################
## Set RPATH to ensure that both installed libraries and executables
## can find their dependencies at runtime
## #############################################################################

set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
list(APPEND CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
list(REMOVE_DUPLICATES CMAKE_INSTALL_RPATH)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

## #################################################################
## Generate compilation database
## #################################################################

set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

## #################################################################
## Build setup
## #################################################################

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)

## ###################################################################
## Test setup
## ###################################################################

include(CTest)
enable_testing()

## ###################################################################
## Dependencies - cmake
## ###################################################################

include(GenerateExportHeader)

## #################################################################
## Dependencies - external
## #################################################################

find_package(GSL REQUIRED)

if (APPLE)
  find_library(OpenMP_LIBRARY
    NAMES omp
  )
  find_path(OpenMP_INCLUDE_DIR
    omp.h
    PATHS
    $ENV{CONDA_PREFIX}/lib/clang/4.0.1/include
  )
  mark_as_advanced(OpenMP_LIBRARY OpenMP_INCLUDE_DIR)
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(OpenMP DEFAULT_MSG
    OpenMP_LIBRARY OpenMP_INCLUDE_DIR)

  if (OpenMP_FOUND)
    set(OpenMP_LIBRARIES ${OpenMP_LIBRARY})
    set(OpenMP_INCLUDE_DIRS ${OpenMP_INCLUDE_DIR})
    set(OpenMP_COMPILE_OPTIONS -Xpreprocessor -fopenmp)

    add_library(OpenMP::OpenMP SHARED IMPORTED)
    set_target_properties(OpenMP::OpenMP PROPERTIES
      IMPORTED_LOCATION ${OpenMP_LIBRARIES}
      INTERFACE_INCLUDE_DIRECTORIES "${OpenMP_INCLUDE_DIRS}"
      INTERFACE_COMPILE_OPTIONS "${OpenMP_COMPILE_OPTIONS}"
    )
  endif()
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lpthread")
else()
  find_package(OpenMP REQUIRED)
endif()

## #############################################################################
## Options
## #############################################################################

option(BUILD_DOC_DOXYGEN "Build Doxygen Documentation." false)
option(ENABLE_COVERAGE "Enable code coverage." false)
option(IGOR_WITH_SEQAN2 "Enable SeqAn2 banded alignment backend" OFF)
option(IGOR_SEQAN2_SIMD "Compile SeqAn2 with SIMD vectorization (requires -msse4 or -mavx2)" OFF)

if(IGOR_WITH_SEQAN2)
  if(DEFINED ENV{CONDA_PREFIX})
    list(APPEND CMAKE_MODULE_PATH "$ENV{CONDA_PREFIX}/share/seqan/cmake")
  endif()
  find_package(SeqAn REQUIRED)
endif()

if(ENABLE_COVERAGE)
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
#    add_compile_options(Wall -Wextra -pedantic-errors) # To activate once it compiles without -pedantic-errors
    add_compile_options(-Wall -Wextra)
  endif()
  if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    add_compile_options(-fprofile-generate)
    add_link_options(-fprofile-generate)
  elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    add_compile_options(-fprofile-instr-generate -fcoverage-mapping)
    add_link_options(-fprofile-instr-generate -fcoverage-mapping)
  elseif(WIN32)
    add_compile_options(/FS)
  else()
    message(FATAL_ERROR "Unsupported system: ${CMAKE_SYSTEM_NAME}")
  endif()
endif(ENABLE_COVERAGE)
#add_subdirectory(tst)
#

## #################################################################
## Formatting
## #################################################################

find_program(CLANG_FORMAT clang-format)
if(CLANG_FORMAT)
  file(GLOB_RECURSE ALL_SOURCE_FILES
    "${CMAKE_SOURCE_DIR}/src/*.cpp" "${CMAKE_SOURCE_DIR}/src/*.h"
    "${CMAKE_SOURCE_DIR}/app/*.cpp" "${CMAKE_SOURCE_DIR}/app/*.h"
    "${CMAKE_SOURCE_DIR}/tst/*.cpp" "${CMAKE_SOURCE_DIR}/tst/*.h"
  )
  add_custom_target(clang-format
    COMMAND ${CLANG_FORMAT} -i -style=file ${ALL_SOURCE_FILES}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Running clang-format on source files"
  )
else()
  message(WARNING "clang-format not found, formatting target will not be available.")
endif()

## Input
## ###################################################################

add_subdirectory(src)
add_subdirectory(app)
add_subdirectory(tst)
if(BUILD_DOC_DOXYGEN)
  add_subdirectory(doc)
endif()

## ###################################################################
## Exporting
## ###################################################################

include(${PROJECT_NAME}Export)

export(PACKAGE ${PROJECT_NAME})

## ###################################################################
## Beautifying
## ###################################################################

mark_as_advanced(${PROJECT_NAME}_VERSION_MAJOR)
mark_as_advanced(${PROJECT_NAME}_VERSION_MINOR)
mark_as_advanced(${PROJECT_NAME}_VERSION_BUILD)

mark_as_advanced(CMAKE_AR)
mark_as_advanced(CMAKE_OSX_ARCHITECTURES)
mark_as_advanced(CMAKE_OSX_DEPLOYMENT_TARGET)
mark_as_advanced(CMAKE_OSX_SYSROOT)

######################################################################
### CMakeLists.txt ends here
