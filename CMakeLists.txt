# This file is part of blas-interface.
# SPDX-Identifier: Apache-2.0
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

cmake_minimum_required(VERSION 3.14)
get_directory_property(is-subproject PARENT_DIRECTORY)

project(
  "blas-interface"
  LANGUAGES "Fortran"
  VERSION "1.0.0"
  DESCRIPTION "Molecular structure store for testing"
)

# Follow GNU conventions for installing directories
include(GNUInstallDirs)

# General configuration information
add_subdirectory("config")

if(NOT TARGET "OpenMP::OpenMP_Fortran" AND WITH_OpenMP)
  find_package("OpenMP" REQUIRED)
endif()

if(WITH_ILP64)
  message(STATUS "Using LAPACK/BLAS ILP64 interface")
endif()

if(NOT TARGET "BLAS::BLAS")
  find_package("custom-blas" REQUIRED)
endif()

# Collect source of the project
set(srcs)
add_subdirectory("src")

# multicharge library target
add_library(
  "${PROJECT_NAME}-lib"
  "${srcs}"
)
set_target_properties(
  "${PROJECT_NAME}-lib"
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME "${PROJECT_NAME}"
  VERSION "${PROJECT_VERSION}"
  SOVERSION "${PROJECT_VERSION_MAJOR}"
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include"
)
target_link_libraries(
  "${PROJECT_NAME}-lib"
  PUBLIC
  "BLAS::BLAS"
  "$<$<BOOL:${WITH_OpenMP}>:OpenMP::OpenMP_Fortran>"
)
target_compile_definitions(
  "${PROJECT_NAME}-lib"
  PUBLIC "-DIK=$<IF:$<BOOL:${WITH_ILP64}>,i8,i4>"
)

target_include_directories(
  "${PROJECT_NAME}-lib"
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}>
)

# Export targets for other projects
add_library("${PROJECT_NAME}" INTERFACE)
target_link_libraries("${PROJECT_NAME}" INTERFACE "${PROJECT_NAME}-lib")
install(
  TARGETS
  "${PROJECT_NAME}"
  "${PROJECT_NAME}-lib"
  EXPORT
  "${PROJECT_NAME}-targets"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
install(
  EXPORT
  "${PROJECT_NAME}-targets"
  NAMESPACE
  "${PROJECT_NAME}::"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
install(
  DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}"
)
# Package license files
install(
  FILES
  "LICENSE"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/${PROJECT_NAME}"
)
