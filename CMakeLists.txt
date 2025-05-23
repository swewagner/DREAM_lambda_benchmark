cmake_minimum_required(VERSION 3.7...3.25)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

## CUSTOMISE

# Define the application name and version.
project(Dream_lambda_benchmark VERSION 1.0
                  DESCRIPTION "A project to benchmark lambda like approaches"
                  LANGUAGES CXX)

## BUILD

# Make Release default build type
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING
         "Choose the type of build, options are: Debug Release RelWithDebInfo"
         FORCE)
endif ()

# Specify the directories where to store  executables
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Messages
string (ASCII 27 Esc)
set (FontBold "${Esc}[1m")
set (FontReset "${Esc}[m")

# Dependency: SeqAn3.
set(SEQAN3_CLONE_DIR "${CMAKE_SOURCE_DIR}/extern/seqan3")
find_package (SeqAn3 REQUIRED HINTS extern/seqan3/build_system)

# Dependency: SeqAn
# Note: includes seqan-config.cmake manually to have
# SEQAN_DEFINITIONS, SEQAN_CXX_FLAGS_LIST, SEQAN_LIBRARIES, SEQAN_INCLUDE_DIRS
# defined. We don't want any seqan version but the submodule one!
set(CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/extern/seqan/util/cmake")
set(CMAKE_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/extern/seqan/include/")
include (${CMAKE_SOURCE_DIR}/extern/seqan/util/cmake/seqan-config.cmake)

# Define cmake configuration flags to configure and build external projects with the same flags as specified for
# this project.
set (APP_TEMPLATE_EXTERNAL_PROJECT_CMAKE_ARGS "")
list (APPEND APP_TEMPLATE_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
list (APPEND APP_TEMPLATE_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
list (APPEND APP_TEMPLATE_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}")
list (APPEND APP_TEMPLATE_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}")

# Dependancy: Mason2
set (mason2_args ${APP_TEMPLATE_EXTERNAL_PROJECT_CMAKE_ARGS})
list (APPEND mason2_args "-DSEQAN_BUILD_SYSTEM=APP:mason2")

include (ExternalProject)
ExternalProject_Add (
    mason2
    PREFIX mason2
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/extern/seqan"
    CMAKE_ARGS "${mason2_args}"
    INSTALL_COMMAND ""
)
unset (mason2_args)

# Dependancy: Lambda
include (ExternalProject)
ExternalProject_Add (
    lambda3
    PREFIX lambda3
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/extern/lambda"
    CMAKE_ARGS "${APP_TEMPLATE_EXTERNAL_PROJECT_CMAKE_ARGS}"
    INSTALL_COMMAND ""
)

# own scripts
add_subdirectory (src)

message (STATUS "${FontBold}You can run `make` to build the application.${FontReset}")
