cmake_minimum_required(VERSION 3.29)
set(CMAKE_CXX_COMPILER clang++)
PROJECT(SignalServer23 LANGUAGES CXX)

#include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
#conan_basic_setup()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Enable compile commands export")

set(CMAKE_CXX_STANDARD 23)

# Check for println support
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "14.0")
    message(FATAL_ERROR "gcc version too old, will not support println")
endif()
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "18.0")
    message(FATAL_ERROR "clang version too old, will not support println")
endif()

find_package(BZip2 REQUIRED)
find_package(ZLIB REQUIRED)

option(ENABLE_CLANG_TIDY "Enable clang-tidy checks" ON)

if(ENABLE_CLANG_TIDY)
  find_program(CLANG_TIDY_PATH clang-tidy REQUIRED)
  function(AddClangTidy target)
    set_target_properties(${target}
      PROPERTIES CXX_CLANG_TIDY
      "${CLANG_TIDY_PATH};--config-file=${CMAKE_SOURCE_DIR}/.clang-tidy"
      )
  endfunction()
endif()

add_subdirectory(src)

## TODO: Add tests
enable_testing()
#add_subdirectory(test)
