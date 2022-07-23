cmake_minimum_required (VERSION 3.13)

project (RSymSQCD VERSION 0.9.0 LANGUAGES CXX)
configure_file(include/version.hpp.in include/version.hpp)

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_STANDARD_REQUIRED True)
set (CMAKE_CXX_EXTENSIONS OFF)            # use -std=c++ instead of -std=gnu

# user should build from his own directory, e.g. ./build
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are not allowed.")
endif()

# build type
set (default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message (STATUS "Setting build type to '${default_build_type}' as none was specified")
  set (CMAKE_BUILD_TYPE "${default_build_type}")
endif()

include(FetchContent)
FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.10.5/json.tar.xz)
FetchContent_MakeAvailable(json)

find_package(Boost REQUIRED COMPONENTS program_options)
find_package (Eigen3)

find_path(LIBCUBA_INCLUDE_DIR NAMES cuba.h)
find_library(LIBCUBA_LIBRARY NAMES libcuba.a)

find_path(LIBLOOPTOOLS_INCLUDE_DIR NAMES clooptools.h)
find_library(LIBLOOPTOOLS_LIBRARY NAMES libooptools.a)

find_path(LIBLHAPDF_INCLUDE_DIR NAMES LHAPDF/LHAPDF.h)
if (APPLE)
   find_library(LIBLHAPDF_LIBRARY NAMES libLHAPDF.dylib)
elseif (UNIX)
   find_library(LIBLHAPDF_LIBRARY NAMES libLHAPDF.so)
endif()

find_path(LIBRK_INCLUDE_DIR NAMES rk/rk.hh)
if(APPLE)
   find_library(LIBRK_LIBRARY NAMES librk.dylib)
elseif(UNIX)
   find_library(LIBRK_LIBRARY NAMES librk.so)
endif()

# find_path(LIBCOLORFULL_INCLUDE_DIRS NAMES ColorFull/Core/Tree_level_gluon_basis.h)
# if(APPLE)
#    find_library(LIBCOLORFULL_LIBRARY NAMES ColorFull/libColorFull.dylib)
# elseif(UNIX)
#    find_library(LIBCOLORFULL_LIBRARY NAMES ColorFull/libColorFull.so)
# endif()

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
endif()

# we link against LoopTools which requires a fortran library
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
   # assume that if we're using GCC or Clang then LoopTools was compiled with GCC
   if(APPLE)
      find_library(FORTRAN_LIBS NAMES "libgfortran.dylib")
   elseif(UNIX)
      find_library(FORTRAN_LIBS NAMES "libgfortran.so.5")
   endif()
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
   set(FORTRAN_LIBS "ifcore;ifport")
endif()

file (GLOB SOURCES
   "src/Li2.cpp"  "src/Process.cpp"  "src/splitting_kernels.cpp"  "XSection*.cpp"
)
add_library(RSymSQCDlib SHARED ${SOURCES})
target_include_directories (
      RSymSQCDlib
   PUBLIC
      include
      ${LIBCUBA_INCLUDE_DIR}
      ${LIBLOOPTOOLS_INCLUDE_DIR}
      ${LIBLHAPDF_INCLUDE_DIR}
      ${LIBRK_INCLUDE_DIR}
      ${Boost_INCLUDE_DIRS}
      # ${LIBCOLORFULL_INCLUDE_DIRS}
)
target_link_libraries(
   RSymSQCDlib
   Boost::program_options
   ${LIBCUBA_LIBRARY}
   ${LIBLOOPTOOLS_LIBRARY}
   ${FORTRAN_LIBS}
   ${LIBLHAPDF_LIBRARY}
   # ${LIBCOLORFULL_LIBRARY}
   Eigen3::Eigen
)

add_executable(RSymSQCD "src/main.cpp" "src/XSection.cpp" "src/XSection_Tree.cpp" "src/XSection_Virt.cpp" "src/XSection_SC.cpp" "src/XSection_HnonC.cpp" "src/Process.cpp")
target_include_directories(
      RSymSQCD
   PRIVATE
      include
      ${CMAKE_CURRENT_BINARY_DIR}/include # location of version.hpp
      ${LIBCUBA_INCLUDE_DIR}
      ${LIBLOOPTOOLS_INCLUDE_DIR}
      ${LIBLHAPDF_INCLUDE_DIR}
      ${LIBRK_INCLUDE_DIR}
      ${Boost_INCLUDE_DIRS}
      # ${LIBCOLORFULL_INCLUDE_DIRS}
)
target_link_libraries(
      RSymSQCD
   PRIVATE
      RSymSQCDlib
      nlohmann_json::nlohmann_json
      Boost::program_options
      ${LIBCUBA_LIBRARY}
      ${LIBLOOPTOOLS_LIBRARY}
      ${FORTRAN_LIBS}
      ${LIBLHAPDF_LIBRARY}
      # ${LIBCOLORFULL_LIBRARY}
      Eigen3::Eigen
)

if(ENABLE_TESTS)
   enable_testing()
   add_subdirectory(test)
endif()
