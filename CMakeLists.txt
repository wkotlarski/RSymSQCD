cmake_minimum_required (VERSION 3.13)

project (RSymSQCD VERSION 0.9.0 LANGUAGES CXX)
configure_file(include/version.hpp.in include/version.hpp)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED True)
# use -std=c++ instead of -std=gnu
set (CMAKE_CXX_EXTENSIONS OFF)

# user should build from his own directory, e.g. ./build
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are not allowed.")
endif()

# build type
set (default_build_type "Release")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message (STATUS "Setting build type to '${default_build_type}' as none was specified")
  set (CMAKE_BUILD_TYPE "${default_build_type}")
endif()

# processor specific optimizations
option(OPTIMIZE_FOR_NATIVE "Optimize for native" OFF)
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
   if (${OPTIMIZE_FOR_NATIVE})
      include(CheckCXXCompilerFlag)
      if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")
         CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
         if(${COMPILER_SUPPORTS_MARCH_NATIVE})
            message (STATUS "Optimizing for native (adding '-march=native' to CMAKE_CXX_FLAGS)")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
         endif()
      elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fp-model fast")
         CHECK_CXX_COMPILER_FLAG("-xHost" COMPILER_SUPPORTS_xHost)
         if(${COMPILER_SUPPORTS_xHost})
            message (STATUS "Optimizing for native (adding -xHost)")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHost")
         endif()
      endif()
   endif()
endif()

# interprocedural optimizations
# IPO crashes icpc
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release" AND NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
   include(CheckIPOSupported)
   check_ipo_supported(RESULT ipoWorks)
   if(${ipoWorks})
      message(STATUS "Enabling interprocedural optimization (IPO/LTO)")
      set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
   endif()
endif()

if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
   add_compile_options(-Wall -Wextra -pedantic)
endif()

include(FetchContent)
if(${CMAKE_VERSION} VERSION_EQUAL "3.24.0" OR ${CMAKE_VERSION} VERSION_GREATER "3.24.0")
   cmake_policy(SET CMP0135 NEW)
endif()
set(JSON_MultipleHeaders ON)
FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz)
FetchContent_MakeAvailable(json)

find_package(Boost REQUIRED COMPONENTS program_options)

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
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
   # assume that if we're using GCC or Clang then LoopTools was compiled with GCC
   if(APPLE)
      find_library(FORTRAN_LIBS NAMES "libgfortran.dylib")
   elseif(UNIX)
      find_library(FORTRAN_LIBS NAMES "libgfortran.so.5")
   endif()
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
   set(FORTRAN_LIBS "ifcore;ifport")
endif()

find_package(Eigen3 REQUIRED)

add_library(process OBJECT "src/models/MRSSM/MRSSM.cpp" "src/models/Sgluons/Sgluons.cpp" "src/models/MSSM/MSSM.cpp")
target_compile_options(process PRIVATE -fPIC)
target_include_directories(
      process
   PUBLIC
      include
      ${EIGEN3_INCLUDE_DIR}
      ${LIBLOOPTOOLS_INCLUDE_DIR}
      ${LIBLHAPDF_INCLUDE_DIR}
)

find_package(TBB REQUIRED)

file (GLOB SOURCES
   "src/Li2.cpp"  "src/splitting_kernels.cpp"  "XSection*.cpp" "colour_me.cpp"
)
add_library(RSymSQCDlib SHARED ${SOURCES})
target_include_directories (
      RSymSQCDlib
   PUBLIC
      include
      ${LIBCUBA_INCLUDE_DIR}
      ${LIBRK_INCLUDE_DIR}
      ${Boost_INCLUDE_DIRS}
      ${LIBCOLORFULL_INCLUDE_DIRS}
)
target_link_libraries(
      RSymSQCDlib
   PUBLIC
      process
      ${LIBCUBA_LIBRARY}
      ${LIBLOOPTOOLS_LIBRARY}
      ${FORTRAN_LIBS}
      ${LIBLHAPDF_LIBRARY}
      ${LIBCOLORFULL_LIBRARY}
      TBB::tbb
)

add_executable(RSymSQCD "src/main.cpp" "src/XSection_Tree.cpp" "src/XSection_Virt.cpp" "src/XSection_SC.cpp" "src/XSection_HnonC.cpp" "src/utils.cpp")
target_include_directories(
      RSymSQCD
   PUBLIC
      include
      ${CMAKE_CURRENT_BINARY_DIR}/include # location of version.hpp
      ${Boost_INCLUDE_DIRS}
)
target_link_libraries(
      RSymSQCD
   PUBLIC
      process
      RSymSQCDlib
      Boost::program_options
      nlohmann_json::nlohmann_json
)

if(ENABLE_TESTS)
   enable_testing()
   add_subdirectory(test)
endif()
