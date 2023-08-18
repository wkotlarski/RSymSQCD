cmake_minimum_required (VERSION 3.13)

project (RSymSQCD VERSION 1.0.0 LANGUAGES CXX)
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

string(REGEX REPLACE "-g" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}" )
string(APPEND CMAKE_CXX_FLAGS_DEBUG " -O1")

# processor specific optimizations
option(OPTIMIZE_FOR_NATIVE "Optimize for native" FALSE)
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
   if (${OPTIMIZE_FOR_NATIVE})
      include(CheckCXXCompilerFlag)
      if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math -fno-finite-math-only")
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
if ("${ENABLE_IPO}" STREQUAL "")
   set(ENABLE_IPO TRUE)
endif()
# IPO crashes icpc
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release" AND "${ENABLE_IPO}" AND NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
   include(CheckIPOSupported)
   check_ipo_supported(RESULT ipoWorks)
   if(${ipoWorks})
      message(STATUS "Enabling interprocedural optimization (IPO/LTO)")
      set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
   endif()
else()
   set(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE)
endif()

if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
   add_compile_options(-Wall -Wextra -pedantic)
endif()

find_package(nlohmann_json 3.0.0 CONFIG)
if (NOT nlohmann_json_FOUND)
   include(FetchContent)
   if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.24.0")
      cmake_policy(SET CMP0135 NEW)
   endif()
   set(JSON_MultipleHeaders ON)
   FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz)
   FetchContent_MakeAvailable(json)
endif()

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
find_package(TBB REQUIRED)
find_package(spdlog REQUIRED)

foreach (MODEL MSSM MRSSM Sgluons)
   add_library(${MODEL} SHARED "src/models/${MODEL}/${MODEL}.cpp" "src/Li2.cpp")
   target_include_directories (
      ${MODEL}
      PUBLIC
         include
         include/Li2
         ${EIGEN3_INCLUDE_DIR}
         ${Boost_INCLUDE_DIRS}
         ${LIBLOOPTOOLS_INCLUDE_DIR}
   )
   target_link_libraries(
      ${MODEL}
         PUBLIC
         ${FORTRAN_LIBS}
         ${LIBLOOPTOOLS_LIBRARY}
         spdlog::spdlog
   )
endforeach()

add_executable(RSymSQCD "src/main.cpp" "src/XSections/XSection_Tree.cpp" "src/XSections/XSection_Virt.cpp" "src/XSections/XSection_SC.cpp" "src/XSections/XSection_HnonC.cpp" "src/utils.cpp" "src/splitting_kernels.cpp")
target_include_directories(
      RSymSQCD
   PUBLIC
      include
      ${LIBCUBA_INCLUDE_DIR}
      ${LIBRK_INCLUDE_DIR}
      ${LIBLOOPTOOLS_INCLUDE_DIR}
      ${LIBLHAPDF_INCLUDE_DIR}
      ${CMAKE_CURRENT_BINARY_DIR}/include # location of version.hpp
      ${Boost_INCLUDE_DIRS}
)
target_link_libraries(
      RSymSQCD
   PUBLIC
      MSSM MRSSM Sgluons
      ${LIBCUBA_LIBRARY}
      ${LIBLHAPDF_LIBRARY}
      ${LIBLOOPTOOLS_LIBRARY}
      Boost::program_options
      nlohmann_json::nlohmann_json
      TBB::tbb
      spdlog::spdlog
)

if(ENABLE_TESTS)
   enable_testing()
   add_subdirectory(test)
endif()
