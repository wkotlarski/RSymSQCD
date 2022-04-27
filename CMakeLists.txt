cmake_minimum_required (VERSION 3.8)
project (RSymSQCD)

set (CMAKE_CXX_STANDARD 17)

include(FetchContent)
FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.10.5/json.tar.xz)
FetchContent_MakeAvailable(json)

#add_subdirectory(tests)

if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
   link_directories(/opt/local/lib)
   link_directories($ENV{LT}/lib)
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
   link_directories($ENV{LT}/lib64)
else()
   message("Couldn't recognize system (Darwin or Linux). Continuing hoping for the best")
endif()

find_package(Boost REQUIRED COMPONENTS program_options)
find_package (Eigen3)

include_directories (
   include
   $ENV{COLORFULL}/include
   $ENV{Cuba}
   $ENV{LT}/include
   $ENV{rk}/include
   ${Boost_INCLUDE_DIRS}
)

link_directories(
   $ENV{COLORFULL}/lib/ColorFull
   $ENV{Cuba}/lib
)

file (GLOB SOURCES
   "src/*.cpp"
   "src/matrix_elements/*.cpp"
   "src/matrix_elements/*/*.cpp"
)

add_executable(RSymSQCD ${SOURCES})

target_link_libraries(
   RSymSQCD
   cuba
   Boost::program_options
   ooptools
   gfortran
   LHAPDF
   ColorFull
   Eigen3::Eigen
   nlohmann_json::nlohmann_json
)
