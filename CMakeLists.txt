cmake_minimum_required (VERSION 3.8)
project (RSymSQCD)

set (CMAKE_CXX_STANDARD 17)

#add_subdirectory(tests)

if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
   link_directories(/opt/local/lib)
   link_directories(${Cuba}/lib ${LT}/lib)
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
   link_directories(${Cuba}/lib64 ${LT}/lib64)
else()
   message("Couldn't recognize system (Darwin or Linux). Continuing hoping for the best")
endif()

include_directories (
   include 
   /Users/Navir/Fizyka/programy/gcc/boost-1.65.0/include ${LHAPDF6}/include 
   ${COLORFULL}/include
   ${Cuba}
   ${LT}/include
   ${rk}/include
)

link_directories(
   /Users/Navir/Fizyka/programy/gcc/boost-1.65.0/lib
   ${COLORFULL}/lib/ColorFull
   ${LHAPDF6}/lib
   /Users/Navir/Fizyka/Programowanie/rust/catani-seymour-dipole-subtraction/target/debug
   /Users/Navir/Fizyka/Programowanie/rust/catani-seymour-dipole-subtraction/target/debug/deps
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
   boost_program_options
   ooptools
   gfortran
   LHAPDF
   ColorFull
   cs_dipole_subtraction
)
