cmake_minimum_required (VERSION 3.8)

project (RSymSQCD CXX)

set (CMAKE_CXX_STANDARD 17)

file (GLOB SOURCES
   "src/*.cpp"
   "src/matrix_elements/*.cpp"
   "src/matrix_elements/*/*.cpp"
)


include_directories (
   include 
   ${LHAPDF6}/include 
   ${COLORFULL}/include
   ${Cuba}
   ${LT}/include
   ${rk}/include
)

link_directories(
   ${COLORFULL}/lib/ColorFull
   ${Cuba}/lib64
   ${LT}/lib64
   ${LHAPDF6}/lib
   $ENV{HOME}/Dokumenty/IdeaProjects/cs_dipole_subtraction/target/debug
   $ENV{HOME}/Dokumenty/IdeaProjects/cs_dipole_subtraction/target/debug/deps
)

# why does it have to be after above include and link?
add_executable(RSymSQCD ${SOURCES})

# locate boost instalation
FIND_PACKAGE(Boost)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
    target_link_libraries(RSymSQCD boost_program_options)
endif()

target_link_libraries(
   RSymSQCD
   cuba
   ooptools
   gfortran
   LHAPDF
   ColorFull
   cs_dipole_subtraction
)

