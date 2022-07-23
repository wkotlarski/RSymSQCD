message(STATUS "Enabling tests")

find_package(GTest REQUIRED)

add_executable(tree_matrix_elements tree_matrix_elements.cpp)
target_link_libraries(tree_matrix_elements GTest::gtest GTest::gtest_main RSymSQCDlib)

add_test("test tree-level matrix elements" tree_matrix_elements)

add_executable(hard_matrix_elements hard_matrix_elements.cpp)
target_link_libraries(hard_matrix_elements GTest::gtest GTest::gtest_main RSymSQCDlib)
add_test("test hard and non-colinear matrix elements" hard_matrix_elements)

find_package(Ruby)

# allow whitespaces in test names
cmake_policy(SET CMP0110 NEW)

add_test(
   NAME
      "pp_suLsuR total xsec in BMP1"
   COMMAND
      RSymSQCD --card ${CMAKE_SOURCE_DIR}/run_cards/pp_suLsuR/pp_suLsuR_BMP1.ini --json-outputfile-name "pp_suLsuR_BMP1.out" --precision-hard=4 --precision-virt=4
)
add_test(NAME "compare outputs for pp_suLsuR in BMP1" COMMAND ${RUBY_EXECUTABLE} "${CMAKE_SOURCE_DIR}/test/compare_json_outputs.rb" "${CMAKE_SOURCE_DIR}/test/ref_outputs/pp_suLsuR_BMP1.ref" "pp_suLsuR_BMP1.out")

add_test(
   NAME
      "pp_suLsuR total xsec in BMP2 with DR"
   COMMAND
      RSymSQCD --card ${CMAKE_SOURCE_DIR}/run_cards/pp_suLsuR/pp_suLsuR_BMP2_DR.ini --json-outputfile-name "pp_suLsuR_BMP2_DR.out" --precision-hard=4 --precision-virt=4
)
add_test(NAME "compare outputs for pp_suLsuR in BMP2 with DR" COMMAND ${RUBY_EXECUTABLE} "${CMAKE_SOURCE_DIR}/test/compare_json_outputs.rb" "${CMAKE_SOURCE_DIR}/test/ref_outputs/pp_suLsuR_BMP2_DR.ref" "pp_suLsuR_BMP2_DR.out")

add_test(
   NAME
      "pp_suLsuR total xsec in BMP3 with DR"
   COMMAND
      RSymSQCD --card ${CMAKE_SOURCE_DIR}/run_cards/pp_suLsuR/pp_suLsuR_BMP3_DR.ini --json-outputfile-name "pp_suLsuR_BMP3_DR.out" --precision-hard=4 --precision-virt=4
)
add_test(NAME "compare outputs for pp_suLsuR in BMP3 with DR" COMMAND ${RUBY_EXECUTABLE} "${CMAKE_SOURCE_DIR}/test/compare_json_outputs.rb" "${CMAKE_SOURCE_DIR}/test/ref_outputs/pp_suLsuR_BMP3_DR.ref" "pp_suLsuR_BMP3_DR.out")
