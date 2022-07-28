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

list(APPEND BORN_RUN_CARDS pp_sqLsqR/pp_sqLsqR_BMP1_LO.ini pp_sqLsqR/pp_sqLsqR_BMP2_LO.ini pp_sqLsqR/pp_sqLsqR_BMP3_LO.ini)
# σLO(LO PDF) column from Tab. 3 of 1707.04557
list(APPEND BORN_RESULTS 13.065515658699582 4.608168157871303 159.82254549998856)
foreach(el IN ZIP_LISTS BORN_RUN_CARDS BORN_RESULTS)
add_test(
   NAME
      "σLO for ${el_0}"
   COMMAND ${RUBY_EXECUTABLE} ${CMAKE_SOURCE_DIR}/test/check_born_xsec.rb ${CMAKE_SOURCE_DIR}/run_cards/${el_0} ${el_1} ${CMAKE_BINARY_DIR}
)
endforeach()

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

add_test(
   NAME
      "pp_suLsuLdagger total xsec in BMP1"
   COMMAND
      RSymSQCD --card ${CMAKE_SOURCE_DIR}/run_cards/pp_suLsuLdagger/pp_suLsuLdagger_BMP1.ini --json-outputfile-name "pp_suLsuLdagger_BMP1.out" --precision-hard=4 --precision-virt=4
)
add_test(NAME "compare outputs for pp_suLsuL* in BMP1" COMMAND ${RUBY_EXECUTABLE} "${CMAKE_SOURCE_DIR}/test/compare_json_outputs.rb" "${CMAKE_SOURCE_DIR}/test/ref_outputs/pp_suLsuLdagger_BMP1.ref" "pp_suLsuLdagger_BMP1.out")

add_test(
   NAME
      "pp_suLsuLdagger total xsec in BMP2 with DR"
   COMMAND
      RSymSQCD --card ${CMAKE_SOURCE_DIR}/run_cards/pp_suLsuLdagger/pp_suLsuLdagger_BMP2_DR.ini --json-outputfile-name "pp_suLsuLdagger_BMP2_DR.out" --precision-hard=4 --precision-virt=4
)
add_test(NAME "compare outputs for pp_suLsuL* in BMP2 with DR" COMMAND ${RUBY_EXECUTABLE} "${CMAKE_SOURCE_DIR}/test/compare_json_outputs.rb" "${CMAKE_SOURCE_DIR}/test/ref_outputs/pp_suLsuLdagger_BMP2_DR.ref" "pp_suLsuLdagger_BMP2_DR.out")

add_test(
   NAME
      "pp_suLsuLdagger total xsec in BMP3 with DR"
   COMMAND
      RSymSQCD --card ${CMAKE_SOURCE_DIR}/run_cards/pp_suLsuLdagger/pp_suLsuLdagger_BMP3_DR.ini --json-outputfile-name "pp_suLsuLdagger_BMP3_DR.out" --precision-hard=4 --precision-virt=4
)
add_test(NAME "compare outputs for pp_suLsuL* in BMP3 with DR" COMMAND ${RUBY_EXECUTABLE} "${CMAKE_SOURCE_DIR}/test/compare_json_outputs.rb" "${CMAKE_SOURCE_DIR}/test/ref_outputs/pp_suLsuLdagger_BMP3_DR.ref" "pp_suLsuLdagger_BMP3_DR.out")

