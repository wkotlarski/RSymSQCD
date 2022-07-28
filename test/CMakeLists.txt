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

list(
   APPEND LO_RUN_CARDS
   pp_sqLsqR/pp_sqLsqR_BMP1_LO.ini pp_sqLsqR/pp_sqLsqR_BMP2_LO.ini pp_sqLsqR/pp_sqLsqR_BMP3_LO.ini
   pp_sqsqdagger/pp_sqsqdagger_BMP1_LO.ini pp_sqsqdagger/pp_sqsqdagger_BMP2_LO.ini pp_sqsqdagger/pp_sqsqdagger_BMP3_LO.ini
)
# σLO(LO PDF) column from Tab. 3 of 1707.04557
list(
   APPEND LO_RESULTS
   pp_sqLsqR_BMP1_LO.json pp_sqLsqR_BMP2_LO.json pp_sqLsqR_BMP3_LO.json
   pp_sqsqdagger_BMP1_LO.json pp_sqsqdagger_BMP2_LO.json pp_sqsqdagger_BMP3_LO.json
)
foreach(el IN ZIP_LISTS LO_RUN_CARDS LO_RESULTS)
add_test(
   NAME
      "LO cross-section for ${el_0}"
   COMMAND ${RUBY_EXECUTABLE} ${CMAKE_SOURCE_DIR}/test/compare_json_outputs.rb ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/run_cards/${el_0} ${CMAKE_SOURCE_DIR}/test/ref_outputs/${el_1}
)
endforeach()

list(APPEND NLO_RUN_CARDS pp_suLsuR/pp_suLsuR_BMP1.ini)
list(APPEND NLO_RESULTS pp_suLsuR_BMP1.json)
foreach(el IN ZIP_LISTS NLO_RUN_CARDS NLO_RESULTS)
add_test(
   NAME
      "NLO cross-section for ${el_0}"
   COMMAND ${RUBY_EXECUTABLE} ${CMAKE_SOURCE_DIR}/test/compare_json_outputs.rb ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/run_cards/${el_0} ${CMAKE_SOURCE_DIR}/test/ref_outputs/${el_1}
)
endforeach()
