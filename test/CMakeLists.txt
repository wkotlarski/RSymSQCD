message(STATUS "Enabling tests")

find_package(GTest REQUIRED)

# allow whitespaces in test names
cmake_policy(SET CMP0110 NEW)

find_package(Ruby)

# test matrix elements
# ------------------------------------------------------------------------------------------------------

add_executable(tree_matrix_elements tree_matrix_elements.cpp)
target_link_libraries(tree_matrix_elements GTest::gtest GTest::gtest_main RSymSQCDlib)
add_test("tree-level matrix elements" tree_matrix_elements)

add_executable(MRSSM_hard_matrix_elements MRSSM_hard_matrix_elements.cpp)
target_link_libraries(MRSSM_hard_matrix_elements GTest::gtest GTest::gtest_main RSymSQCDlib)
add_test("hard and non-colinear matrix elements" MRSSM_hard_matrix_elements)

add_executable(Sgluons_hard_matrix_elements Sgluons_hard_matrix_elements.cpp)
target_link_libraries(Sgluons_hard_matrix_elements GTest::gtest GTest::gtest_main RSymSQCDlib)
add_test("hard and non-colinear matrix elements" Sgluons_hard_matrix_elements)

# test LO cross sections
# ------------------------------------------------------------------------------------------------------

list(
   APPEND LO_RUN_CARDS
   pp_sqLsqR/pp_sqLsqR_BMP1_LO.ini pp_sqLsqR/pp_sqLsqR_BMP2_LO.ini pp_sqLsqR/pp_sqLsqR_BMP3_LO.ini
   pp_sqLsqR+cc/pp_sqLsqR+cc_BMP1_LO.ini pp_sqLsqR+cc/pp_sqLsqR+cc_BMP2_LO.ini pp_sqLsqR+cc/pp_sqLsqR+cc_BMP3_LO.ini
   pp_sqsqdagger/pp_sqsqdagger_BMP1_LO.ini pp_sqsqdagger/pp_sqsqdagger_BMP2_LO.ini pp_sqsqdagger/pp_sqsqdagger_BMP3_LO.ini
   pp_glglbar/pp_glglbar_BMP1_LO.ini pp_glglbar/pp_glglbar_BMP2_LO.ini pp_glglbar/pp_glglbar_BMP3_LO.ini
)
# checked against σLO(LO PDF) column from Tab. 3 of 1707.04557 (pp_sqLsqR and pp_sqsqdagger)
# and against Sebastian Liebschner's private code (pp_sqLsqR+cc)
list(
   APPEND LO_RESULTS
   pp_sqLsqR_BMP1_LO.json pp_sqLsqR_BMP2_LO.json pp_sqLsqR_BMP3_LO.json
   pp_sqLsqR+cc_BMP1_LO.json pp_sqLsqR+cc_BMP2_LO.json pp_sqLsqR+cc_BMP3_LO.json
   pp_sqsqdagger_BMP1_LO.json pp_sqsqdagger_BMP2_LO.json pp_sqsqdagger_BMP3_LO.json
   pp_glglbar_BMP1_LO.json pp_glglbar_BMP2_LO.json pp_glglbar_BMP3_LO.json
)
foreach(el IN ZIP_LISTS LO_RUN_CARDS LO_RESULTS)
add_test(
   NAME
      "LO cross-section for ${el_0}"
   COMMAND
      ${RUBY_EXECUTABLE} ${CMAKE_SOURCE_DIR}/test/compare_json_outputs.rb
         ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/run_cards/${el_0} ${CMAKE_SOURCE_DIR}/test/ref_outputs/LO/${el_1} 3
)
endforeach()

# test NLO cross sections
# ------------------------------------------------------------------------------------------------------

list(
   APPEND NLO_RUN_CARDS
   # pp -> suL suR
   pp_suLsuR/pp_suLsuR_BMP1.ini pp_suLsuR/pp_suLsuR_BMP2_DR.ini pp_suLsuR/pp_suLsuR_BMP3_DR.ini
   # pp -> suL suL*
   pp_suLsuLdagger/pp_suLsuLdagger_BMP1.ini pp_suLsuLdagger/pp_suLsuLdagger_BMP2_DR.ini
   # pp -> sqL sqR
   pp_sqLsqR/pp_sqLsqR_BMP1.ini pp_sqLsqR/pp_sqLsqR_BMP2_DR.ini pp_sqLsqR/pp_sqLsqR_BMP3_DR.ini
   # pp -> sqL sqR + sqL* sqR*
   pp_sqLsqR+cc/pp_sqLsqR+cc_BMP1.json pp_sqLsqR+cc/pp_sqLsqR+cc_BMP2_DR.json pp_sqLsqR+cc/pp_sqLsqR+cc_BMP3_DR.json
)
list(
   APPEND NLO_RESULTS
   pp_suLsuR_BMP1.json pp_suLsuR_BMP2_DR.json pp_suLsuR_BMP3_DR.json
   pp_suLsuLdagger_BMP1.json pp_suLsuLdagger_BMP2_DR.json
   pp_sqLsqR_BMP1.json pp_sqLsqR_BMP2_DR.jso pp_sqLsqR_BMP3_DR.json
   pp_sqLsqR+cc_BMP1.json pp_sqLsqR+cc_BMP2_DR.json pp_sqLsqR+cc_BMP3_DR.json
)
foreach(el IN ZIP_LISTS NLO_RUN_CARDS NLO_RESULTS)
add_test(
   NAME
      "NLO cross-section for ${el_0}"
      COMMAND
         ${RUBY_EXECUTABLE} ${CMAKE_SOURCE_DIR}/test/compare_json_outputs.rb
            ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/run_cards/${el_0} ${CMAKE_SOURCE_DIR}/test/ref_outputs/NLO/${el_1} 3
)
endforeach()
