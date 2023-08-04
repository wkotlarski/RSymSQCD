message(STATUS "Enabling tests")

find_package(GTest REQUIRED)

# allow whitespaces in test names
cmake_policy(SET CMP0110 NEW)

find_program(Ruby REQUIRED NAMES ruby)

# test matrix elements
# ------------------------------------------------------------------------------------------------------

add_executable(tree_matrix_elements tree_matrix_elements.cpp)
target_include_directories(tree_matrix_elements PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(tree_matrix_elements GTest::gtest GTest::gtest_main MRSSM)
add_test("tree-level matrix elements" tree_matrix_elements)

add_executable(MRSSM_hard_matrix_elements MRSSM_hard_matrix_elements.cpp)
target_include_directories(MRSSM_hard_matrix_elements PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(MRSSM_hard_matrix_elements GTest::gtest GTest::gtest_main MRSSM)
add_test("MRSSM hard and non-colinear matrix elements" MRSSM_hard_matrix_elements)

add_executable(MRSSM_virtual_matrix_elements MRSSM_virtual_matrix_elements.cpp)
target_include_directories(MRSSM_virtual_matrix_elements PUBLIC ${PROJECT_SOURCE_DIR}/include ${LIBLOOPTOOLS_INCLUDE_DIR})
target_link_libraries(MRSSM_virtual_matrix_elements GTest::gtest GTest::gtest_main MRSSM ${FORTRAN_LIBS} ${LIBLOOPTOOLS_LIBRARY})
add_test("MRSSM virtual matrix elements" MRSSM_virtual_matrix_elements)

add_executable(Sgluons_hard_matrix_elements Sgluons_hard_matrix_elements.cpp)
target_include_directories(Sgluons_hard_matrix_elements PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(Sgluons_hard_matrix_elements GTest::gtest GTest::gtest_main Sgluons)
add_test("Sgluons hard and non-colinear matrix elements" Sgluons_hard_matrix_elements)

add_executable(MRSSM_tree_xsections MRSSM_tree_xsections.cpp)
target_include_directories(MRSSM_tree_xsections PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(MRSSM_tree_xsections GTest::gtest GTest::gtest_main MRSSM)
add_test("MRSSM tree xsections" MRSSM_tree_xsections)

# test LO cross sections
# ------------------------------------------------------------------------------------------------------

list(
   APPEND LO_RUN_CARDS
   pp_sqLsqR/pp_sqLsqR_BMP1_LO pp_sqLsqR/pp_sqLsqR_BMP2_LO pp_sqLsqR/pp_sqLsqR_BMP3_LO
   pp_sqLsqR+cc/pp_sqLsqR+cc_BMP1_LO pp_sqLsqR+cc/pp_sqLsqR+cc_BMP2_LO pp_sqLsqR+cc/pp_sqLsqR+cc_BMP3_LO
   pp_sqsqdagger/pp_sqsqdagger_BMP1_LO pp_sqsqdagger/pp_sqsqdagger_BMP2_LO pp_sqsqdagger/pp_sqsqdagger_BMP3_LO
   pp_glglbar/pp_glglbar_BMP1_LO pp_glglbar/pp_glglbar_BMP2_LO pp_glglbar/pp_glglbar_BMP3_LO
)
# checked against σLO(LO PDF) column from Tab. 3 of 1707.04557 (pp_sqLsqR and pp_sqsqdagger)
# and against Sebastian Liebschner's private code (pp_sqLsqR+cc)
list(
   APPEND LO_RESULTS
   pp_sqLsqR_BMP1_LO pp_sqLsqR_BMP2_LO pp_sqLsqR_BMP3_LO
   pp_sqLsqR+cc_BMP1_LO pp_sqLsqR+cc_BMP2_LO pp_sqLsqR+cc_BMP3_LO
   pp_sqsqdagger_BMP1_LO pp_sqsqdagger_BMP2_LO pp_sqsqdagger_BMP3_LO
   pp_glglbar_BMP1_LO pp_glglbar_BMP2_LO pp_glglbar_BMP3_LO
)
foreach(el IN ZIP_LISTS LO_RUN_CARDS LO_RESULTS)
add_test(
   NAME
      "LO cross-section for ${el_0}"
   COMMAND
      ${RUBY_EXECUTABLE} ${CMAKE_SOURCE_DIR}/test/compare_json_outputs.rb
         $<TARGET_FILE:RSymSQCD> ${CMAKE_SOURCE_DIR}/run_cards/${el_0}.ini ${CMAKE_SOURCE_DIR}/test/ref_outputs/LO/${el_1}.json
)
endforeach()

# test NLO cross sections
# ------------------------------------------------------------------------------------------------------

# we test 16 cross-sections
# 3 BMPs x 5 processes ([suLsuR, suLsuLdagger, sqLsqR, sqLsqR+cc, sqsqdagger])
#    + suLsuR for BMP1 with alternative scale choice (muR != muF)
list(
   APPEND NLO_RUN_CARDS
   # pp -> suL suR
   pp_suLsuR/pp_suLsuR_BMP1 pp_suLsuR/pp_suLsuR_BMP2_DR pp_suLsuR/pp_suLsuR_BMP3_DR
   # pp -> suL suR with mu_r != mu_F
   pp_suLsuR/pp_suLsuR_BMP1-scale_variation
   # pp -> suL suL*
   pp_suLsuLdagger/pp_suLsuLdagger_BMP1 pp_suLsuLdagger/pp_suLsuLdagger_BMP2_DR pp_suLsuLdagger/pp_suLsuLdagger_BMP3_DR
   # pp -> sqL sqR
   pp_sqLsqR/pp_sqLsqR_BMP1 pp_sqLsqR/pp_sqLsqR_BMP2_DR pp_sqLsqR/pp_sqLsqR_BMP3_DR
   # pp -> sqL sqR + sqL* sqR*
   pp_sqLsqR+cc/pp_sqLsqR+cc_BMP1 pp_sqLsqR+cc/pp_sqLsqR+cc_BMP2_DR pp_sqLsqR+cc/pp_sqLsqR+cc_BMP3_DR
   # pp -> sqL sq'L* + sqR sq'R*
   pp_sqsqdagger/pp_sqsqdagger_BMP1 pp_sqsqdagger/pp_sqsqdagger_BMP2_DR pp_sqsqdagger/pp_sqsqdagger_BMP3_DR
)
list(
   APPEND NLO_RESULTS
   # pp -> suL suR
   pp_suLsuR_BMP1 pp_suLsuR_BMP2_DR pp_suLsuR_BMP3_DR
   # pp -> suL suR with mu_r != mu_F
   pp_suLsuR_BMP1-scale_variation
   # pp -> suL suL*
   pp_suLsuLdagger_BMP1 pp_suLsuLdagger_BMP2_DR pp_suLsuLdagger_BMP3_DR
   # pp -> sqL sqR
   pp_sqLsqR_BMP1 pp_sqLsqR_BMP2_DR pp_sqLsqR_BMP3_DR
   # pp -> sqL sqR + sqL* sqR*
   pp_sqLsqR+cc_BMP1 pp_sqLsqR+cc_BMP2_DR pp_sqLsqR+cc_BMP3_DR
   # pp -> sqL sq'L* + sqR sq'R*
   pp_sqsqdagger_BMP1 pp_sqsqdagger_BMP2_DR pp_sqsqdagger_BMP3_DR
)
foreach(el IN ZIP_LISTS NLO_RUN_CARDS NLO_RESULTS)
add_test(
   NAME
      "NLO cross-section for ${el_0}"
      COMMAND
         ${RUBY_EXECUTABLE} ${CMAKE_SOURCE_DIR}/test/compare_json_outputs.rb
            $<TARGET_FILE:RSymSQCD> ${CMAKE_SOURCE_DIR}/run_cards/${el_0}.ini ${CMAKE_SOURCE_DIR}/test/ref_outputs/NLO/${el_1}.json
)
endforeach()
