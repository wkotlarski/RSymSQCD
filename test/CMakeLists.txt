message(STATUS "Enabling tests")

find_package(GTest REQUIRED)

add_executable(tree_matrix_elements tree_matrix_elements.cpp)
target_link_libraries(tree_matrix_elements GTest::gtest GTest::gtest_main RSymSQCDlib)

add_test(TreeMatrixElements tree_matrix_elements)

add_executable(hard_matrix_elements hard_matrix_elements.cpp)
target_link_libraries(hard_matrix_elements GTest::gtest GTest::gtest_main RSymSQCDlib)
add_test(HardMatrixElements hard_matrix_elements)

add_test(NAME RSymSQCD COMMAND RSymSQCD --card ${CMAKE_SOURCE_DIR}/run_cards/pp_suLsuR/pp_suLsuR_BMP1.ini --json-outputfile-name "pp_suLsuR_BMP1.out" --precision-hard=4 --precision-sc=7 --precision-virt=5)
add_test(compare_pp_suLsuR_BMP1_outputs ${CMAKE_COMMAND} -E compare_files pp_suLsuR_BMP1.out ${CMAKE_SOURCE_DIR}/test/ref_outputs/pp_suLsuR_BMP1.ref)
