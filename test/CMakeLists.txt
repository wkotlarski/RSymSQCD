message(STATUS "Enabling tests")

find_package(GTest REQUIRED)

add_executable(foo foo.cpp)
target_link_libraries(foo GTest::gtest GTest::gtest_main)

add_test(AllTestsInFoo foo)
