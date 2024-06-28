enable_testing()

add_executable( common-test
    ${PROJECT_SOURCE_DIR}/test/common/common_test.cpp
    )

target_link_libraries(
    common-test
    GTest::gtest_main
    )

include(GoogleTest)
gtest_discover_tests(common-test)
