enable_testing()

add_executable( common-test
    ${PROJECT_SOURCE_DIR}/test/common/common_test.cpp
    )

target_link_libraries(common-test PRIVATE
    GTest::gtest_main
    amcset-common
    )

target_include_directories( common-test PRIVATE
    ${PROJECT_SOURCE_DIR}/common/include
    )

include(GoogleTest)
gtest_discover_tests(common-test)
