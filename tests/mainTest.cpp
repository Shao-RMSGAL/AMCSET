/**
 * @file mainTest.cpp
 * @brief Unit tests for the main function.
 */

#ifndef MAINTEST_CXX
#define MAINTEST_CXX

#include <gtest/gtest.h>
#include <sstream>

#include "eSRIM.h" // Include the header file containing your main function

// This suite of tests ensures that the main function behaves as expected


// This is a sample test that should always pass
TEST(MainTest, OneEqualsZero) {
    EXPECT_EQ(0, 0);
}

// No arguments test
TEST(MainFunctionTest, RunProgramWithNoArguments) {
    const char* argv[] = {"eSRIM"};
    int argc = sizeof(argv) / sizeof(argv[0]);

    EXPECT_EQ(startESRIM(argc, const_cast<char**>(argv)), 0);
}

// Test for providing a custom filename
TEST(MainFunctionTest, RunProgramWithCustomFilename) {
    const char* argv[] = {"eSRIM", "-f", "custom_settings.txt"};
    int argc = sizeof(argv) / sizeof(argv[0]);

    std::string input_string = "y";
    std::istringstream stringStream(input_string); 

    EXPECT_EQ(startESRIM(argc, const_cast<char**>(argv), stringStream), 0);
}

// Test for recording execution time
TEST(MainFunctionTest, RunProgramWithTimeOption) {
    const char* argv[] = {"eSRIM", "-t"};
    int argc = sizeof(argv) / sizeof(argv[0]);

    EXPECT_EQ(startESRIM(argc, const_cast<char**>(argv)), 0);
}

// Test for displaying example settings file
TEST(MainFunctionTest, RunProgramWithSettingsOption) {
    const char* argv[] = {"eSRIM", "-s"};
    int argc = sizeof(argv) / sizeof(argv[0]);

    EXPECT_EQ(startESRIM(argc, const_cast<char**>(argv)), 0);
}

// Test for displaying active settings
TEST(MainFunctionTest, RunProgramWithDisplayOption) {
    const char* argv[] = {"eSRIM", "-d"};
    int argc = sizeof(argv) / sizeof(argv[0]);

    EXPECT_EQ(startESRIM(argc, const_cast<char**>(argv)), 0);
}

// Test for displaying help message
TEST(MainFunctionTest, RunProgramWithHelpOption) {
    const char* argv[] = {"eSRIM", "-h"};
    int argc = sizeof(argv) / sizeof(argv[0]);

    EXPECT_EQ(startESRIM(argc, const_cast<char**>(argv)), 0);
}

// Test for displaying progress of the simulation
TEST(MainFunctionTest, RunProgramWithProgressOption) {
    const char* argv[] = {"eSRIM", "-p"};
    int argc = sizeof(argv) / sizeof(argv[0]);

    EXPECT_EQ(startESRIM(argc, const_cast<char**>(argv)), 0);
}

// Test for combining multiple arguments
TEST(MainFunctionTest, RunProgramWithMultipleArguments) {
    const char* argv[] = {"eSRIM", "-f", "custom_settings.txt", "-t", "-d"};
    int argc = sizeof(argv) / sizeof(argv[0]);
    
    const std::string input_string = "y";
    std::stringstream stringStream(input_string); 

    EXPECT_EQ(startESRIM(argc, const_cast<char**>(argv), stringStream), 0);
}

// Test for unusual/unrecognized/malformed arguments
TEST(MainFunctionTest, RunProgramWithUnrecognizedArguments) {
    const char* argv[] = {"eSRIM", "-x", "-y", "-z"};
    int argc = sizeof(argv) / sizeof(argv[0]);

    EXPECT_EQ(startESRIM(argc, const_cast<char**>(argv)), 0);
    ASSERT_EQ(0,1);
}

#endif
