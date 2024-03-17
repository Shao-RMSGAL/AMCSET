/**
 * @file mainFunctionTest.cpp
 * @brief Unit tests for the main function.
 */

#ifndef MAINFUNCTIONTEST_CXX
#define MAINFUNCTIONTEST_CXX

#include <gtest/gtest.h>
#include <sstream>
#include <filesystem>

#include "eSRIM.h" // Include the header file containing your main function

// Refer to https://google.github.io/googletest/
// for proper test-case writing practices.

// This suite of tests ensures that the main function behaves as expected

// Write a class MainFunctionTest that inherits from testing::Test
class MainFunctionTest : public testing::Test {
    protected:
    

        // This function will be called before each test
        static void SetUpTestSuite() {
            // Delete the settings file if it exists and confirm that it was deleted
            
            
            if(std::filesystem::exists("settings.txt")) {
                std::filesystem::remove("settings.txt");
            }

            // Confirm that there is no file called "non-existing-settings.txt". If there is, delete it
            if(std::filesystem::exists("non-existing-settings.txt")) {
                std::filesystem::remove("non-existing-settings.txt");
            }

            const char* argv[] = {"eSRIM", nullptr};
            const int argc = sizeof(argv) / sizeof(argv[0]) - 1;

            std::istringstream stdIn("y");
            std::ostringstream stdOut;
            std::ostringstream stdErr;
            int result = startESRIM(argc, argv, stdIn, stdOut, stdErr);

            // There should be output to standard error or standard output
            EXPECT_FALSE(stdOut.str().empty());
            EXPECT_FALSE(stdErr.str().empty());
            EXPECT_FALSE(stdIn.good() && stdIn.tellg() == 0);
            ASSERT_EQ(result, 0);
        }

        // This function will be called after each test
        static void TearDownTestSuite() {
            // Reset contents of the settings file
            // SetUpTestSuite();

            // Delete all the files in the /output directory.
            for(const auto& entry : std::filesystem::directory_iterator("output")) {
                std::filesystem::remove(entry.path());
            }
        }

        void SetUp() override {

        }

        void TearDown() override {

        }
};

// This is a sample test that should always pass
TEST_F(MainFunctionTest, OneEqualsZero) {
    ASSERT_NE(1, 0);
}

// No arguments test
TEST_F(MainFunctionTest, RunProgramWithNoArguments) {
    const char* argv[] = {"eSRIM", nullptr};
    const int argc = sizeof(argv) / sizeof(argv[0]) - 1;

    std::istringstream stdIn("y");
    std::ostringstream stdOut;
    std::ostringstream stdErr;
    int result = startESRIM(argc, argv, stdIn, stdOut, stdErr);

    // There should be no output to standard error or standard output
    EXPECT_TRUE(stdOut.str().empty());
    EXPECT_TRUE(stdErr.str().empty());
    EXPECT_TRUE(stdIn.good() && stdIn.tellg() == 0);
    ASSERT_EQ(result, 0);
}

// // Test an unrecognized flag.
// TEST_F(MainFunctionTest, RunProgramWithUnrecognizedFlag) {
//     const char* argv[] = {"eSRIM", "-xxxx", nullptr};
//     const int argc = sizeof(argv) / sizeof(argv[0]) - 1;

//     std::istringstream stdIn;
//     std::ostringstream stdOut;
//     std::ostringstream stdErr;
//     int result = startESRIM(argc, argv, stdIn, stdOut, stdErr);

//     // There should be output to standard error
//     EXPECT_TRUE(stdOut.str().empty());
//     EXPECT_TRUE(stdErr.str().empty());
//     EXPECT_TRUE(stdIn.good() && stdIn.tellg() == 0);
//     // ASSERT_NE(result, 0);
// }

// Test the -h flag
// TEST_F(MainFunctionTest, RunProgramWithHelpFlag) {
//     const char* argv[] = {"eSRIM", "-h", nullptr};
//     const int argc = sizeof(argv) / sizeof(argv[0]) - 1;

//     std::istringstream stdIn;
//     std::ostringstream stdOut;
//     std::ostringstream stdErr;
//     int result = startESRIM(argc, argv, stdIn, stdOut, stdErr);

//     // There should be output to standard output
//     std::cout << "Standard output:\n" 
//               << stdOut.str()
//               << std::endl;
//     std::cout << "Standard error:\n"
//               << stdErr.str()
//               << std::endl;
//     EXPECT_FALSE(stdOut.str().empty());
//     EXPECT_TRUE(stdErr.str().empty());
//     EXPECT_TRUE(stdIn.good() && stdIn.tellg() == 0);
//     ASSERT_EQ(result, 0);
// }

// // Test the --help flag
// TEST_F(MainFunctionTest, RunProgramWithDoubleDashHelpFlag) {
//     const char* argv[] = {"eSRIM", "--help"};
//     const int argc = sizeof(argv) / sizeof(argv[0]);

//     std::istringstream stdIn;
//     std::ostringstream stdOut;
//     std::ostringstream stdErr;
//     int result = startESRIM(argc, const_cast<char**>(argv), stdIn, stdOut, stdErr);

//     // There should be output to standard output
//     std::cout << "Standard output:\n" 
//               << stdOut.str()
//               << std::endl;
//     std::cout << "Standard error:\n"
//               << stdErr.str()
//               << std::endl;
//     EXPECT_FALSE(stdOut.str().empty());
//     EXPECT_TRUE(stdErr.str().empty());
//     EXPECT_TRUE(stdIn.good() && stdIn.tellg() == 0);
//     ASSERT_EQ(result, 0);
// }

// // // Test the -f flag with an existing settings file passed
// TEST_F(MainFunctionTest, RunProgramWithExistingSettingsFile) {

//     const char* argv[] = {"eSRIM", "-f", "settings.txt"};
//     const int argc = sizeof(argv) / sizeof(argv[0]);

//     std::istringstream stdIn;
//     std::ostringstream stdOut;
//     std::ostringstream stdErr;
//     int result = startESRIM(argc, const_cast<char**>(argv), stdIn, stdOut, stdErr);

//     // There should be no output to standard error or standard output
//     std::cout << "Standard output:\n" 
//               << stdOut.str()
//               << std::endl;
//     std::cout << "Standard error:\n"
//               << stdErr.str()
//               << std::endl;
//     EXPECT_TRUE(stdOut.str().empty());
//     EXPECT_TRUE(stdErr.str().empty());
//     EXPECT_TRUE(stdIn.good() && stdIn.tellg() == 0);
//     ASSERT_EQ(result, 0);
// }

// // Test the -f flag with a non-existing settings file passed and "n" from the standard input
// TEST_F(MainFunctionTest, RunProgramWithNonExistingSettingsFile) {
//     const char* argv[] = {"eSRIM", "-f", "non-existing-settings.txt"};
//     const int argc = sizeof(argv) / sizeof(argv[0]);

//     std::istringstream stdIn("n\n");
//     std::ostringstream stdOut;
//     std::ostringstream stdErr;
//     int result = startESRIM(argc, const_cast<char**>(argv), stdIn, stdOut, stdErr);

//     // There should be output to standard error and standard output
//     std::cout << "Standard output:\n" 
//               << stdOut.str()
//               << std::endl;
//     std::cout << "Standard error:\n"
//               << stdErr.str()
//               << std::endl;
//     EXPECT_FALSE(stdOut.str().empty());
//     EXPECT_FALSE(stdErr.str().empty());
//     EXPECT_FALSE(stdIn.good() && stdIn.tellg() == 1);
//     ASSERT_EQ(result, 0);
// }

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

#endif
