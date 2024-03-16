#ifndef MAINTEST_CXX
#define MAINTEST_CXX

#include <gtest/gtest.h>
#include "main.h" // Include the header file containing your main function

// Buggy, causes linker error
// TEST(MainFunctionTest, NormalExit) {
//     // Call your main function here
//     // Replace "YourMainFunction" with the actual name of your main function
//     int argc = 1;  // or appropriate argc value
//     char* argv[] = { (char*)"build/bin/eSRIM", nullptr };  // or appropriate argv values
//     EXPECT_EQ(startESRIM(argc, argv), 0);  // Assuming 0 indicates normal exit
// }

// Trivial test for demonstration
TEST(TrivialTest, ZeroEqualsZero) {
    EXPECT_EQ(0, 0); 
}

TEST(TrivialTest, OneEqualsZero) {
    // This test is expected to fail
    EXPECT_EQ(1, 0); 
}

#endif  