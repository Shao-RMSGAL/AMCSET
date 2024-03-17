/******************************************************************************
 * Filename:        main.cpp
 * Project:         eSRIM
 * Author:          Nathaniel Thomas
 * Date Created:    March 15, 2024
 * Date Modified:   March 15, 2024
 * File Version:    1.0
 * Group:           Dr. Shao's RMSLCF Group
 *
 * Description:
 * This is the main file for the eSRIM code developed by Dr. Shao and converted
 * to the C++ programming language
 *
 * Revision History
 * 1.0:
 * - File created to account for testing requirements
 *****************************************************************************/
#include "eSRIM.h"

// Main program  start
int main(const int argc, const char* argv[]) {
    return startESRIM(argc, argv);
}

// TODO: TEMPORARY DEBUGGING BELOW. DELETE WHEN DONE.

// #ifdef DEBUG_MODE
// #define DEBUG_PRINT(x) std::cout << x << std::endl;
// #else
// #define DEBUG_PRINT(x)
// #endif
// #include <sstream>
// #include <filesystem>

// // Main program  start
// int main() {
//     // Delete the "settings.txt" file if it exists
//     if(std::filesystem::exists("settings.txt")) {
//         std::filesystem::remove("settings.txt");
//     }
//     const char* argv[] = {"eSRIM", nullptr};
//     std::istringstream stdIn("y\n");
//     std::ostringstream stdOut;
//     std::ostringstream stdErr;
//     const int argc = sizeof(argv) / sizeof(argv[0]) - 1;
//     int result=  startESRIM(argc, argv, stdIn, stdOut, stdErr);
//     DEBUG_PRINT("Standard output:\n" << stdOut.str() << std::endl);
//     DEBUG_PRINT("Standard error:\n" << stdErr.str() << std::endl);
//     DEBUG_PRINT("Result: " << result << std::endl);
//     // Delete every file in the output directory
//     for(const auto& entry : std::filesystem::directory_iterator("output")) {
//         std::filesystem::remove(entry.path());
//     }
//     return result;
// }