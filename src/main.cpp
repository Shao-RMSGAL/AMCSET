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
int main(/*const int argc, const char* argv[]*/) {
    const char* argv[] = {"eSRIM", "-xxxx", nullptr};
    const int argc = sizeof(argv) / sizeof(argv[0]) - 1;
    return startESRIM(argc, argv);
}