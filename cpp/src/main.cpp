/******************************************************************************
 * Filename:        main.cpp
 * Project:         eSRIM
 * Author:          Nathaniel Thomas
 * Date Created:    March 8, 2024
 * Date Modified:   March 8, 2024
 * File Version:    1.0
 * Group:           Dr. Shao's RMSLCF Group
 *
 * Description:
 * This is the main file for the eSRIM code developed by Dr. Shao and converted
 * to the C++ programming language
 *
 * Revision History
 * 1.0:
 * - File created and initial program structure made
 *
 *****************************************************************************/

// Header files
#include <iostream> 
#include <cmath>


// Local header files
#include "eSRIM_classes.h"


// Debug compiler directives 
#ifdef DEBUG_MODE
#define DEBUG_PRINT(x) std::cout << x << std::endl;
#else
#define DEBUG_PRINT(x)
#endif

// Main program  start
int main() {
    DEBUG_PRINT("Test");

    Bombardment sim;

    return 0;
}