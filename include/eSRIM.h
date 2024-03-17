/******************************************************************************
 * Filename:        eSRIM.h
 * Project:         eSRIM
 * Author:          Nathaniel Thomas
 * Date Created:    March 15, 2024
 * Date Modified:   March 15, 2024
 * File Version:    1.1
 * Group:           Dr. Shao's RMSLCF Group
 *
 * Description:
 * This is the main header file for the eSRIM code developed by Dr. Shao and converted
 * to the C++ programming language by Nathaniel Thomaas
 *
 * Revision History
 * 1.0:
 * - File created and initial program structure made
 *****************************************************************************/
#ifndef eSRIM_H
#define eSRIM_H

#include <iostream>

// Main program  start
int startESRIM(
    const int argc,
    const char* argv[],
    std::istream& inputStream = std::cin,
    std::ostream& outputStream = std::cout,
    std::ostream& errorStream = std::cerr
        );

#endif