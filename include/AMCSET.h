/******************************************************************************
 * Filename:        AMCSET.h
 * Project:         AMCSET
 * Author:          Nathaniel Thomas
 * Date Created:    March 15, 2024
 * Date Modified:   April 30, 2024
 * File Version:    1.1
 * Group:           Dr. Shao's RMSLCF Group
 *
 * Description:
 * This is the main header file for the AMCSET code developed by Dr. Shao and converted
 * to the C++ programming language by Nathaniel Thomaas
 *
 * Revision History
 * 1.0:
 * - File created and initial program structure made
 *****************************************************************************/
#ifndef AMCSET_H
#define AMCSET_H

#include <iostream>

// Main program  start
int startAMCSET(
    const int argc,
    const char* argv[],
    std::istream& inputStream = std::cin,
    std::ostream& outputStream = std::cout,
    std::ostream& errorStream = std::cerr
        );

#endif