/******************************************************************************
 * Filename:        utilities.cpp
 * Project:         eSRIM
 * Author:          Nathaniel Thomas
 * Date Created:    March 12, 2024
 * Date Modified:   March 12, 2024
 * File Version:    1.0
 * Group:           Dr. Shao's RMSLCF Group
 *
 * Description:
 * This file contains useful utilities and constants used by the 
 * main function in the eSRIM program
 * 
 *
 * Revision History
 * 1.0:
 * - File created
 *
 *****************************************************************************/
// Compiler directives
#ifndef UTILITIES_CXX
#define UTILITIES_CXX

// Library includes
#include <cstdlib>
#include <cstddef>
#include <iostream>

// Local includes
#include "constants.h"
#include "utilities.h"


// Function to parse command line arguments
Arguments parseCommandLine(int argc, char* argv[]) {
    Arguments settings;
    settings.filename = "settings.txt"; // Default filename
    settings.time = false; // Default time flag

    // Check for help option
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            std::cout << "eSRIM Version 1.0\n"<<
                        "Author:\t\tNathaniel Thomas\n" << 
                        "Contact:\tnathaniel@swbell.net\n" << 
                        "Release date:\tMarch 12, 2024\n\n" <<
                        "This is a program for simulating the Stopping and Range of Ions in Matter (SRIM), with electron bombardment simulation capabilities.\n" <<
                        "Usage: ./eSRIM [options][paths...]\n\n" <<
                        "Options\n" <<
                        " -f --filename <filename>\tRead settings for eSRIM from <filename>. [Default=\"settings.txt\"]\n" <<
                        " -t --time\t\t\tRecord execution time and output to standard output.\n"<<
                        " -h --help\t\t\tDisplay this help message.\n" <<
                        std::flush;
            std::exit(EXIT_SUCCESS);
        }
    }

    // Parse other options
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-f" || arg == "--filename") {
            if (i + 1 < argc) {
                settings.filename = argv[i + 1];
                ++i; // Skip next argument
            } else {
                std::cerr << "Error: Missing argument for filename option" << std::endl;
                std::exit(EXIT_FAILURE);
            }
        } else if (arg == "-t" || arg == "--time") {
            settings.time = true;
        } else {
            std::cerr << "Error: Unrecognized flag '" << arg << "'" << std::endl;
            std::cerr << "Use -h or --help for usage information" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    return settings;
}

#endif