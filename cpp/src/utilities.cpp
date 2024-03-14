/******************************************************************************
 * Filename:        utilities.cpp
 * Project:         eSRIM
 * Author:          Nathaniel Thomas
 * Date Created:    March 12, 2024
 * Date Modified:   March 13, 2024
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
 * See main.cpp for other revision history
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
    settings.displaySettings = false; // Default time flag

    // Check for help option
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            std::cout << helpMessage << std::endl;
            std::exit(EXIT_SUCCESS);
        }
    }

    // Check for settings display option
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-s" || arg == "--settings") {
            std::cout << settingsMessage << std::endl;
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
        } else if (arg == "-d" || arg == "--display") {
            settings.displaySettings = true;
        } else if (arg == "-p" || arg == "--progress") {
            settings.progress = true;
        } else {
            std::cerr << "Error: Unrecognized flag '" << arg << "'" << std::endl;
            std::cerr << "Use -h or --help for usage information" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    return settings;
}

bool promptContinue() {
    char answer;
    bool isValidInput = false;

    do {
        std::cout << "Continue? (y/n) ";
        std::cin >> answer;

        // Clear input buffer to handle incorrect input
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        if (answer == 'y' || answer == 'Y') {
            isValidInput = true;
            return true;
        } else if (answer == 'n' || answer == 'N') {
            isValidInput = true;
            return false;
        } else {
            std::cout << "Invalid input. Please enter 'y' or 'n'." << std::endl;
        }
    } while (!isValidInput);

    return false; // This line should never be reached, but added for completeness
}

void clearLine() {
    constexpr int width = 80;
    std::cout << "\r";
    for (size_t i = 0; i < width; ++i) {
        std::cout << " ";
    }
    std::cout << "\r";
    std::cout.flush();
}

#endif