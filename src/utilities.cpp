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
#include <fstream>

// Local includes
#include "constants.h"
#include "utilities.h"


// Function to parse command line arguments
Arguments parseCommandLine(int argc, char* argv[]) {
    Arguments settings;
    settings.filename = "settings.txt"; // Default filename
    settings.time = false; // Default time flag
    settings.displaySettings = false; // Default time flag
    settings.progress = false; // Default time flag

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
            std::cout << "Settings\n\n"
                      << "The settings file is a text file, usually named \"settings.txt\", "
                      << "but you can use a custom settings filename if you pass it via the command line."
                      << "Here is an example file of the default settings:\n\n"
                      << settingsMessage
                      << std::endl;
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
    constexpr int width = 120;
    std::cout << "\r";
    for (size_t i = 0; i < width; ++i) {
        std::cout << " ";
    }
    std::cout << "\r";
    std::cout.flush();
}

void checkHardwareThreads(std::shared_ptr<InputFields> &input)
{
    unsigned int cores = std::thread::hardware_concurrency();

    if (cores < input->getNumThreads())
    {
        const std::string ANSI_COLOR_GREEN = "\033[1;33m";
        const std::string ANSI_COLOR_RESET = "\033[0m";
        std::cerr << ANSI_COLOR_GREEN
                  << "Warning: Number of threads ("
                  << input->getNumThreads()
                  << ") exceeds number of hardware threads ("
                  << std::thread::hardware_concurrency()
                  << "). Continuing may cause unexpected behavior on this system."
                  << " Consider settings the \"numThreads\" setting to a value less"
                  << " than or equal to "
                  << std::thread::hardware_concurrency()
                  << "."
                  << ANSI_COLOR_RESET
                  << "\nUse -s for settings help. Use -h for usage help."
                  << std::endl;
        if (!promptContinue())
        {
            std::exit(EXIT_SUCCESS);
        }
    }
}

void checkDisplayOption(Arguments &arguments, std::shared_ptr<InputFields> &input)
{
    if (arguments.displaySettings)
    {
        const std::string ANSI_COLOR_GREEN = "\033[1;32m";
        const std::string ANSI_COLOR_RESET = "\033[0m";
        std::cout
            << ANSI_COLOR_GREEN
            << "Stimulation settings:\n\n"
            << ANSI_COLOR_RESET
            << input->printInputFields()
            << std::endl;
    }
}

void writeSettingsToFile() {
    // Open the file "settings.txt" for writing
    std::ofstream outputFile(Defaults::settingsFilename);

    if (outputFile.is_open()) {
        outputFile << settingsMessage;
        outputFile.close();
        
        std::cout << "Created new \"" << Defaults::settingsFilename << "\" file.\n";
    } else {
        std::cerr << "Error opening " << Defaults::settingsFilename << " for writing.\n";
    }
}

#endif