/******************************************************************************
 * Filename:        utilities.h
 * Project:         eSRIM
 * Author:          Nathaniel Thomas
 * Date Created:    March 12, 2024
 * Date Modified:   March 16, 2024
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
 *****************************************************************************/
#ifndef UTILITIES_H
#define UTILITIES_H

// Library includes
#include <cstdlib>
#include <cstddef>
#include <iostream>

// Local includes
#include "eSRIM_classes.h"

// Define a struct to store settings
struct Arguments {
    std::string filename;
    bool time;
    bool displaySettings;
    bool progress;
};

// Define a custom exception that carries an exit status
class ExitException : public std::exception {
public:
    ExitException(int status);
    int status;
};

class InputFields;


// This class is for handler USER IO, output file data is handled by the Simulation::writeData() function
class IOHandler {
    private:
        // Private constructor to prevent external instantiation
        IOHandler(
            std::shared_ptr<InputFields> input,
            std::istream& inputStream,
            std::ostream& outputStream,
            std::ostream& errorStream);

        // Static instance pointer for a single IOHandler instance
        static std::shared_ptr<IOHandler> instance;

        // Constants

        const std::string_view inputNotInitializedMessage = "InputFields input has not been initialized.";

        const std::string_view helpMessage = 
        R"(eSRIM Version 1.0
        Author:         Nathaniel Thomas
        Contact:        nathaniel@swbell.net
        Release date:   March 12, 2024

        This is a program for simulating the Stopping and Range of Ions in Matter (SRIM), with electron bombardment simulation capabilities.

        Usage: ./eSRIM [options][paths...]

        Options
            -f --filename <filename>    Read settings for eSRIM from <filename>. [Default="settings.txt"]
            -t --time                   Record execution time and output to standard output.
            -s --settings               Display an example settings file.   
            -d --display                Output the active settings to the standard output.
            -h --help                   Display this help message.
            -p --progress               Display the progress of the simulation while it is running to the standard output.
        )";

        std::string settingsMessage;

        // Private Members
        std::shared_ptr<InputFields> input;
        std::istream& inputStream;
        std::ostream& outputStream;
        std::ostream& errorStream;
        Arguments arguments;
        bool argumentsSet = false;

    public:

        // Instantiators
        static std::shared_ptr<IOHandler> getInstance();  
        
        static std::shared_ptr<IOHandler> getInstance(
            std::shared_ptr<InputFields> input,
            std::istream& inputStream = std::cin,
            std::ostream& outputStream = std::cout,
            std::ostream& errorStream = std::cerr);

        // Getters
        std::shared_ptr<InputFields> getInput() const;
        Arguments& getArguments();
        std::istream& getInputStream() const;
        std::ostream& getOutputStream() const;
        std::ostream& getErrorStream() const;

        // Setters
        void setInput(std::shared_ptr<InputFields> input);
        // No setters for reference members, as they cannot be rebound.

        // Functions
        void parseCommandLine(const int argc, const char* argv[]); // Add scope resolution operator
        bool promptContinue();
        void clearLine();
        void checkHardwareThreads();
        void checkDisplayOption(); // Add scope resolution operator
        void writeSettingsToFile(const std::string& settingsFilename);
        bool parseSetting(const std::string& line, std::string& name, std::string& value);
        bool stringToBool(const std::string& str);
        void handleFileOpenError(const std::string& filename);
        void readSettingsFromFile();
        void handleNoOutputDirectory(const std::filesystem::path& outputPath);
        void handleOutputFileOpenError(
            OutputType outputType,
            std::filesystem::path outputPath,
            std::string filename,
            int bombardmentID);
        bool areArgumentsSet() const;

        // Templated function
        template<typename T>
        IOHandler& operator<<(const T& message) {
            errorStream << message;
            return *this;
        }

        IOHandler& operator<<(std::ostream& (*func)(std::ostream&)) {
            func(errorStream);
            return *this;
        }
};

#endif