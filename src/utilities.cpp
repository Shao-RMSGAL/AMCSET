/******************************************************************************
 * Filename:        utilities.cpp
 * Project:         AMCSET
 * Author:          Nathaniel Thomas
 * Date Created:    March 12, 2024
 * Date Modified:   April 30, 2024
 * File Version:    1.0
 * Group:           Dr. Shao's RMSLCF Group
 *
 * Description:
 * This file contains useful utilities and constants used by the 
 * main function in the AMCSET program
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

std::shared_ptr<IOHandler> IOHandler::instance = nullptr;

IOHandler::IOHandler(
        std::shared_ptr<InputFields> input,
        std::istream& inputStream,
        std::ostream& outputStream,
        std::ostream& errorStream)
        : input(input),
        inputStream(inputStream),
        outputStream(outputStream),
        errorStream(errorStream),
        arguments({}), 
        argumentsSet(false) {
    std::stringstream settingsMessageStream;
    settingsMessageStream << "electronMode=" << ((Defaults::type == ELECTRON) ? "true" : "false") << std::endl;
    settingsMessageStream << "electronEnergy(keV)=" << Defaults::electronEnergy << std::endl;
    settingsMessageStream << "electronStoppingEnergy(keV)=" << Defaults::electronStoppingEnergy << std::endl;
    settingsMessageStream << "electronScreeningParametersFilename=\"" << Defaults::electronScreeningParametersFilename << "\"" << std::endl;
    settingsMessageStream << "numElecScreeningPotentialElements=" << Defaults::numElecScreeningPotentialElements << std::endl;
    settingsMessageStream << "mottScatteringParametersFilename=\"" << Defaults::mottScatteringParametersFilename << "\"" << std::endl;
    settingsMessageStream << "numMottScatteringPotentialElements=" << Defaults::numMottElements << std::endl;
    settingsMessageStream << "numAngleDivisors=" << Defaults::numAngleDivisors << std::endl;
    settingsMessageStream << "numFlyingDistances=" << Defaults::numFlyingDistances << std::endl;
    settingsMessageStream << "enableDamageCascade=" << ((Defaults::enableDamageCascade) ? "true" : "false") << std::endl;
    settingsMessageStream << "ionCharge(e)=" << Defaults::ionCharge << std::endl;
    settingsMessageStream << "ionEnergy(keV)=" << Defaults::ionEnergy << std::endl;
    settingsMessageStream << "ionMass(amu)=" << Defaults::ionMass << std::endl;
    settingsMessageStream << "substrateDisplacementEnergy(keV)=" << Defaults::ionDisplacementEnergy << std::endl;
    settingsMessageStream << "ionStoppingEnergy(keV)=" << Defaults::ionStoppingEnergy << std::endl;
    settingsMessageStream << "logSingleDisplacement=" << ((Defaults::logSingleDisplacement) ? "true" : "false") << std::endl;
    settingsMessageStream << "inputDirectoryName=\"" << Defaults::inputDirectory << "\"" << std::endl;
    settingsMessageStream << "outputCoordinateFilename=\"" << Defaults::outputCoordinateFilename << "\"" << std::endl;
    settingsMessageStream << "outputFileEndMarker=\"" << Defaults::outputFileEnd << "\"" << std::endl;
    settingsMessageStream << "outputFileExtension=\"" << Defaults::outputFileExtension << "\"" << std::endl;
    settingsMessageStream << "outputDirectory=\"" << Defaults::outputDirectory << "\"" << std::endl;
    settingsMessageStream << "substrateCharge(e)=" << Defaults::substrateCharge << std::endl;
    settingsMessageStream << "substrateMass(amu)=" << Defaults::substrateMass << std::endl;
    settingsMessageStream << "logEndOfFlyingDistanceOnly=" << ((Defaults::logEndOfFlyingDistanceOnly) ? "true" : "false") << std::endl;
    settingsMessageStream << "logStoppingPointOnly=" << ((Defaults::logStoppingPointOnly) ? "true" : "false") << std::endl;
    settingsMessageStream << "simulationCount=" << Defaults::simulationCount << std::endl;
    settingsMessageStream << "numThreads=" << Defaults::numThreads;
    settingsMessage = settingsMessageStream.str();
};

std::shared_ptr<IOHandler> IOHandler::getInstance()
{
    if(!instance) {
        throw std::logic_error("IOHandler instance has not been initialized.");
    }
    return instance;
}

std::shared_ptr<IOHandler> IOHandler::getInstance(
    std::shared_ptr<InputFields> input,
    std::istream& inputStream,
    std::ostream& outputStream,
    std::ostream& errorStream)
{
    if(!instance) {
        instance = std::shared_ptr<IOHandler>(new IOHandler(
            input,
            inputStream,
            outputStream,
            errorStream));
    }
    return instance;
}

// Getters
std::shared_ptr<InputFields> IOHandler::getInput() const { return input; }
Arguments& IOHandler::getArguments() { return arguments; }
std::istream& IOHandler::getInputStream() const { return inputStream; }
std::ostream& IOHandler::getOutputStream() const { return outputStream; }
std::ostream& IOHandler::getErrorStream() const { return errorStream; }

// Setters
void IOHandler::setInput(std::shared_ptr<InputFields> input) { this->input = input; }

// Function to parse command line arguments
void IOHandler::parseCommandLine(const int argc, const char* argv[]) {

    DEBUG_PRINT("-----------------------Progress Checking 116-----------------------");

    Arguments& settings = arguments;
    
    settings.filename = "settings.txt"; // Default filename
    settings.time = false; // Default time flag
    settings.displaySettings = false; // Default time flag
    settings.progress = false; // Default time flag

    // No arguments passed
    if(argc == 1) {
        DEBUG_PRINT("-----------------------Progress Checking 127-----------------------");
        return;
    }


    // Check for help option
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            outputStream << helpMessage << std::endl;
            throw ExitException(EXIT_SUCCESS);
        }
    }

    DEBUG_PRINT("-----------------------Progress Checking 141-----------------------");

    // Check for settings display option
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-s" || arg == "--settings") {
            outputStream << "Settings\n\n"
                    << "The settings file is a text file, usually named \"settings.txt\", "
                    << "but you can use a custom settings filename if you pass it via the command line."
                    << "Here is an example file of the default settings:\n\n"
                    << settingsMessage
                    << std::endl;
            throw ExitException(EXIT_SUCCESS);
        }
    }

    DEBUG_PRINT("-----------------------Progress Checking 157-----------------------");
    
    // Parse other options
    for (int i = 1; i < argc; ++i) {
        DEBUG_PRINT("Val: " << argc);
        std::string arg = argv[i];
        if (arg == "-f" || arg == "--filename") {
            if (i + 1 < argc) {
                settings.filename = argv[i + 1];
                ++i; // Skip next argument
            } else {
                errorStream << "Error: Missing argument for filename option" << std::endl;
                throw ExitException(EXIT_FAILURE);
            }
        } else if (arg == "-t" || arg == "--time") {
            settings.time = true;
        } else if (arg == "-d" || arg == "--display") {
            settings.displaySettings = true;
        } else if (arg == "-p" || arg == "--progress") {
            settings.progress = true;
        } else {
            errorStream << "Error: Unrecognized flag '" << arg << "'" << std::endl;
            errorStream << "Use -h or --help for usage information" << std::endl;
            throw ExitException(EXIT_FAILURE);
        }
    }

    argumentsSet = true;
}

bool IOHandler::promptContinue() {
    char answer;

    while(true) {
        outputStream << "Continue? (y/n) ";
        inputStream >> answer;

        // Clear input buffer to handle incorrect input
        inputStream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        if (answer == 'y' || answer == 'Y') {
            return true;
        } else if (answer == 'n' || answer == 'N') {
            return false;
        } else {
            outputStream << "Invalid input. Please enter 'y' or 'n'." << std::endl;
        }
    }
}

void IOHandler::clearLine() {
    constexpr int width = 120;
    outputStream << "\r";
    for (size_t i = 0; i < width; ++i) {
        outputStream << " ";
    }
    outputStream << "\r";
    outputStream.flush();
}

void IOHandler::checkHardwareThreads()
{
    if(input == nullptr) {
        throw std::logic_error(std::string(inputNotInitializedMessage));
    }

    DEBUG_PRINT("-----------------------Progress Checking 223-----------------------");


    size_t cores = std::thread::hardware_concurrency();
    DEBUG_PRINT("-----------------------Progress Checking 227-----------------------");
    size_t threads = input->getNumThreads();
    
    if (cores < threads)
    {
        DEBUG_PRINT("cores: " << cores << " threads: " << threads);
        DEBUG_PRINT("-----------------------Progress Checking 232-----------------------");

        const std::string ANSI_COLOR_GREEN = "\033[1;33m";
        const std::string ANSI_COLOR_RESET = "\033[0m";
        errorStream << ANSI_COLOR_GREEN
                << "Warning: Number of threads ("
                << threads
                << ") exceeds number of hardware threads ("
                << cores
                << "). Continuing may cause unexpected behavior on this system."
                << " Consider settings the \"numThreads\" setting to a value less"
                << " than or equal to "
                << cores
                << "."
                << ANSI_COLOR_RESET
                << "\nUse -s for settings help. Use -h for usage help."
                << std::endl;
        if (!promptContinue())
        {
            DEBUG_PRINT("-----------------------Progress Checking 251-----------------------");

            throw ExitException(EXIT_SUCCESS);
        }
    }
    DEBUG_PRINT("-----------------------Progress Checking 256-----------------------");

}

void IOHandler::checkDisplayOption()
{
    if(input == nullptr) {
        throw std::logic_error(std::string(inputNotInitializedMessage));
    }
    if (arguments.displaySettings)
    {
        const std::string ANSI_COLOR_GREEN = "\033[1;32m";
        const std::string ANSI_COLOR_RESET = "\033[0m";
        outputStream
            << ANSI_COLOR_GREEN
            << "Stimulation settings:\n\n"
            << ANSI_COLOR_RESET
            << input->printInputFields()
            << std::endl;
    }
}

void IOHandler::writeSettingsToFile(const std::string& settingsFilename) {
    // Open the file "settings.txt" for writing
    std::ofstream outputFile(settingsFilename);

    if (outputFile.is_open()) {
        outputFile << settingsMessage;
        outputFile.close();
        outputStream << "Created new \"" << Defaults::settingsFilename << "\" file.\n";
    } else {
        errorStream << "Error opening " 
            << Defaults::settingsFilename
            << " for writing.\n"
            << " Try again?";
        if (promptContinue()) {
            writeSettingsToFile(Defaults::settingsFilename);
        } else {
            throw ExitException(EXIT_SUCCESS);
        }
    }
}

bool IOHandler::parseSetting(const std::string& line, std::string& name, std::string& value) {
        std::istringstream iss(line);
    if (std::getline(iss, name, '=') && std::getline(iss, value)) {
        return true;
    }
    return false;
} 

bool IOHandler::stringToBool(const std::string& str) {
    return str == "True" || str == "true" || str == "1";
}

void IOHandler::handleFileOpenError(const std::string& settingsFilename) {

    DEBUG_PRINT("-----------------------Progress Checking 302-----------------------");

    if(input == nullptr) {
                DEBUG_PRINT("-----------------------Progress Checking 305-----------------------");
        throw std::logic_error(std::string(inputNotInitializedMessage));
    }

    DEBUG_PRINT("-----------------------Progress Checking 309-----------------------");

    // errorStream = std::cerr;

    DEBUG_PRINT("File that could not be opened: " << settingsFilename);
    errorStream << "Error: Unable to open file \""
                <<  settingsFilename << "\". ";
            // Check if the default settings file exists

    DEBUG_PRINT("-----------------------Progress Checking 315-----------------------");

    DEBUG_PRINT("Tryng to open file:" << Defaults::settingsFilename);
    std::ifstream defaultSettingsFile(Defaults::settingsFilename);
    DEBUG_PRINT("File opened:" << Defaults::settingsFilename);


    if (defaultSettingsFile.is_open()) {
                DEBUG_PRINT("-----------------------Progress Checking 317-----------------------");
        outputStream << "The file \""
                    << Defaults::settingsFilename
                    << "\" has been found. Continue with this file? Use -s for settings help" 
                    << std::endl;
    } else {
                DEBUG_PRINT("-----------------------Progress Checking 323-----------------------");
        errorStream << "Creating default settings file \""
                        << Defaults::settingsFilename
                        << "\". Use -s for settings help."
                        << std::endl;
    }

    if(promptContinue()) {
                DEBUG_PRINT("-----------------------Progress Checking 331-----------------------");
        arguments.filename = Defaults::settingsFilename;
        if(!defaultSettingsFile.is_open()) { 
                    DEBUG_PRINT("-----------------------Progress Checking 334-----------------------");
            writeSettingsToFile(arguments.filename);
        } else {
                    DEBUG_PRINT("-----------------------Progress Checking 337-----------------------");
            defaultSettingsFile.close();
        }
        return;
    } else {
        throw ExitException(EXIT_SUCCESS);
    }
}

void IOHandler::readSettingsFromFile() {

    DEBUG_PRINT("-----------------------Progress Checking 339-----------------------");

    if(input == nullptr) {
                DEBUG_PRINT("-----------------------Progress Checking 342-----------------------");
        throw std::logic_error(std::string(inputNotInitializedMessage));
    }

    DEBUG_PRINT("-----------------------Progress Checking 346-----------------------");

    std::ifstream file(arguments.filename);

    while (!file.is_open()) {
        DEBUG_PRINT("-----------------------Progress Checking 351-----------------------");
        handleFileOpenError(arguments.filename);
        file.open(arguments.filename);
        return;
    }
    
    std::string line;
    while (std::getline(file, line)) {
        DEBUG_PRINT("-----------------------Progress Checking 558-----------------------");
        DEBUG_PRINT("Line: " << line);
        std::string name, value;
        if (parseSetting(line, name, value)) {
            DEBUG_PRINT("Name: " << name << " Value: " << value);
            if (name == "electronMode") {
                input->setType(stringToBool(value) ? ELECTRON : ION);
            } else if (name == "electronEnergy(keV)") {
                input->setEnergy(std::stod(value));
            } else if (name == "electronStoppingEnergy(keV)") {
                input->setElectronStoppingEnergy(std::stod(value));
            } else if (name == "electronScreeningParametersFilename") {
                value.erase(std::remove(value.begin(), value.end(), '\"'), value.end());
                input->setElectronScreeningParametersFilename(value);
            } else if (name == "numElecScreeningPotentialElements") {
                input->setNumElecScreeningPotentialElements(std::stoi(value));
            } else if (name == "mottScatteringParametersFilename") {
                value.erase(std::remove(value.begin(), value.end(), '\"'), value.end());
                input->setMottScatteringParametersFilename(value);
            } else if (name == "numMottScatteringPotentialElements") {
                input->setNumMottScatteringPotentialElements(std::stoi(value));
            } else if (name == "numAngleDivisors") {
                input->setNumAngleDivisors(std::stoi(value));
            } else if (name == "numFlyingDistances") {
                input->setNumFlyingDistances(std::stoi(value));
            } else if (name == "enableDamageCascade") {
                input->setEnableDamageCascade(stringToBool(value));
            } else if (name == "ionCharge(e)") {
                input->setCharge(input->getType() == ELECTRON ? Constants::electronCharge : std::stod(value));
            } else if (name == "ionEnergy(keV)") {
                input->setEnergy(input->getType() == ELECTRON ? Defaults::electronEnergy : std::stod(value));
            } else if (name == "ionMass(amu)") {
                input->setMass(input->getType() == ELECTRON ? Constants::electronMass : std::stod(value));
            } else if (name == "substrateDisplacementEnergy(keV)") {
                input->setIonDisplacementEnergy(std::stod(value));
            } else if (name == "ionStoppingEnergy(keV)") {
                input->setIonStoppingEnergy(std::stod(value));
            } else if (name == "logSingleDisplacement") {
                input->setLogSingleDisplacement(stringToBool(value));
            } else if (name == "inputDirectoryName") {
                value.erase(std::remove(value.begin(), value.end(), '\"'), value.end());
                input->setInputDirectoryName(value);
            } else if (name == "outputCoordinateFilename") {
                value.erase(std::remove(value.begin(), value.end(), '\"'), value.end());
                input->setOutputCoordinateFilename(value);
            } else if (name == "outputFileEndMarker") {
                value.erase(std::remove(value.begin(), value.end(), '\"'), value.end());
                input->setOutputFileEndMarker(value);
            } else if (name == "outputFileExtension") {
                value.erase(std::remove(value.begin(), value.end(), '\"'), value.end());
                input->setOutputFileExtension(value);
            } else if (name == "outputDirectory") {
                value.erase(std::remove(value.begin(), value.end(), '\"'), value.end());
                input->setOutputDirectory(value);
            } else if (name == "substrateCharge(e)") {
                input->setSubstrateCharge(std::stod(value));
            } else if (name == "substrateMass(amu)") {
                input->setSubstrateMass(std::stod(value));
            } else if (name == "logEndOfFlyingDistanceOnly") {
                input->setLogEndOfFlyingDistanceOnly(stringToBool(value));
            } else if (name == "logStoppingPointOnly") {
                input->setLogStoppingPointOnly(stringToBool(value));
            } else if (name == "simulationCount") {
                input->setSimulationCount(std::stoi(value));
            } else if (name == "numThreads") {
                input->setNumThreads(std::stoi(value));
            } else
                // Handle other `tings here if needed
                errorStream << "Warning: Unknown setting name '" << name << "'" << std::endl;
            }
        }
    file.close();
}

void IOHandler::handleNoOutputDirectory(const std::filesystem::path& outputPath) {
    const std::string ANSI_COLOR_RED = "\033[1;31m";
    const std::string ANSI_COLOR_RESET = "\033[0m";
    errorStream << "Unable to access output directory " << outputPath.string()
        << ".\nCreating directory "
        << outputPath.string()
        << ".\n"
        << ANSI_COLOR_RED << "Warning: Cancellation at this point will lead to loss of simulation data." << ANSI_COLOR_RESET
        << std::endl;
    if(promptContinue()) {
        fs::create_directories(outputPath);
    } else {
        exit(EXIT_FAILURE);
    }
}

void IOHandler::handleOutputFileOpenError(
        OutputType outputType,
        std::filesystem::path outputPath,
        std::string filename,
        int bombardmentID) {

    if(input == nullptr) {
        throw std::logic_error(std::string(inputNotInitializedMessage));
    }

    if(outputType != ENDOFFILE) {
        errorStream << "Failed to open output file "
                    << (outputPath / filename).string() 
                    << ".\nData will be lost for simulation "
                    << bombardmentID
                    << "."
                    << std::endl;
        return;
    } else {
        errorStream << "Failed to open output file "
                    << (outputPath / filename).string()
                    << ".\nFile end marker \""
                    << input->getOutputFileEndMarker() << "\""
                    << " could not be appended." 
                    << std::endl;
    }
}

bool IOHandler::areArgumentsSet() const { return argumentsSet; };

#endif