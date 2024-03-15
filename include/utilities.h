/******************************************************************************
 * Filename:        utilities.h
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
 * See main.cpp for other revision history
 *****************************************************************************/
#ifndef UTILITIES_H
#define UTILITIES_H

// Library includes
#include <cstdlib>
#include <cstddef>
#include <iostream>

// Local includes
#include "constants.h"
#include "eSRIM_classes.h"

// Define a struct to store settings
struct Arguments {
    std::string filename;
    bool time;
    bool displaySettings;
    bool progress;
};

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

const std::string_view settingsMessage = 
R"(electronMode=true
electronEnergy(keV)=10000
electronStoppingEnergy(keV)=1;
electronScreeningParametersFilename="electron_screeening_potentials.csv"
numElecScreeningPotentialElements=92
mottScatteringParametersFilename="mott_scattering_parameters.csv"
numMottScatteringPotentialElements=118
numAngleDivisors=1000
numFlyingDistances=1000
enableDamageCascade=false
ionCharge(e)=14
ionEnergy(keV)=50
ionMass(amu)=27.97692653442
substrateDisplacementEnergy(keV)=0.04
ionStoppingEnergy(keV)=0.04   
logSingleDisplacement=false
inputDirectoryName="input"
outputCoordinateFilename="coordinateOutput"
outputFileEndMarker="End of file"
outputFileExtension=".csv"
outputDirectory="output"
substrateCharge(e)=26
substrateMass(amu)=55.9349363
settingsFilename="settings.txt"
logEndOfFlyingDistanceOnly=false
logStoppingPointOnly=true
simulationCount=1
numThreads=8
)";

// Function to parse command line arguments
Arguments parseCommandLine(int argc, char* argv[]);

// Function to prompt user to continue when an error is encountered
bool promptContinue();

// Function for clearing the standard output
void clearLine();

void checkHardwareThreads(std::shared_ptr<InputFields> &input);

void checkDisplayOption(Arguments &arguments, std::shared_ptr<InputFields> &input);

void writeSettingsToFile();

#endif