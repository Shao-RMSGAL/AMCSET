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
 *
 *****************************************************************************/
#ifndef UTILITIES_H
#define UTILITIES_H

// Library includes
#include <cstdlib>
#include <cstddef>
#include <iostream>

// Local includes
#include "constants.h"

// Define a struct to store settings
struct Arguments {
    std::string filename;
    bool time;
};

const std::string_view helpMessage = 
R"(eSRIM Version 1.0
Author:         Nathaniel Thomas
Contact:        nathaniel@swbell.net\n
Release date:   March 12, 2024\n\n 

This is a program for simulating the Stopping and Range of Ions in Matter (SRIM), with electron bombardment simulation capabilities.

Usage: ./eSRIM [options][paths...]

Options
    -f --filename <filename>\tRead settings for eSRIM from <filename>. [Default=\"settings.txt\"]
    -t --time\t\t\tRecord execution time and output to standard output.
    -h --help\t\t\tDisplay this help message.
    -s --settings\t\tDisplay an example settings file.
)";

const std::string_view settingsMessage = 
R"(Settings
    The settings file is a text file, usually named settings.txt, but you can use a custom settings filename if you pass it via the command line. Here is an example file of the default settings:

electronMode=true
electronEnergy(keV)=1000
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
simulationCount=1
)";

// Function to parse command line arguments
Arguments parseCommandLine(int argc, char* argv[]);

#endif