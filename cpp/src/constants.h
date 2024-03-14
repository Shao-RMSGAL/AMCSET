/******************************************************************************
 * Filename:        constants.cpp
 * Project:         eSRIM
 * Author:          Nathaniel Thomas
 * Date Created:    March 8, 2024
 * Date Modified:   March 13, 2024
 * File Version:    1.0
 * Group:           Dr. Shao's RMSLCF Group
 *
 * Description:
 * This file contains useful constants to be used in the eSRIM program.
 * 
 *
 * Revision History
 * 1.0:
 * - File created and constants added
 * See main.cpp for other revision history
 * 
 *****************************************************************************/
#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

enum OutputType{
    COORDINATE,
    ENDOFFILE
};

enum ParticleType{
    ION,
    ELECTRON,
    SUBSTRATE
};

namespace Constants {
    constexpr double avogadrosNumber = 6.02214076e23; // atoms/mol
    constexpr double cmPerAngstrom = 1e-8; // cm/Angstrom
    constexpr double electronCharge = -1; // e
    constexpr double electronChargeSI = 1.602176634e-19; // C
    constexpr double electronMass = 9.1093837015E-31; // kg
    constexpr double keVToeV = 1000; // eV/keV
    constexpr size_t numElecScreeningPotentialParams = 6;
    constexpr size_t numMottjParam = 5; 
    constexpr size_t numMottkParam = 6;
    constexpr double plancksConstant = 6.62607015e-34; // J·s
    constexpr double speedOfLight = 299792458; // m/s

    // Dependent values
    constexpr double mec2 = electronMass
                            *speedOfLight*speedOfLight
                            /electronChargeSI
                            /keVToeV; // keV
}

namespace Defaults {
    constexpr double electronEnergy = 10000; // keV, default electron energy
    constexpr const char* electronScreeningParametersFilename = "electron_screeening_potentials.csv";
    constexpr double electronStoppingEnergy = 1; // keV
    constexpr bool enableDamageCascade = true;
    constexpr double ionCharge = 14; // e, silicon
    constexpr size_t ionDisplacementEnergy = 0.04; // keV
    constexpr size_t ionStoppingEnergy = 0.04; // keV
    constexpr double ionEnergy = 50; // keV, default ion energy
    constexpr double ionMass = 27.97692653442; // amu, silicon 29
    constexpr double range = 80000000; // angstrom
    constexpr bool logEndOfFlyingDistanceOnly = false;
    constexpr bool logSingleDisplacement = false; // Disable to shrink output file size
    constexpr bool logStoppingPointOnly = false;
    constexpr size_t numAngleDivisors = 1000;
    constexpr size_t numElecScreeningPotentialElements = 92;
    constexpr size_t numFlyingDistances = 1000;
    constexpr size_t numMottElements = 118; 
    constexpr const char* inputDirectory =  "input";
    constexpr const char* mottScatteringParametersFilename = "mott_scattering_parameters.csv";
    constexpr const char* outputCoordinateFilename =  "coordinateOutput";
    constexpr const char* outputFileEnd = "End of file";
    constexpr const char* outputFileExtension = ".csv";
    constexpr const char* outputDirectory =  "output";
    constexpr ParticleType type = ELECTRON;
    constexpr double substrateCharge = 26; // e, iron
    constexpr double substrateMass = 55.9349363; // amu, iron 56
    constexpr const char* settingsFilename = "settings.txt";
    constexpr size_t simulationCount = 1;
    constexpr double windowRange = Defaults::range;

    // Dependent values
    constexpr double substrateDensity = 7.874
                                        /Defaults::substrateMass
                                        *Constants::avogadrosNumber
                                        *(   Constants::cmPerAngstrom
                                            *Constants::cmPerAngstrom
                                            *Constants::cmPerAngstrom
                                        ); // g/cm^3 -> atoms/angstrom^3
}

#endif