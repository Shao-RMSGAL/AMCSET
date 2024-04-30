/******************************************************************************
 * Filename:        constants.cpp
 * Project:         AMCSET
 * Author:          Nathaniel Thomas
 * Date Created:    March 8, 2024
 * Date Modified:   April 30, 2024
 * File Version:    1.0
 * Group:           Dr. Shao's RMSLCF Group
 *
 * Description:
 * This file contains useful constants to be used in the AMCSET program.
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

// Options to add
// electronMode
// electronEnergy
// electronStoppingEnergy
// electronScreeningParametersFilename
// numElecScreeningPotentialElements
// mottScatteringParametersFilename
// numMottScatteringPotentialElements
// numAngleDivisors
// numFlyingDistances
// enableDamageCascade
// ionCharge
// ionEnergy
// ionMass
// substrateDisplacementEnergy
// ionStoppingEnergy
// logSingleDisplacement
// inputDirectoryName
// outputCoordinateFilename
// outputFileEndMarker
// outputFileExtension
// outputDirectory
// substrateCharge
// substrateMass
// logEndOfFlyingDistanceOnly
// logStoppingPointOnly
// simulationCount
// numThreads

enum InputOptionType{
    ELECTRON_MODE,
    ELECTRON_ENERGY,
    ELECTRON_STOPPING_ENERGY,
    ELECTRON_SCREENING_PARAMETERS_FILENAME,
    NUM_ELEC_SCREENING_POTENTIAL_ELEMENTS,
    MOTT_SCATTERING_PARAMETERS_FILENAME,
    NUM_MOTT_SCATTERING_POTENTIAL_ELEMENTS,
    NUM_ANGLE_DIVISORS,
    NUM_FLYING_DISTANCES,
    ENABLE_DAMAGE_CASCADE,
    ION_CHARGE,
    ION_ENERGY,
    ION_MASS,
    SUBSTRATE_DISPLACEMENT_ENERGY,
    ION_STOPPING_ENERGY,
    LOG_SINGLE_DISPLACEMENT,
    INPUT_DIRECTORY_NAME,
    OUTPUT_COORDINATE_FILENAME,
    OUTPUT_FILE_END_MARKER,
    OUTPUT_FILE_EXTENSION,
    OUTPUT_DIRECTORY,
    SUBSTRATE_CHARGE,
    SUBSTRATE_MASS,
    LOG_END_OF_FLYING_DISTANCE_ONLY,
    LOG_STOPPING_POINT_ONLY,
    SIMULATION_COUNT,
    NUM_THREADS
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
    constexpr bool enableDamageCascade = false;
    constexpr double ionCharge = 14; // e, silicon
    constexpr double ionDisplacementEnergy = 0.04; // keV
    constexpr double ionStoppingEnergy = 0.04; // keV
    constexpr double ionEnergy = 50; // keV, default ion energy
    constexpr double ionMass = 27.97692653442; // amu, silicon 29
    constexpr double range = 80000000; // angstrom
    constexpr bool logEndOfFlyingDistanceOnly = false;
    constexpr bool logSingleDisplacement = false; // Disable to shrink output file size
    constexpr bool logStoppingPointOnly = true;
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
    constexpr ParticleType type = ION;
    constexpr double substrateCharge = 26; // e, iron
    constexpr double substrateMass = 55.9349363; // amu, iron 56
    constexpr const char* settingsFilename = "settings.txt";
    constexpr size_t simulationCount = 1;
    constexpr bool progressChecking = false;
    constexpr size_t numThreads = 1;

    // Dependent values
    constexpr double windowRange = Defaults::range;
    constexpr double substrateDensity = 7.874
                                        /Defaults::substrateMass
                                        *Constants::avogadrosNumber
                                        *(   Constants::cmPerAngstrom
                                            *Constants::cmPerAngstrom
                                            *Constants::cmPerAngstrom
                                        ); // g/cm^3 -> atoms/angstrom^3
}

#endif