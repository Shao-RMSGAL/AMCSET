/******************************************************************************
 * Filename:        constants.cpp
 * Project:         eSRIM
 * Author:          Nathaniel Thomas
 * Date Created:    March 8, 2024
 * Date Modified:   March 8, 2024
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
 *
 *****************************************************************************/
#include <cmath>

#pragma once


enum ParticleType{
    ION,
    ELECTRON
};

namespace Constants {
    constexpr double avogadrosNumber = 6.02214076e23; // atoms/mol
    constexpr double cmPerAngstrom = 1e-8; // cm/Angstrom
    constexpr double ionDisplacementEnergy = 0.04; // keV
    constexpr double electronMass = 9.1093837015; // kg
    constexpr double electronCharge = -1; // e
    constexpr double plancksConstant = 6.62607015e-34; // J·s
    constexpr double speedOfLight = 299792458; // m/s
    constexpr double ionStoppingEnergy = 0.04; // keV
}

namespace Defaults {
    constexpr bool enableDamageCascade = false;
    constexpr double electronEnergy = 10000; // keV, default electron energy
    constexpr double ionEnergy = 50; // keV, default ion energy
    constexpr double ionCharge = 14; // e, silicon
    constexpr double ionMass = 27.97692653442; // amu, silicon 29
    constexpr double range = 80000000; // angstrom
    constexpr double substrateCharge = 26; // e, iron
    constexpr double substrateMass = 55.9349363; // amu, iron 56
    constexpr double substrateDensity = 7.874
                                        /Defaults::substrateMass
                                        *Constants::avogadrosNumber
                                        *(   Constants::cmPerAngstrom
                                            *Constants::cmPerAngstrom
                                            *Constants::cmPerAngstrom
                                        ); // g/cm^3 -> atoms/angstrom^3
    constexpr ParticleType type = ION;
    constexpr double windowRange = Defaults::range;
    constexpr size_t simulationCount = 1; // Temporarily set to 1 for ease of development
}