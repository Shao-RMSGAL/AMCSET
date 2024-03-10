/******************************************************************************
 * Filename:        main.cpp
 * Project:         eSRIM
 * Author:          Nathaniel Thomas
 * Date Created:    March 8, 2024
 * Date Modified:   March 8, 2024
 * File Version:    1.0
 * Group:           Dr. Shao's RMSLCF Group
 *
 * Description:
 * This is the main file for the eSRIM code developed by Dr. Shao and converted
 * to the C++ programming language
 *
 * Revision History
 * 1.0:
 * - File created and initial program structure made
 *
 *****************************************************************************/

// Header files
#include <cmath>


// Local header files
#include "eSRIM_classes.h"


// Main program  start
int main() {
    DEBUG_PRINT("Run starting...");

    InputFields input;

    DEBUG_PRINT("Initializing simulation using default values...");
    Simulation simulation(input);
    Particle::seedRandomGenerator();

    simulation.initiate();

    simulation.writeCoordinateData();

    return 0;
}