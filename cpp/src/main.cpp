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
int main(int argc, char* argv[]) {

    size_t numSimul;;

    if(argc != 2) {
        // std::cerr << "Usage: " << argv[0] << " <integer_argument>" << std::endl;
        // return 1;
        numSimul = Defaults::simulationCount;
    } else {
        numSimul = std::stoi(argv[1]);
    }

    // Comment this out if you want ion mode
    InputFields input(Constants::electronCharge, Defaults::electronEnergy,
            Constants::electronMass, numSimul,
            Defaults::substrateCharge, Defaults::substrateDensity,
            Defaults::substrateMass, ELECTRON, Defaults::range);

    // Comment this out if you want electron mode
    // InputFields input(Defaults::ionCharge, Defaults::ionEnergy,
    //         Defaults::ionMass, numSimul,
    //         Defaults::substrateCharge, Defaults::substrateDensity,
    //         Defaults::substrateMass, ION, Defaults::range);

    auto start = std::chrono::high_resolution_clock::now();

    #ifdef DEBUG_MODE
    const char* mode = (input.getType() == ION) ? "Ion" : "Electron";
    #endif

    DEBUG_PRINT("Run starting...");

    DEBUG_PRINT("Initializing simulation in " << mode << " mode..");

    if(input.getType() == ELECTRON) {
        Electron::readParametersAndInitialize();
    }

    std::shared_ptr<Simulation> simulation = std::make_shared<Simulation>(input);
    Particle::seedRandomGenerator();

    simulation->initiate();

    simulation.reset();

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration = end - start;
    std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;

    return 0;
}