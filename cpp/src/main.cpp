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

    if(argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <integer_argument>" << std::endl;
        return 1;
    }

    size_t numSimul = std::stoi(argv[1]);

    DEBUG_PRINT("Run starting...");

    auto start = std::chrono::high_resolution_clock::now();

    InputFields input(Defaults::ionCharge, Defaults::ionEnergy,
                      Defaults::ionMass, numSimul,
                      Defaults::substrateCharge, Defaults::substrateDensity,
                      Defaults::substrateMass, Defaults::type, Defaults::range);

    DEBUG_PRINT("Initializing simulation using default values...");
    std::shared_ptr<Simulation> simulation = std::make_shared<Simulation>(input);
    Particle::seedRandomGenerator();

    simulation->initiate();

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration = end - start;
    std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;

    return 0;
}