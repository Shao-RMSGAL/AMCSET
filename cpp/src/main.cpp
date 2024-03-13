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
#ifndef MAIN_CXX
#define MAIN_CXX

// Local header files
#include "eSRIM_classes.h"
#include "utilities.h"


// Main program  start
int main(int argc, char* argv[]) {
    Arguments arguments = parseCommandLine(argc, argv);
    std::shared_ptr<InputFields> input = InputFields::getInstance(arguments.filename);
    input->readSettingsFromFile();
    
    auto start = std::chrono::high_resolution_clock::now();

    if(input->getType() == ELECTRON) {
        Electron::readParametersAndInitialize(input);
    }

    std::shared_ptr<Simulation> simulation = std::make_shared<Simulation>(input);
    Particle::seedRandomGenerator();

    simulation->initiate();
    simulation.reset();
    
    if(arguments.time) {
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;
    }

    return 0;
}

#endif