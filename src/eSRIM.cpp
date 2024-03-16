/******************************************************************************
 * Filename:        eSRIM.cpp
 * Project:         eSRIM
 * Author:          Nathaniel Thomas
 * Date Created:    March 8, 2024
 * Date Modified:   March 15, 2024
 * File Version:    1.1
 * Group:           Dr. Shao's RMSLCF Group
 *
 * Description:
 * This is the main file for the eSRIM code developed by Dr. Shao and converted
 * to the C++ programming language
 *
 * Revision History
 * 1.0:
 * - File created and initial program structure made
 * 1.1: (March 13, 2024)
 * - Added File I/O
 * 1.2: (March 15, 2024)
 * - Renamed from main.cpp to eSRIM.cpp
 *****************************************************************************/
#ifndef eSRIM_CXX
#define eSRIM_CXX

// Local header files
#include "eSRIM.h"
#include "utilities.h"
// #include "eSRIM_classes.h" // No need, included in "main.h"

// Main program  start
int startESRIM(int argc, char* argv[], std::istream& inputStream) {

    // Pre-simulation code. Input and error checking, as well as signal handling. 
    Arguments arguments = parseCommandLine(argc, argv);
    std::shared_ptr<InputFields> input = InputFields::getInstance(arguments.filename);
    input->readSettingsFromFile(inputStream);

    if(arguments.progress) {
        input->setProgressChecking(true);
    }

    checkDisplayOption(arguments, input);

    checkHardwareThreads(input);

    // Primary simulation section
    auto start = std::chrono::high_resolution_clock::now();

    if(input->getType() == ELECTRON) {
        try{
            Electron::readParametersAndInitialize(input);
        } catch(std::ios_base::failure& e) {
            std::cerr << e.what() << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    std::shared_ptr<Simulation> simulation = std::make_shared<Simulation>(input);
    Particle::seedRandomGenerator();

    simulation->initiate();

    // Simulation complete. Cleanup code.

    simulation.reset();
    
    if(arguments.time) {
        clearLine();
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std::cout   << "Execution time: "
                    << duration.count()
                    << " seconds"
                    << std::endl;
    }

    return 0;
}

#endif