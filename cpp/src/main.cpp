/******************************************************************************
 * Filename:        main.cpp
 * Project:         eSRIM
 * Author:          Nathaniel Thomas
 * Date Created:    March 8, 2024
 * Date Modified:   March 13, 2024
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

    if(arguments.progress) {
        input->setProgressChecking(true);
    } 

    if(arguments.displaySettings) {
        const std::string ANSI_COLOR_GREEN = "\033[1;32m";
        const std::string ANSI_COLOR_RESET = "\033[0m";
        std::cout 
        << ANSI_COLOR_GREEN
        << "Stimulation settings:\n\n"
        << ANSI_COLOR_RESET
        << input->printInputFields() 
        << std::endl;
    }

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