/******************************************************************************
 * Filename:        AMCSET.cpp
 * Project:         AMCSET
 * Author:          Nathaniel Thomas
 * Date Created:    March 8, 2024
 * Date Modified:   April 30, 2024
 * File Version:    1.1
 * Group:           Dr. Shao's RMSLCF Group
 *
 * Description:
 * This is the main file for the AMCSET code developed by Dr. Shao and converted
 * to the C++ programming language
 *
 * Revision History
 * 1.0:
 * - File created and initial program structure made
 * 1.1: (March 13, 2024)
 * - Added File I/O
 * 1.2: (March 15, 2024)
 * - Renamed from main.cpp to AMCSET.cpp
 *****************************************************************************/
#ifndef AMCSET_CXX
#define AMCSET_CXX

// Local header files
#include "AMCSET.h"
#include "utilities.h"

ExitException::ExitException(int status) : status(status) {};

// Main program  start
int startAMCSET(
        const int argc, const char* argv[],
        std::istream& inputStream,
        std::ostream& outputStream,
        std::ostream& errorStream) {

    try {

        // Pre-simulation code. Input and error checking, as well as signal handling. 

        // IOHandler and InputFields contain pointers to each other. Because of this,
        // one must be given a temporary "nullptr" as the pointer to the other
        // before being properly initialized.
        std::shared_ptr<IOHandler> ioHandler = IOHandler::getInstance(
            std::shared_ptr<InputFields>(nullptr) ,
            inputStream,
            outputStream,
            errorStream);
        
        DEBUG_PRINT("-----------------------Progress Checking 51-----------------------");


        // Arguments must be parsed before passing to InputFields
        ioHandler->parseCommandLine(argc, argv);

        DEBUG_PRINT("-----------------------Progress Checking 54-----------------------");

        std::shared_ptr<InputFields> input = InputFields::getInstance(ioHandler);
        ioHandler->setInput(input); // Update the input pointer

        DEBUG_PRINT("-----------------------Progress Checking 59-----------------------");

        ioHandler->readSettingsFromFile();

        DEBUG_PRINT("-----------------------Progress Checking 63-----------------------");

        if(ioHandler->getArguments().progress) {
            input->setProgressChecking(true);
        }

        DEBUG_PRINT("-----------------------Progress Checking 69-----------------------");

        ioHandler->checkDisplayOption();
        
        DEBUG_PRINT("-----------------------Progress Checking 74-----------------------");

        
        ioHandler->checkHardwareThreads();

        DEBUG_PRINT("-----------------------Progress Checking 75-----------------------");

        // Primary simulation section
        auto start = std::chrono::high_resolution_clock::now();

        if(input->getType() == ELECTRON) {
            try{
                Electron::readParametersAndInitialize(input);
            } catch(std::ios_base::failure& e) {
                errorStream << e.what() << std::endl;
                throw ExitException(EXIT_FAILURE);
            }
        }

        DEBUG_PRINT("-----------------------Progress Checking 88-----------------------");

        std::shared_ptr<Simulation> simulation = std::make_shared<Simulation>(input);
        Particle::seedRandomGenerator();

        DEBUG_PRINT("-----------------------Progress Checking 93-----------------------");

        simulation->initiate();

        DEBUG_PRINT("-----------------------Progress Checking 97-----------------------");

        // Simulation complete. Cleanup code.

        simulation.reset();
        
        DEBUG_PRINT("-----------------------Progress Checking 103-----------------------");

        if(ioHandler->getArguments().time) {
            ioHandler->clearLine();
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = end - start;
            *ioHandler  << "Execution time: "
                        << std::to_string(duration.count())
                        << " seconds"
                        << std::endl;
        }

        DEBUG_PRINT("-----------------------Progress Checking 115-----------------------");

        return 0;
    } catch (const ExitException& e) {
        DEBUG_PRINT("-----------------------Progress Checking 119-----------------------");
        return e.status;
    }
}

#endif