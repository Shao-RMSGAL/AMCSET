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

ExitException::ExitException(int status) : status(status) {};

// Main program  start
int startESRIM(
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
        
        // Arguments must be parsed before passing to InputFields
        ioHandler->parseCommandLine(argc, argv);

        std::shared_ptr<InputFields> input = InputFields::getInstance(ioHandler);
        ioHandler->setInput(input); // Update the input pointer

        ioHandler->readSettingsFromFile();

        if(ioHandler->getArguments().progress) {
            input->setProgressChecking(true);
        }

        ioHandler->checkDisplayOption();
        ioHandler->checkHardwareThreads();

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

        std::shared_ptr<Simulation> simulation = std::make_shared<Simulation>(input);
        Particle::seedRandomGenerator();

        simulation->initiate();

        // Simulation complete. Cleanup code.

        simulation.reset();
        
        if(ioHandler->getArguments().time) {
            ioHandler->clearLine();
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = end - start;
            *ioHandler  << "Execution time: "
                        << std::to_string(duration.count())
                        << " seconds"
                        << std::endl;
        }

        return 0;
    } catch (const ExitException& e) {
        return e.status;
    }
}

#endif