/******************************************************************************
 * Filename:        main.h
 * Project:         eSRIM
 * Author:          Nathaniel Thomas
 * Date Created:    March 15, 2024
 * Date Modified:   March 15, 2024
 * File Version:    1.1
 * Group:           Dr. Shao's RMSLCF Group
 *
 * Description:
 * This is the main header file for the eSRIM code developed by Dr. Shao and converted
 * to the C++ programming language by Nathaniel Thomaas
 *
 * Revision History
 * 1.0:
 * - File created and initial program structure made
 *****************************************************************************/
#include "utilities.h"

#ifndef MAIN_H
#define MAIN_H

// Main program  start
int startESRIM(int argc, char* argv[]);

#endif

void checkHardwareThreads(std::shared_ptr<InputFields> &input);

void checkDisplayOption(Arguments &arguments, std::shared_ptr<InputFields> &input);
