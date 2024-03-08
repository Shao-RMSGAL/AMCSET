/******************************************************************************
 * Filename:        eSRIM_classes.cpp
 * Project:         eSRIM
 * Author:          Nathaniel Thomas
 * Date Created:    March 8, 2024
 * Date Modified:   March 8, 2024
 * File Version:    1.0
 * Group:           Dr. Shao's RMSLCF Group
 *
 * Description:
 * This file contains the class declarations for the eSRIM program. 
 * implementations are in eSRIM_classes.cpp
 * 
 *
 * Revision History
 * 1.0:
 * - File created and class implementations added
 *
 *****************************************************************************/


// Local includes
#include "eSRIM_classes.h"


void Particle::add_coordinate(Coordinate coordinate) {
    coordinate_vector.push_back(coordinate);
};

const std::vector<Coordinate>& Particle::get_coordinates() const {
    return coordinate_vector;
};