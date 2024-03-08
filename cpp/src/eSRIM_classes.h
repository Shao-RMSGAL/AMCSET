/******************************************************************************
 * Filename:        eSRIM_classes.h
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
 * - File created and class declarations added
 *
 *****************************************************************************/


// Includes
#include <vector>

// Local header files
#include "constants.h"

// Class and function declarations (to be moved to a header file in a later 
// version)
struct Coordinate {
    double x;
    double y;
    double z;
};


class Particle {
    protected:
        std::vector<Coordinate> coordinate_vector;
        double mass; // in angstrom
        double charge; // in e
    public:

        // Modifiers
        virtual void add_coordinate(Coordinate coordinate);

        // Accessors
        virtual const std::vector<Coordinate>& get_coordinates() const;
};


class Ion : public Particle {
    private:

    public:

};


class Electron : public Particle {
    private:
        static constexpr double mass = Constants::electronMass;

    public:
};


class Bombardment {
    private:
        std::vector<Particle> particles;
};


class Simulation {
    private:

    public:
        
};