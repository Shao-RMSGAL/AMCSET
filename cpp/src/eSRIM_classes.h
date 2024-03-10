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
#include <iostream> 
#include <memory>
#include <random>
#include <tuple>
#include <vector>

// Local header files
#include "constants.h"


// Debug compiler directives 
#ifdef DEBUG_MODE
#define DEBUG_PRINT(x) std::cout << x << std::endl;
#else
#define DEBUG_PRINT(x)
#endif


// Class and function declarations (to be moved to a header file in a later 
// version)
struct Coordinate {
    double x;
    double y;
    double z;
};

struct Velocity{
    double zAngle;
    double xAngle;
    double energy;
};


class InputFields {
    private:
        const double charge;
        const double energy;
        const double mass;
        const size_t simulationCount;
        const double substrateCharge;
        const double substrateDensity;
        const double substrateMass;
        const ParticleType type;
        const double range;

    public:
        // Constructors
        InputFields();
        InputFields(double charge, double energy,
                    double mass, size_t simulationCount,
                    double substrateCharge, double substrateDensity,
                    double substrateMass, ParticleType type, double range);

        // Accessors
        double getCharge() const;
        double getEnergy() const;
        double getMass() const;
        size_t getSimulationCount() const;
        double getSubstrateCharge() const;
        double getSubstrateDensity() const;
        double getSubstrateMass() const;
        ParticleType getType() const;
        double getRange() const;
};

class Bombardment;

class Particle {
    protected:
        std::vector<Coordinate> coordinate_vector;
        Coordinate coordinate;
        Velocity velocity;
        double &energy = velocity.energy;
        const double charge;
        const double mass;
        const double substrateCharge;
        const double substrateDensity;
        const double substrateMass;
        const ParticleType type;
        const double range;
        static std::mt19937 randomGenerator;
        static std::uniform_real_distribution<double> randomDistribution;
        std::weak_ptr<Bombardment> bombardment;

    public:
        // Constructor
        Particle(Coordinate coordinate, Velocity velocity, InputFields& input,
                std::weak_ptr<Bombardment> bombardment);

        Particle(Coordinate coordinate, Velocity velocity, double charge,
                double mass, double substrateCharge, double substrateDensity,
                double substrateMass, ParticleType type, double range,
                std::weak_ptr<Bombardment> bombardment); 

        // Destructor
        virtual ~Particle();

        // Modifiers
        virtual void addCoordinate(Coordinate coordinate);

        // Accessors
        virtual const std::vector<Coordinate>& getCoordinates() const;
        static double random();

        // Member functions
        virtual void fire();
        double atomicSpacing() const;
        std::tuple<Velocity, Velocity> relativeToAbsoluteVelocity(double angle,
            double targetAngle, double targetEnergy);
        double sign(double x);

        // Static functions
        static void seedRandomGenerator();
};  


class Ion : public Particle {
    public:
        // Constructor
        Ion(Coordinate coordinate, Velocity velocity, InputFields& input,
            std::weak_ptr<Bombardment> bombardment);

        Ion(Coordinate coordinate, Velocity velocity, double charge,
                double mass, double substrateCharge, double substrateDensity,
                double substrateMass, ParticleType type, double range,
                std::weak_ptr<Bombardment> bombardment); 

        // Actions
        void fire() override;
        double F(double X, double COLUMBIAVK, double AU);
        double DF(double X, double COLUMBIAVK, double AU);
        std::tuple<Velocity,Velocity> recoilEnergyAndVelocity();
        double electronicStoppingEnergy();
};


class Substrate : public Ion {

    public:
        // Constructors
        Substrate(Coordinate coordinate, Velocity velocity, double charge,
            double mass, double density, double range,
            std::weak_ptr<Bombardment> bombardment);
};


class Electron : public Particle {
    private:
        static constexpr double mass = Constants::electronMass;

    public:
};

class Bombardment : public std::enable_shared_from_this<Bombardment> {
    private:
        std::vector<std::unique_ptr<Particle>> particles;
    
    public:
        // Constructors
        Bombardment();

        // Destructor
        ~Bombardment();

        // Accessor
        const std::vector<std::unique_ptr<Particle>>& getParticles() const;

        void initiate(InputFields& inputFields);
        void addParticle(std::unique_ptr<Particle> particle);
};

class Simulation {
    private:
        InputFields inputs;
        std::vector<std::shared_ptr<Bombardment>> bombardments;
    public:
        // Constructors
        Simulation();
        Simulation(InputFields input);

        // Destructors
        ~Simulation();

        // Accessors
        InputFields getInputs();
        size_t getBombardmentSize();
        size_t getBombardmentCapacity();

        // Modifiers
        void addParticle();

        // Functions
        void initiate();
        void writeCoordinateData();
};