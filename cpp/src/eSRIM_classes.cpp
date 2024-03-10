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


// Header files
#include <chrono>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <typeinfo>


// Local includes
#include "eSRIM_classes.h"

// Namespaces
namespace fs = std::filesystem;

// Particle functions and static member initializers
std::mt19937 Particle::randomGenerator;
std::uniform_real_distribution<double> Particle::randomDistribution(0,1);

void Particle::seedRandomGenerator() {
    randomGenerator.seed(std::chrono::system_clock::now().time_since_epoch().count());
};

Particle::Particle(Coordinate coordinate, Velocity velocity, InputFields& input,
    std::weak_ptr<Bombardment> bombardment)
    : coordinate(coordinate), velocity(velocity), charge(input.getCharge()),
    mass(input.getMass()), substrateCharge(input.getSubstrateCharge()),
    substrateDensity(input.getSubstrateDensity()),
    substrateMass(input.getSubstrateMass()),
    type(input.getType()), range(input.getRange()), bombardment(bombardment) {};

Particle::Particle(Coordinate coordinate, Velocity velocity, double charge,
    double mass, double substrateCharge, double substrateDensity,
    double substrateMass, ParticleType type, double range,
    std::weak_ptr<Bombardment> bombardment)
    : coordinate(coordinate), velocity(velocity), charge(charge), mass(mass),
    substrateCharge(substrateCharge), substrateDensity(substrateDensity),
    substrateMass(substrateMass), type(type), range(range),
    bombardment(bombardment) {};

Particle::~Particle() {};

void Particle::addCoordinate(Coordinate coordinate) {
    // DEBUG_PRINT("Pushing back coordinate...");
    coordinate_vector.push_back(coordinate);
};

const std::vector<Coordinate>& Particle::getCoordinates() const {
    return coordinate_vector;
};

double Particle::random() {
    return randomDistribution(randomGenerator);
};

void Particle::fire() {
    throw std::logic_error("Particle base class fire() function called.");
    return;
};

double Particle::atomicSpacing() const {
    return std::pow(substrateDensity, -1.0/3.0);
};

std::tuple<Velocity, Velocity> Particle::relativeToAbsoluteVelocity(
        double THETA1RELATIVE, double THETA2RELATIVE, double targetEnergy) {
    double ALPHA1RELATIVE = random()*2.0*M_PI;
    double ALPHA2RELATIVE = ALPHA1RELATIVE + M_PI;

    double& THETAO = velocity.zAngle;
    double& ALPHAO = velocity.xAngle;

    // Angle
    double X1 = std::sin(THETA1RELATIVE) * std::cos(ALPHA1RELATIVE);
    double Y1 = std::sin(THETA1RELATIVE) * std::sin(ALPHA1RELATIVE);
    double Z1 = std::cos(THETA1RELATIVE);

    double Y0 = Y1 * std::cos(THETAO) + Z1*std::sin(THETAO);
    double Z0 = -Y1*std::sin(THETAO)+Z1*std::cos(THETAO);
    double X0 = X1;

    double Z = Z0;
    double X = X0*std::sin(ALPHAO)+Y0*std::cos(ALPHAO);
    double Y = -X0*std::cos(ALPHAO)+Y0*std::sin(ALPHAO);

    double THETA1;
    if (Z < 0.0) {
        THETA1 = M_PI+std::atan(std::sqrt(X*X + Y*Y)/Z);
    } else if (Z == 0.0) {
        THETA1 = M_PI/2.0;
    } else {
        THETA1 = std::atan(std::sqrt(X*X+Y*Y)/Z);
    } 

    double ALPHA1;
    if(std::sin(THETA1)!= 0.0) {
        if(X < 0.0) {
            ALPHA1 = M_PI+std::atan(Y/X); 
        } else if(X == 0.0) {
            ALPHA1 = M_PI - sign(Y)*M_PI/2.0; 
        } else {
            ALPHA1 = std::atan(Y/X); 
        }
    } else {
        ALPHA1 = 0.0;
    }

    // Target

    X1 = std::sin(THETA2RELATIVE)*std::cos(ALPHA2RELATIVE);
    Y1 = std::sin(THETA2RELATIVE)*std::sin(ALPHA2RELATIVE);
    Z1 = std::cos(THETA2RELATIVE);

    Y0 = Y1*std::cos(THETAO)+Z1*std::sin(THETAO);
    Z0 = -Y1*std::sin(THETAO)+Z1*std::cos(THETAO);
    X0 = X1;

    Z = Z0;
    X = X0 * std::sin(ALPHAO) + Y0 * std::cos(ALPHAO);
    Y = -X0 * std::cos(ALPHAO) + Y0 * std::sin(ALPHAO);

    double THETA2;
    if(Z < 0.0) {
        THETA2 = M_PI+std::atan(std::sqrt(X*X+Y*Y)/Z);
    } else if(Z == 0.0) {
        THETA2 = M_PI/2;
    } else {
        THETA2 = std::atan(std::sqrt(X*X+Y*Y)/Z);
    }

    double ALPHA2;
    if(std::sin(THETA2) != 0.0) {
        if(X < 0.0) {
            ALPHA2 = M_PI+std::atan(Y/X);
        }
        if(X == 0.0) {
             ALPHA2 = M_PI - 0.5*sign(Y)*M_PI;
            }
        if(X > 0.0) {
            ALPHA2 = std::atan(Y/X);
        }
    } else {
        ALPHA2 = 0.0;
    }
    Velocity newVelocity = {THETA1, ALPHA1, energy};
    Velocity targetVelocity = {THETA2, ALPHA2, targetEnergy};

    return std::make_tuple(newVelocity, targetVelocity);
};

double Particle::sign(double x) {
    if(x < 0.0) {
        return -1.0;
    } else if(x == 0.0) {
        return 0.0;
    } else {
        return 1;
    }
};

//  Ion functions
Ion::Ion(Coordinate coordinate, Velocity velocity, InputFields& input,
    std::weak_ptr<Bombardment> bombardment)
    : Particle(coordinate, velocity, input, bombardment) {
    DEBUG_PRINT("Ion made");
};

Ion::Ion(Coordinate coordinate, Velocity velocity, double charge,
    double mass, double substrateCharge, double substrateDensity,
    double substrateMass, ParticleType type, double range,
    std::weak_ptr<Bombardment> bombardment)
    : Particle(coordinate, velocity, charge, mass, substrateCharge, 
    substrateDensity, substrateMass, type, range, bombardment) {
    DEBUG_PRINT("Ion made");
}; 


void Ion::fire() {
    DEBUG_PRINT("Mass: " << mass);
    size_t fireCount = 0;
    while(velocity.energy > Constants::ionStoppingEnergy) {
        DEBUG_PRINT("Energy " << fireCount << ":\t" << velocity.energy);
        velocity.energy -= electronicStoppingEnergy();
        
        if(velocity.energy < Constants::ionStoppingEnergy) {
            return;
        }

        Velocity newVelocity;
        Velocity targetVelocity;

        std::tie(newVelocity, targetVelocity) = recoilEnergyAndVelocity();

        // DEBUG_PRINT("New velocity zAngle: " << newVelocity.zAngle);
        // DEBUG_PRINT("New velocity xAngle: " << newVelocity.xAngle);
        // DEBUG_PRINT("New velocity energy: " << newVelocity.energy);
        // DEBUG_PRINT("Target velocity zAngle: " << targetVelocity.zAngle);
        // DEBUG_PRINT("Target velocity xAngle: " << targetVelocity.xAngle);
        // DEBUG_PRINT("Target velocity energy: " << targetVelocity.energy);
        // DEBUG_PRINT("Random number: " << random());

        Coordinate newCoordinate = {
                coordinate.x + atomicSpacing()*std::sin(velocity.zAngle)
                                *std::cos(velocity.xAngle),
                coordinate.y + atomicSpacing()*std::sin(velocity.zAngle)
                                *std::sin(velocity.xAngle),
                coordinate.z + atomicSpacing()*std::cos(velocity.zAngle)
            };

        // TODO: Uncomment this code to enable damage cascades
        // if(energy - targetVelocity.energy > Constants::ionDisplacementEnergy) {
        //     std::unique_ptr<Substrate> targetParticle(new Substrate(newCoordinate,
        //         targetVelocity, substrateCharge, substrateMass, 
        //         substrateDensity, range, bombardment));

        //     DEBUG_PRINT("Attemtping lock...");
        //     if(auto locked = bombardment.lock()) {    
        //         locked->addParticle(std::move(targetParticle));
        //     } else {
        //         std::cerr << "Failed to create knock-on particle" << std::endl;
        //     }
        //     DEBUG_PRINT("Knock-off substrate ion created");
        // }

        DEBUG_PRINT("x, y, z: " << newCoordinate.x << ", "<< newCoordinate.y << ", "<< newCoordinate.z );
        addCoordinate(coordinate);
        coordinate = newCoordinate;
        velocity = newVelocity;
        fireCount++;
    }
};

double Ion::F(double X, double COLUMBIAVK, double AU) {
    if(X == 0.0) {
        return 0.0;
    } else {
        return COLUMBIAVK*X*(0.35*std::exp(-0.3/X/AU)+0.55*std::exp(-1.2/X/AU)
                +0.1*std::exp(-6.0/X/AU));
    }
};

double Ion::DF(double X, double COLUMBIAVK, double AU) {
    return F(X, COLUMBIAVK, AU)/X+COLUMBIAVK/X*(0.35*std::exp(-0.3/X
            /AU)*0.3/AU+0.55*std::exp(-1.2/X/AU)*1.2/AU+0.1*std::exp(-6.0/X/
            AU)*6.0/AU);
};

std::tuple<Velocity,Velocity> Ion::recoilEnergyAndVelocity() {
    // Initial parameter
    double P = std::sqrt(random()/(M_PI*std::pow(substrateDensity, 2.0/3.0)));
    double COLUMBIAVK = 0.0143992*charge*substrateCharge;
    double MC = mass*substrateMass/(mass+substrateMass);
    double INVLAB = std::sqrt(energy*2.0/mass);
    double EC = 0.5*MC*INVLAB*INVLAB;
    double AU = 0.8854*0.529/std::pow(std::sqrt(mass)+std::sqrt(substrateCharge), 2.0/3.0);
    double ELINHARD = EC*AU/COLUMBIAVK;
    // Find RMIN for different energy
    double AA = P*P;
    if(AA == 0) {
        AA = .00001;
    }
    double BB = COLUMBIAVK/EC;
    double CC = -1;
    double COLUMRMIN = 0.5/AA*(-BB+std::sqrt(BB*BB-4*AA*CC));
    int CALTIME = 1;
    double RMINTRY1;
    double RMINTRY2;
    do {
        RMINTRY1 = COLUMRMIN;
        double DV = std::abs(-DF(RMINTRY1, COLUMBIAVK, AU) /EC-2*P*P*RMINTRY1);
        if(std::abs(DV) < .000001) {
            DV = 0.1;
        }
        RMINTRY2 = RMINTRY1+(1.0-F(RMINTRY1, COLUMBIAVK, AU)/EC-P*P*
            RMINTRY1*RMINTRY1)/DV;
        COLUMRMIN = RMINTRY2;
        CALTIME = CALTIME + 1;
        if(CALTIME > 10000) {
            break;
        }
    } while (std::abs(RMINTRY2 - RMINTRY1) <= .00001);
    double RMIN = 0.5*(RMINTRY2+RMINTRY1);
    // Calculate deflection angle
    double RBIERSACK = 2.0*(EC-F(RMIN, COLUMBIAVK, AU))/RMIN*RMIN
        /DF(RMIN, COLUMBIAVK, AU);
    double BBIERSACK = P/AU;
    double ROBIERSACK = 1.0/(RMIN*AU);
    double RCBIERSACK = RBIERSACK/AU;
    double C1BIERSACK = 0.6743;
    double C2BIERSACK = 0.009611;
    double C3BIERSACK = 0.005175;
    double C4BIERSACK = 10.0;
    double C5BIERSACK = 6.314;
    double ALTHABIERSACK = 1.0+1.0/std::sqrt(C1BIERSACK*ELINHARD);
    double BELTABIERSACK = (C2BIERSACK+std::sqrt(ELINHARD))
        /(C3BIERSACK+std::sqrt(ELINHARD));
    double GAMABIERSACK = (C4BIERSACK+ELINHARD)/(C5BIERSACK+ELINHARD);
    double ABIERSACK = 2.0*ALTHABIERSACK*ELINHARD
        *std::pow(BBIERSACK, BELTABIERSACK);
    double GBIERSACK = GAMABIERSACK*1.0
        /(std::sqrt((1.0+ABIERSACK*ABIERSACK))-ABIERSACK);
    double DELTABIERSACK = ABIERSACK*(ROBIERSACK-BBIERSACK)/(1.0+GBIERSACK);
    double CALPHA1;
    if(P == 0.0) {
        CALPHA1 = M_PI;
    } else {
        CALPHA1 = 2.0*std::atan(std::sqrt((ROBIERSACK+RCBIERSACK)
            *(ROBIERSACK+RCBIERSACK)/((BBIERSACK+RCBIERSACK+DELTABIERSACK)
            *(BBIERSACK+RCBIERSACK+DELTABIERSACK))-1.0));
    }
    // Defelction angle in lab coordinator
    double COSPLUSMASS = std::cos(CALPHA1)+mass/substrateMass;
    double THETA1RELATIVE;
    if(COSPLUSMASS < 0.0) {
        THETA1RELATIVE = M_PI+std::atan(std::sin(CALPHA1)/COSPLUSMASS);
    } else if (COSPLUSMASS == 0.0) {
        THETA1RELATIVE = M_PI/2.0;
    } else {
        THETA1RELATIVE = std::atan(std::sin(CALPHA1)/COSPLUSMASS);
    }

    // Calculate target direction and energy
    double RE = 4.0*EC*MC/substrateMass*std::sin(CALPHA1/2.0)*std::sin(CALPHA1
                /2.0);
    double THETA2RELATIVE = (M_PI-CALPHA1)/2.0;

    Velocity newVelocity;
    Velocity targetVelocity;

    std::tie(newVelocity, targetVelocity) = 
        relativeToAbsoluteVelocity(THETA1RELATIVE, THETA2RELATIVE, RE);

    return std::make_tuple(newVelocity, targetVelocity);
};

double Ion::electronicStoppingEnergy() {
    double stoppingEnergy = 1.212*std::pow(charge, 7.0/6.0)*substrateCharge
                /(std::pow((std::pow(charge, 2.0/3.0) + std::pow(substrateCharge, 2.0/3.0))
                , 3.0/2.0)*std::sqrt(mass))*std::sqrt(energy*1000);
    double energyLoss = 1.59*atomicSpacing()*substrateDensity*stoppingEnergy/1000;
    return energyLoss;
};


// Substrate functions
Substrate::Substrate(Coordinate coordinate, Velocity velocity, double charge,
    double mass, double density, double range, std::weak_ptr<Bombardment> bombardment) :
    Ion(coordinate, velocity, charge, mass, mass, density, mass, ION, range,
        bombardment) {};

//  Electron functions


// Input fields functions
InputFields::InputFields() : charge(Defaults::ionCharge),
energy(Defaults::ionEnergy), mass(Defaults::ionMass),
simulationCount(Defaults::simulationCount),
substrateCharge(Defaults::substrateCharge),
substrateDensity(Defaults::substrateDensity),
substrateMass(Defaults::substrateMass), type(ION), range(Defaults::range) {};

InputFields::InputFields(double charge, double energy,
double mass, size_t simulationCount,
double substrateCharge, double substrateDensity,
double substrateMass, ParticleType type, double range)
: charge(charge), energy(energy), mass(mass), 
simulationCount(simulationCount),  substrateCharge(substrateCharge),
substrateDensity(substrateDensity), substrateMass(substrateMass),
type(type), range(range) {};

double InputFields::getCharge() const { return charge; };
double InputFields::getEnergy() const { return energy; };
double InputFields::getMass() const { return mass; };
size_t InputFields::getSimulationCount() const { return simulationCount; };
double InputFields::getSubstrateCharge() const { return substrateCharge; };
double InputFields::getSubstrateDensity() const { return substrateDensity; };
double InputFields::getSubstrateMass() const { return substrateMass; };
ParticleType InputFields::getType() const { return type; };
double InputFields::getRange() const { return range; };

// Bombardment functions
Bombardment::Bombardment() : particles(0) {
    DEBUG_PRINT("Bombardment created");
};

Bombardment::~Bombardment() {
    DEBUG_PRINT("Destructing bombardment...");
};

const std::vector<std::unique_ptr<Particle>>& Bombardment::getParticles() const {
    return particles;
};

void Bombardment::initiate(InputFields& inputFields) {
    DEBUG_PRINT("Initiating bombardment...");
    std::unique_ptr<Particle> initialParticle;
    Coordinate initialCoordinate = {0.0, 0.0, 0.0};
    Velocity initialVelocity = {0.0, 0.0, inputFields.getEnergy()};

    if(inputFields.getType() == ION) {
        DEBUG_PRINT("Making Ion...");
        initialParticle = std::make_unique<Ion>(initialCoordinate, initialVelocity,
            inputFields, shared_from_this());
    }

    DEBUG_PRINT("Call addParticle()");
    addParticle(std::move(initialParticle));
};

void Bombardment::addParticle(std::unique_ptr<Particle> particle) {
    DEBUG_PRINT("Pushing back particle...");
    particle->fire();
    particles.push_back(std::move(particle));
}

// Simulation functions
Simulation::Simulation() : inputs(InputFields()), 
    bombardments(Defaults::simulationCount) {
    DEBUG_PRINT("Creating simulation (Default)");
    for(std::shared_ptr<Bombardment> bombardment : bombardments) {
        bombardment = std::make_shared<Bombardment>();
    }
};

Simulation::Simulation(InputFields inputs)
: inputs(inputs), bombardments(inputs.getSimulationCount()) {
    DEBUG_PRINT("Creating simulation");
    for(std::shared_ptr<Bombardment>& bombardment : bombardments) {
        DEBUG_PRINT("Allocating bombardment");
        bombardment = std::make_shared<Bombardment>();
    }
};

Simulation::~Simulation() {
    DEBUG_PRINT("Destructing simulation");
};

InputFields Simulation::getInputs() { return inputs; };
size_t Simulation::getBombardmentSize() { return bombardments.capacity(); };
size_t Simulation::getBombardmentCapacity() { return bombardments.capacity(); };

void Simulation::initiate() {
    DEBUG_PRINT("Initiating simulation");
    for(std::shared_ptr<Bombardment> bombardment : bombardments) {
        DEBUG_PRINT("Initializing bombardment...");
        bombardment->initiate(inputs);
    };
};

void Simulation::writeCoordinateData() {
    DEBUG_PRINT("Writing coordinate data...");
    std::filesystem::path outputPath = fs::path(".") / "output";
    fs::create_directories(outputPath);
    std::string coordinateOutputFilename = "coordinateOutput.csv";

    std::ofstream coordinateOutputFile(outputPath / coordinateOutputFilename);
    if(coordinateOutputFile.is_open()) {
        coordinateOutputFile << "Bombardment,Particle,x,y,z" << std::endl;

        size_t bombardmentIndex = 0;
        size_t particleIndex = 0;

        for(std::shared_ptr<Bombardment>& bombardment : bombardments) {
            DEBUG_PRINT("Bombardment " << bombardmentIndex);
            for(const std::unique_ptr<Particle>& particle : bombardment->getParticles()) {
                DEBUG_PRINT("Particle " << particleIndex);
                DEBUG_PRINT("Size:\t" << particle->getCoordinates().size());
                for(const Coordinate coordinate : particle->getCoordinates()) {
                    coordinateOutputFile << bombardmentIndex << ","
                                         << particleIndex << ","
                                         << coordinate.x << ","
                                         << coordinate.y << ","
                                         << coordinate.z << std::endl;
                }
                particleIndex++;
            }
            particleIndex = 0;
            bombardmentIndex++;
        }

        coordinateOutputFile.close();
    } else {
        std::cerr << "Error: Unable to open " << coordinateOutputFilename 
                  << std::endl;
    };
};