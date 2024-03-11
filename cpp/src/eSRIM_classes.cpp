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
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <typeinfo>


// Local includes
#include "eSRIM_classes.h"

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
    // DEBUG_PRINT("Pushing back coordinate");
    coordinate_vector.emplace_back(coordinate);
};

const Coordinate& Particle::getCoordinate() const {
    return coordinate;
};

const std::vector<Coordinate>& Particle::getCoordinates() const {
    return coordinate_vector;
};

double Particle::random() {
    return randomDistribution(randomGenerator);
};

size_t Particle::getNumParticles() {
    if(auto locked = bombardment.lock()) {
        // DEBUG_PRINT("Bombardment lock");
        return locked->getParticles().size();
    } else {
        DEBUG_PRINT("Could not get bombardment size");
        return 0;
    }
}

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
    } else if (Z > 0.0) {
        THETA1 = std::atan(std::sqrt(X*X+Y*Y)/Z);
    } else {
        THETA1 = M_PI/2.0;
    } 

    double ALPHA1;
    if(std::sin(THETA1)!= 0.0) {
        if(X < 0.0) {
            ALPHA1 = M_PI+std::atan(Y/X); 
        } else if(X > 0.0) {
            ALPHA1 = std::atan(Y/X); 
        } else {
            ALPHA1 = M_PI - sign(Y)*M_PI/2.0; 
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
    } else if(Z > 0.0) {
        THETA2 = std::atan(std::sqrt(X*X+Y*Y)/Z);
    } else {
        THETA2 = M_PI/2;
        
    }

    double ALPHA2;
    if(std::sin(THETA2) != 0.0) {
        if(X < 0.0) {
            ALPHA2 = M_PI+std::atan(Y/X);
        } else  if(X > 0.0) {
            ALPHA2 = std::atan(Y/X);
        } else {
            ALPHA2 = M_PI - 0.5*sign(Y)*M_PI;
        }
    } else {
        ALPHA2 = 0.0;
    }
    Velocity newVelocity = {THETA1, ALPHA1, energy};
    Velocity targetVelocity = {THETA2, ALPHA2, targetEnergy};

    // if(std::isnan(THETA1) || std::isnan(THETA2) || std::isnan(ALPHA1) || std::isnan(ALPHA2)) {
    //     DEBUG_PRINT("zAngle: " << THETA2
    //                 << ", xAngle: " << ALPHA2);
    // }

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
    // DEBUG_PRINT("Ion made");
};

Ion::Ion(Coordinate coordinate, Velocity velocity, double charge,
    double mass, double substrateCharge, double substrateDensity,
    double substrateMass, ParticleType type, double range,
    std::weak_ptr<Bombardment> bombardment)
    : Particle(coordinate, velocity, charge, mass, substrateCharge, 
    substrateDensity, substrateMass, type, range, bombardment) {
    // DEBUG_PRINT("Ion made");
}; 

const size_t& Ion::getDepth() const { 
    DEBUG_PRINT("Ion depth function called");
    return coordinate.depth;
};

void Ion::fire() {
    if(velocity.energy <= Constants::ionStoppingEnergy) {
        // DEBUG_PRINT("Particle " << id  << " does not have enough energy");
        // addCoordinate(coordinate);
        return;
    } else {
        // DEBUG_PRINT("Firing particle " << coordinate.id);
        while(velocity.energy > Constants::ionStoppingEnergy) {
            // DEBUG_PRINT("Energy " << fireCount << ":\t" << velocity.energy);
            velocity.energy -= electronicStoppingEnergy();

            Velocity newVelocity;
            Velocity targetVelocity;

            std::tie(newVelocity, targetVelocity) = recoilEnergyAndVelocity();

            Coordinate newCoordinate = {
                    coordinate.x + atomicSpacing()*std::sin(velocity.zAngle)
                                    *std::cos(velocity.xAngle),
                    coordinate.y + atomicSpacing()*std::sin(velocity.zAngle)
                                    *std::sin(velocity.xAngle),
                    coordinate.z + atomicSpacing()*std::cos(velocity.zAngle),
                    coordinate.depth
                };

            if(newCoordinate.z < 0) {
                DEBUG_PRINT("Sputter")
                return;
            }

            if(Defaults::enableDamageCascade) {
                if(energy - targetVelocity.energy > Constants::ionDisplacementEnergy) {
                    // DEBUG_PRINT("Attempting lock...");
                    Coordinate targetCoordinate = newCoordinate;
                    targetCoordinate.depth++;
                    if(auto locked = bombardment.lock()) {    
                        std::unique_ptr<Substrate> targetParticle;
                        targetParticle = std::make_unique<Substrate>(targetCoordinate,
                            targetVelocity, substrateCharge, substrateMass, 
                            substrateDensity, range, bombardment);
                        locked->addParticle(std::move(targetParticle));
                    } else {
                        DEBUG_PRINT("Failed to create knock-on particle");
                    }
                    // DEBUG_PRINT("Knock-off substrate ion created");
                }
            }

            // DEBUG_PRINT("x, y, z: " << newCoordinate.x << ", "<< newCoordinate.y << ", "<< newCoordinate.z );
            addCoordinate(coordinate);
            coordinate = newCoordinate;
            velocity = newVelocity;
        }
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
    if(P != 0.0) {
        CALPHA1 = 2.0*std::atan(std::sqrt((ROBIERSACK+RCBIERSACK)
            *(ROBIERSACK+RCBIERSACK)/((BBIERSACK+RCBIERSACK+DELTABIERSACK)
            *(BBIERSACK+RCBIERSACK+DELTABIERSACK))-1.0));
    } else {
        CALPHA1 = M_PI;
    }
    // Defelction angle in lab coordinator
    double COSPLUSMASS = std::cos(CALPHA1)+mass/substrateMass;
    double THETA1RELATIVE;
    if(COSPLUSMASS > 0.0) {
        THETA1RELATIVE = std::atan(std::sin(CALPHA1)/COSPLUSMASS);
    } else if (COSPLUSMASS < 0.0) {
        THETA1RELATIVE = M_PI+std::atan(std::sin(CALPHA1)/COSPLUSMASS);
    } else {
        THETA1RELATIVE = M_PI/2.0;
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
    double mass, double density, double range,
    std::weak_ptr<Bombardment> bombardment) :
    Ion(coordinate, velocity, charge, mass, mass, density, mass, SUBSTRATE, range,
        bombardment) {};

const size_t& Substrate::getDepth() const { return depth; };

//  Electron functions
Electron::Electron(Coordinate initialCoordinate, Velocity initialVelocity, InputFields input,
        std::weak_ptr<Bombardment> bombardment)
        : Particle(initialCoordinate, initialVelocity,
        Constants::electronCharge, Constants::electronMass,
        input.getSubstrateCharge(), input.getSubstrateDensity(),
        input.getSubstrateMass(), ELECTRON, input.getRange(),
        bombardment) {
    DEBUG_PRINT("Electron made");
};

void Electron::readParameters() {
    DEBUG_PRINT("Reading Mott scattering values");
};

void Electron::fire() {
    DEBUG_PRINT("Firing electron!");
    // size_t numParticles = getNumParticles();
    while(velocity.energy > Constants::electronStoppingEnergy) {
        DEBUG_PRINT("Energy: " << energy << " keV");
        
    }   
};

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
Bombardment::Bombardment(std::weak_ptr<Simulation> simulation, size_t id)
    : particles(0), simulation(simulation), id(id) {
    DEBUG_PRINT("Bombardment created");
};

Bombardment::~Bombardment() {
    DEBUG_PRINT("Bombardment " << id << " Destroyed")
};

const std::vector<std::unique_ptr<Particle>>& Bombardment::getParticles() const {
    return particles;
};

void Bombardment::initiate(InputFields& inputFields) {
    DEBUG_PRINT("Initiating bombardment " << id);
    std::unique_ptr<Particle> initialParticle;
    Coordinate initialCoordinate = {0.0, 0.0, 0.0, 0};
    Velocity initialVelocity = {0.0, 0.0, inputFields.getEnergy()};

    if(inputFields.getType() == ION) {
        // DEBUG_PRINT("Making Ion...");
        initialParticle = std::make_unique<Ion>(initialCoordinate, initialVelocity,
            inputFields, shared_from_this());
    } else {
        initialParticle = std::make_unique<Electron>(initialCoordinate, initialVelocity, 
            inputFields, shared_from_this());
    }

    addParticle(std::move(initialParticle));

    if(auto locked = simulation.lock()) {
        locked->writeData(COORDINATE, getParticles(), id);
    } else {
        std::cerr << "Failed to initiate write" << std::endl;
    }
};

void Bombardment::addParticle(std::unique_ptr<Particle> particle) {
    // DEBUG_PRINT("Pushing back particle...");
    particles.emplace_back(std::move(particle));
    particles.back()->fire();
}

// Simulation functions
// Simulation::Simulation() : inputs(InputFields()), 
//     bombardments(Defaults::simulationCount) {
//     DEBUG_PRINT("Creating simulation (Default)");
//     for(std::shared_ptr<Bombardment> bombardment : bombardments) {
//         bombardment = std::make_shared<Bombardment>();
//     }
// };

void Bombardment::popParticle() {
    particles.pop_back();
};

Simulation::Simulation(InputFields inputs)
: inputs(inputs), bombardments(inputs.getSimulationCount()),
  outputPath(fs::path(".") / Defaults::outputDirectory) {
    DEBUG_PRINT("Creating simulation");
    // for(std::shared_ptr<Bombardment>& bombardment : bombardments) {
};

Simulation::~Simulation() {
    DEBUG_PRINT("Destructing simulation");
};

InputFields Simulation::getInputs() { return inputs; };
size_t Simulation::getBombardmentSize() { return bombardments.capacity(); };
size_t Simulation::getBombardmentCapacity() { return bombardments.capacity(); };

void Simulation::initiate() {

    checkOutputFiles();

    for(size_t i = 0; i < inputs.getSimulationCount(); i++) {
        // DEBUG_PRINT("Allocating bombardment");
        bombardments.at(i) = std::make_shared<Bombardment>(shared_from_this(), i);
    }

    const size_t simulation_count = inputs.getSimulationCount();
    threads.reserve(simulation_count); // Reserve space to avoid unnecessary reallocations

    DEBUG_PRINT("Initiating simulation");
    for(size_t i = 0; i < simulation_count; i++) {
        // DEBUG_PRINT("Initializing bombardment...");
        threads.emplace_back([this, i]() {
            bombardments.at(i)->initiate(inputs);
            // std::cout << "Use before: " << bombardments.at(i).use_count() << std::endl;
            bombardments.at(i).reset();
            // std::cout << "Use after: " << bombardments.at(i).use_count() << std::endl;
        });
    }

    for(auto& thread : threads) {
        if (thread.joinable()) {
            thread.join();
        }
    }

    writeData(ENDOFFILE, {}, 0);
}

void Simulation::writeData(OutputType outputType,
    const std::vector<std::unique_ptr<Particle>>& particles = {},
    size_t simulationID = 0) {
    
    // Acquire lock
    std::lock_guard<std::mutex> lock(fileLock);

    DEBUG_PRINT("Simulation " << simulationID << " Writing");

    // Check if output directory exists, create it if not
    if (!fs::exists(outputPath)) {
        fs::create_directories(outputPath);
    }

    // Open output file
    std::string filename = std::string(Defaults::outputCoordinateFilename);
    filename = filename + std::string(Defaults::outputFileExtension);
    std::ofstream outputFile(outputPath / filename, std::ios_base::app);
    if (!outputFile.is_open()) {
        std::cerr << "Failed to open output file." << std::endl;
        return;
    }

    // Write headers if the file is empty
    if (outputFile.tellp() == 0) {
        outputFile << "simulationID,particleID,depth,x,y,z\n";
    }

    size_t particleID = 0;

    // Write data to file
    switch (outputType) {
        case OutputType::COORDINATE:
            for (const auto& particle : particles){
                if(particle->getCoordinates().empty()) {
                    if(Defaults::logSingleDisplacement)
                    {
                        outputFile << simulationID << ','
                                << particleID << ','
                                << particle->getCoordinate().depth << ','
                                << particle->getCoordinate().x << ','
                                << particle->getCoordinate().y << ','
                                << particle->getCoordinate().z << '\n';
                        particleID++;
                    } else {
                        continue;
                    }
                } else {
                    for (const auto& coordinate : particle->getCoordinates()) {
                        // DEBUG_PRINT("Writing coordinate " << coordinate.id);
                        outputFile << simulationID << ','
                                    << particleID << ','
                                    << coordinate.depth << ','
                                    << coordinate.x << ','
                                    << coordinate.y << ','
                                    << coordinate.z << '\n';
                    }
                    particleID++;
                }
            }
            break;
        
        // Add cases for other OutputType values as needed
        case OutputType::ENDOFFILE:
            outputFile << Defaults::outputFileEnd;
            break; 
        default:
            std::cerr << "Unsupported output type." << std::endl;
            break;
    }

    // Close output file
    outputFile.close();
};

bool Simulation::fileIsWritten(const fs::path filename) {
    
    const int bufferSize = strlen(Defaults::outputFileEnd);

    std::ifstream file(filename, std::ios::ate);
    if (!file.is_open()) {
        return false;
    }
    // Check if file size is smaller than the phrase length
    if (file.tellg() < bufferSize) {
        return false;
    }
    
    // Move file pointer to bufferSize characters before the end
    file.seekg(-bufferSize, std::ios::end);
    // Read the last bufferSize characters
    std::string lastChars(bufferSize, '\0');
    file.read(&lastChars[0], bufferSize);
    return lastChars == std::string(Defaults::outputFileEnd);
};

void Simulation::renameFileWithTimestamp(const std::string& filename) {
    std::string timestamp = getCurrentDateTime();
    std::string newFilename = filename + "_" + timestamp
        + Defaults::outputFileExtension;
    fs::path newFilePath = outputPath/newFilename;
    std::string oldFilename = filename + Defaults::outputFileExtension;
    fs::path filePath = outputPath/oldFilename;
    fs::rename(filePath, newFilePath);
}

std::string Simulation::getCurrentDateTime() {
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    std::string timestamp = std::ctime(&now_c);
    timestamp.pop_back();
    std::replace(timestamp.begin(), timestamp.end(), ' ', '_');
    std::replace(timestamp.begin(), timestamp.end(), ':', '-');
    return timestamp;
};

void Simulation::checkOutputFiles() {
    fs::path filepath;
    std::string filename;
    for(int i = COORDINATE; i  < ENDOFFILE; i++) {
        switch(i) {
            case COORDINATE:
                filename = std::string(Defaults::outputCoordinateFilename);
                filename = filename + std::string(Defaults::outputFileExtension);
                filepath = outputPath/filename;
                if(fileIsWritten(filepath)) {
                    filename = std::string(Defaults::outputCoordinateFilename);
                    renameFileWithTimestamp(filename);
                }
                break;
            default:
                break;
        }
    }
};