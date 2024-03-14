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
#ifndef ESRIM_CLASSES_CXX
#define ESRIM_CLASSES_CXX

// Header files
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <typeinfo>


// Local includes
#include "eSRIM_classes.h"
#include "utilities.h"

// Particle functions and static member initializers
std::mt19937 Particle::randomGenerator;
std::uniform_real_distribution<double> Particle::randomDistribution(
    0+std::numeric_limits<double>::epsilon(),
    1-std::numeric_limits<double>::epsilon());
std::mutex Particle::randomGeneratorLock;

void Particle::seedRandomGenerator() {
    std::lock_guard<std::mutex> lock(randomGeneratorLock);
    randomGenerator.seed(
        std::chrono::system_clock::now().time_since_epoch().count());
};

Particle::Particle(
    const Coordinate& coordinate,
    const Velocity& velocity,
    std::shared_ptr<InputFields> input,
    std::weak_ptr<Bombardment> bombardment)
    : coordinate(coordinate),
    velocity(velocity),
    enableDamageCascade(input->getEnableDamageCascade()),
    charge(input->getCharge()),
    mass(input->getMass()),
    substrateCharge(input->getSubstrateCharge()),
    substrateDensity(input->getSubstrateDensity()),
    substrateMass(input->getSubstrateMass()),
    type(input->getType()),
    range(input->getRange()),
    ionDisplacementEnergy(input->getIonDisplacementEnergy()),
    ionStoppingEnergy(input->getIonStoppingEnergy()),
    input(input),
    bombardment(bombardment) {};

Particle::Particle(
    const Coordinate& coordinate,
    const Velocity& velocity,
    bool enableDamageCascade,
    double charge,
    double mass,
    double substrateCharge,
    double substrateDensity,
    double substrateMass,
    ParticleType type,
    double range,
    double ionDisplacementEnergy,
    double ionStoppingEnergy,
    std::shared_ptr<InputFields> input,
    std::weak_ptr<Bombardment> bombardment)
    : coordinate(coordinate),
    velocity(velocity),
    enableDamageCascade(enableDamageCascade),
    charge(charge),
    mass(mass),
    substrateCharge(substrateCharge),
    substrateDensity(substrateDensity),
    substrateMass(substrateMass),
    type(type),
    range(range),
    ionDisplacementEnergy(ionDisplacementEnergy),
    ionStoppingEnergy(ionStoppingEnergy),
    input(input),
    bombardment(bombardment) {};

Particle::~Particle() {};

void Particle::addCoordinate(const Coordinate& coordinate) {
    coordinate_vector.emplace_back(coordinate);
};

const Coordinate& Particle::getCoordinate() const {
    return coordinate;
};

const std::vector<Coordinate>& Particle::getCoordinates() const {
    return coordinate_vector;
};

double Particle::random() {
    std::lock_guard<std::mutex> lock(randomGeneratorLock);
    return randomDistribution(randomGenerator);
};

size_t Particle::getNumParticles() const {
    if(auto locked = bombardment.lock()) {
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

inline double Particle::atomicSpacing() const {
    return std::pow(substrateDensity, -1.0/3.0);
};

void Particle::relativeToAbsoluteVelocity (
        Velocity& newVelocity,
        Velocity& targetVelocity) const noexcept {
    const double& THETA1RELATIVE = newVelocity.zAngle;
    const double& THETA2RELATIVE = targetVelocity.zAngle;

    const double ALPHA1RELATIVE = random()*2.0*M_PI;
    const double ALPHA2RELATIVE = ALPHA1RELATIVE + M_PI;

    const double& THETAO = velocity.zAngle;
    const double& ALPHAO = velocity.xAngle;

    // Angle
    const double X1 = std::sin(THETA1RELATIVE) * std::cos(ALPHA1RELATIVE);
    const double Y1 = std::sin(THETA1RELATIVE) * std::sin(ALPHA1RELATIVE);
    const double Z1 = std::cos(THETA1RELATIVE);

    const double Y0 = Y1 * std::cos(THETAO) + Z1*std::sin(THETAO);
    const double Z0 = -Y1*std::sin(THETAO)+Z1*std::cos(THETAO);
    const double X0 = X1;

    const double Z = Z0;
    const double X = X0*std::sin(ALPHAO)+Y0*std::cos(ALPHAO);
    const double Y = -X0*std::cos(ALPHAO)+Y0*std::sin(ALPHAO);

    const double THETA1 = 
        (Z > 0.0) ? std::atan(std::sqrt(X*X+Y*Y)/Z)
        : (Z < 0.0) ? M_PI+std::atan(std::sqrt(X*X + Y*Y)/Z) 
        : M_PI/2.0;

    const double ALPHA1 = 
        (std::sin(THETA1)!= 0.0) ? 
            ((X < 0.0) ? M_PI+std::atan(Y/X)
            : (X > 0.0) ? std::atan(Y/X) 
            : M_PI - sign(Y)*M_PI/2.0)
        : 0.0;

    // Target
    const double X3 = std::sin(THETA2RELATIVE)*std::cos(ALPHA2RELATIVE);
    const double Y3 = std::sin(THETA2RELATIVE)*std::sin(ALPHA2RELATIVE);
    const double Z3 = std::cos(THETA2RELATIVE);

    const double Y2 = Y3*std::cos(THETAO)+Z3*std::sin(THETAO);
    const double Z2 = -Y3*std::sin(THETAO)+Z3*std::cos(THETAO);
    const double X2 = X3;

    const double Z5 = Z2;
    const double X5 = X2 * std::sin(ALPHAO) + Y2 * std::cos(ALPHAO);
    const double Y5 = -X2 * std::cos(ALPHAO) + Y2 * std::sin(ALPHAO);

    const double THETA2 = 
        (Z5 < 0.0) ? M_PI+std::atan(std::sqrt(X5*X5+Y5*Y5)/Z5)
        : (Z5 > 0.0) ? std::atan(std::sqrt(X5*X5+Y5*Y5)/Z5)
        : M_PI/2;

    const double ALPHA2 = 
        (std::sin(THETA2) != 0.0) ?
            ((X5 < 0.0) ? M_PI+std::atan(Y5/X5)
            : (X5 > 0.0) ? std::atan(Y5/X5)
            : M_PI - 0.5*sign(Y5)*M_PI)
        : 0.0;
    newVelocity.zAngle = THETA1;
    newVelocity.xAngle = ALPHA1;
    targetVelocity.zAngle = THETA2;
    targetVelocity.xAngle = ALPHA2;
};

inline double Particle::sign(const double x) const noexcept{
    return (x > 0.0) ? 1.0 : (x < 0.0) ? -1.0 : 0.0;
};

inline Coordinate Particle::calculateNewCoordinate() {
    throw std::logic_error(
        "Base class Particle::calculateNewCoordinate() called");
    Coordinate dummy;
    return dummy;
}

//  Ion functions
Ion::Ion(
    const Coordinate& coordinate,
    const Velocity& velocity,
    std::shared_ptr<InputFields> input,
    std::weak_ptr<Bombardment> bombardment)
    : Particle(coordinate, velocity, input, bombardment) {};

Ion::Ion(
    const Coordinate& coordinate,
    const Velocity& velocity,       
    bool enableDamageCascade,
    double charge,
    double mass,
    double substrateCharge,
    double substrateDensity,
    double substrateMass,
    ParticleType type,
    double range,
    double ionDisplacementEnergy,
    double ionStoppingEnergy,
    std::shared_ptr<InputFields> input,
    std::weak_ptr<Bombardment> bombardment)
    : Particle(
        coordinate,
        velocity,
        enableDamageCascade,
        charge,
        mass,
        substrateCharge, 
        substrateDensity,
        substrateMass,
        type,
        range,
        ionDisplacementEnergy,
        ionStoppingEnergy,
        input,
        bombardment) {}; 

const size_t& Ion::getDepth() const { 
    DEBUG_PRINT("Ion depth function called");
    return coordinate.depth;
};

void Particle::createSubstrateKnockon(
        const Coordinate& newCoordinate,
        const Velocity& targetVelocity) {
    
        Coordinate targetCoordinate = newCoordinate;
        targetCoordinate.depth++;
        if(auto locked = bombardment.lock()) {    
            std::unique_ptr<Substrate> targetParticle;
            targetParticle = std::make_unique<Substrate>(
                targetCoordinate,
                targetVelocity,
                input,
                bombardment);
            locked->addParticle(std::move(targetParticle));
        } else {
            DEBUG_PRINT("Failed to create knock-on particle");
        }
    
};

void Ion::fire() {
    Velocity targetVelocity;
    Velocity newVelocity;
    if(energy < input->getIonStoppingEnergy()
            && !input->getLogSingleDisplacement()) {
        return;
    }
    while(energy > ionStoppingEnergy) {
        velocity.energy -= electronicStoppingEnergy();
        recoilEnergyAndVelocity(newVelocity, targetVelocity);
        velocity = newVelocity;
        coordinate = calculateNewCoordinate();
        if(coordinate.z < 0.0) {
            DEBUG_PRINT("Ion Sputter");
            return;
        }
        if(enableDamageCascade
                && energy - targetVelocity.energy > ionDisplacementEnergy) {
            createSubstrateKnockon(coordinate, targetVelocity);
        }

        if(!input->getLogStoppingPointOnly()) {
            addCoordinate(coordinate);
        }
    }
    if(input->getLogStoppingPointOnly()) {
        addCoordinate(coordinate);
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

void Ion::recoilEnergyAndVelocity(
        Velocity& newVelocity,
        Velocity& targetVelocity) {
    // Initial parameter
    double P = std::sqrt(random()/(M_PI*std::pow(substrateDensity, 2.0/3.0)));
    double COLUMBIAVK = 0.0143992*charge*substrateCharge;
    double MC = mass*substrateMass/(mass+substrateMass);
    double INVLAB = std::sqrt(energy*2.0/mass);
    double EC = 0.5*MC*INVLAB*INVLAB;
    double AU = 0.8854*0.529/std::pow(std::sqrt(mass)
                +std::sqrt(substrateCharge), 2.0/3.0);
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

    newVelocity.zAngle = THETA1RELATIVE;
    newVelocity.energy = energy - RE;
    targetVelocity.zAngle = THETA2RELATIVE;
    targetVelocity.energy = RE;
    
    relativeToAbsoluteVelocity(newVelocity, targetVelocity);
};

double Ion::electronicStoppingEnergy() {
    double stoppingEnergy = 1.212*std::pow(charge, 7.0/6.0)*substrateCharge
                /(std::pow((std::pow(charge, 2.0/3.0) \
                + std::pow(substrateCharge, 2.0/3.0))
                , 3.0/2.0)*std::sqrt(mass))*std::sqrt(energy*1000);
    double energyLoss = 1.59*atomicSpacing()*substrateDensity*stoppingEnergy/1000;
    return energyLoss;
};

inline Coordinate Ion::calculateNewCoordinate() {
    return {
        coordinate.x + atomicSpacing()*std::sin(velocity.zAngle)
                        *std::cos(velocity.xAngle),
        coordinate.y + atomicSpacing()*std::sin(velocity.zAngle)
                        *std::sin(velocity.xAngle),
        coordinate.z + atomicSpacing()*std::cos(velocity.zAngle),
        coordinate.depth
    };
};

// Substrate functions
Substrate::Substrate(
    const Coordinate& coordinate,
    const Velocity& velocity,
    std::shared_ptr<InputFields> input,
    std::weak_ptr<Bombardment> bombardment)
    : Ion(coordinate,
        velocity,
        input->getEnableDamageCascade(),
        input->getSubstrateCharge(),
        input->getSubstrateMass(),   
        input->getSubstrateCharge(),
        input->getSubstrateDensity(),
        input->getSubstrateMass(),
        SUBSTRATE,
        input->getRange(),
        input->getIonDisplacementEnergy(),
        input->getIonStoppingEnergy(),
        input,
        bombardment) {};

const size_t& Substrate::getDepth() const { return depth; };

//  Electron functions
std::vector<std::array<std::array<
    double,
    Constants::numMottkParam>,
    Constants::numMottjParam>> 
    Electron::mottScatteringParams;

std::vector<std::array<
    double,
    Constants::numElecScreeningPotentialParams>>
    Electron::elecScreeningParams;

std::vector<double> Electron::divisionAngles;

Electron::Electron(
    const Coordinate& initialCoordinate,
    const Velocity& initialVelocity,
    std::shared_ptr<InputFields> input,
    std::weak_ptr<Bombardment> bombardment)
    : Particle(
        initialCoordinate,
        initialVelocity,
        input,
        bombardment),
    randomCrossSections(input->getNumFlyingDistances()),
    flyingDistances(input->getNumFlyingDistances()),
    scatteringAngles(input->getNumFlyingDistances()),
    partialCrossSections(input->getNumAngleDivisors()),
    elasticEnergyLoss(input->getNumFlyingDistances()),
    correctionFactor(1.0),
    electronStoppingEnergy(input->getElectronStoppingEnergy()) {};

inline double Electron::findAngle(size_t i, std::shared_ptr<InputFields> input) {
    // Value is offset from 1 to avoid checking a zero-angle
    return -M_PI/(9.0)*(1-std::pow(10.0,(static_cast<double>(i+1)
        /input->getNumAngleDivisors())));   
};

void Electron::readParametersAndInitialize(
        std::shared_ptr<InputFields> inputFields) {
    DEBUG_PRINT("Reading Mott scattering values");
    fs::path currentWorkingDirectory = fs::current_path();
    fs::path inputDirectory(inputFields->getInputDirectoryName());
    fs::path mottParametersPath = 
        currentWorkingDirectory
        / inputDirectory
        / inputFields->getMottScatteringParametersFilename();

    std::string input;
    std::string parameter;
    std::istringstream stringStream;

    std::ifstream mottParametersFile(mottParametersPath);

    if(!mottParametersFile.is_open()) {
        throw std::ios_base::failure("Could not open file "
            +mottParametersPath.string() );
    }
    
    mottScatteringParams.resize(inputFields->getNumMottScatteringPotentialElements());

    for(size_t i = 0; i < inputFields->getNumMottScatteringPotentialElements(); i++) {
        for(size_t j = 0; j < Constants::numMottjParam; j++) {
            std::getline(mottParametersFile, input);
            stringStream.clear();
            stringStream.str(input);

            for(size_t k = 0; k < Constants::numMottkParam; k++) {
                std::getline(stringStream, parameter, ',');
                mottScatteringParams[i][j][k] = std::stod(parameter);
            }    
            
        }
    }

    mottParametersFile.close();

    fs::path elecParametersPath = 
        currentWorkingDirectory
        / inputDirectory
        / inputFields->getElectronScreeningParametersFilename();


    std::ifstream elecParametersFile(elecParametersPath);

    if(!elecParametersFile.is_open()) {
        throw std::ios_base::failure("Could not open file "
            +elecParametersPath.string());
    }

    elecScreeningParams.resize(inputFields->getNumElecScreeningPotentialElements());

    for(size_t i = 0; i < inputFields->getNumElecScreeningPotentialElements(); i++) {
        std::getline(elecParametersFile, input);
        stringStream.clear();
        stringStream.str(input);
        for(size_t j = 0; j < Constants::numElecScreeningPotentialParams; j++){
            std::getline(stringStream, parameter, ',');
            elecScreeningParams[i][j] = std::stod(parameter);
        }
    }

    elecParametersFile.close();

    DEBUG_PRINT("Calculating division angles");

    divisionAngles.reserve(inputFields->getNumAngleDivisors());

    for(size_t i = 0; i < inputFields->getNumAngleDivisors(); i++) {
        divisionAngles.push_back(findAngle(i, inputFields));
    }
    
};

void Electron::fire() { 
    DEBUG_PRINT("Firing electron!");
    double totalMottCrossSection;
    double totalElasticEnergyLoss;
    double ionizationEnergyLoss;
    double bremsstrahlung;
    double newEnergy;
    bool firstIteration = true;
    double flightGroupLength;
    Velocity newVelocity = {0.0,0.0,0.0};
    Velocity targetVelocity = {0.0,0.0,0.0};
    while(velocity.energy > electronStoppingEnergy) {
        totalMottCrossSection = getMottTotalCrossSection();
        flightGroupLength = generateFlyingDistances(totalMottCrossSection);
        generateRandomCrossSections();
        generateIntegratedCrossSections(totalMottCrossSection);
        totalElasticEnergyLoss = getElasticEnergyLoss();
        ionizationEnergyLoss = flightGroupLength*getIonizationEnergyLoss();
        bremsstrahlung = flightGroupLength*getBremsstrahlung();
        newEnergy = energy - totalElasticEnergyLoss - ionizationEnergyLoss -
            bremsstrahlung;

        // TODO: Include logic for explicit and implicit methods
        if(firstIteration) {
            correctionFactor = (energy + newEnergy)/2.0/energy;
            firstIteration = false;
            continue;
        } 
        correctionFactor = (energy + newEnergy)/2.0/energy;
        newVelocity.energy = newEnergy;

        for(size_t i = 0; i < input->getNumFlyingDistances(); i++) {

            newVelocity.zAngle = scatteringAngles[i];
            targetVelocity.zAngle = (M_PI-scatteringAngles[i])/2.0;
            targetVelocity.energy = elasticEnergyLoss[i];

            relativeToAbsoluteVelocity(newVelocity,targetVelocity);
            velocity = newVelocity;
            coordinate = calculateNewCoordinate(i);

            if(coordinate.z < 0.0) {
                DEBUG_PRINT("Electron Sputter");
                return;
            }

            if(input->getEnableDamageCascade()
                && targetVelocity.energy > ionDisplacementEnergy) {
                createSubstrateKnockon(coordinate, targetVelocity);
            }

            if(!input->getLogEndOfFlyingDistanceOnly() 
                    && !input->getLogStoppingPointOnly()) {
                addCoordinate(coordinate);
            }
            
        }
        if(input->getLogEndOfFlyingDistanceOnly()
                && !input->getLogStoppingPointOnly()) {
            addCoordinate(coordinate);
        }
    }   
    if(input->getLogStoppingPointOnly()) {
        addCoordinate(coordinate);
    }
};

double Electron::getMottDifferentialCrossSection(double Theta) {
    const double correctedEnergy = energy*correctionFactor;
    double beta = std::sqrt(1.0-1.0/((correctedEnergy/Constants::mec2 + 1.0)
        *(correctedEnergy/Constants::mec2 + 1.0)));
    const double betaSquared = beta*beta;
    constexpr double betaAVE = 0.7181287;
    const double gamma = correctedEnergy/Constants::mec2+1.0;
    constexpr double elecmass = 1.0; // specific for a.u. unit,  only for Feq calcu.
    constexpr double Vc = 137.0;  //speed of light for a.u. unit, only for Feq calcu.
    constexpr double re1 = 2.817938E-13;  // in the unit of cm)
    const double oneMinusCos = 1.0-std::cos(Theta);
    const double betaQuarted = betaSquared*betaSquared;
    const double DSC = (substrateCharge*re1)*(substrateCharge*re1)*(1.0-betaSquared)/
        betaQuarted/(oneMinusCos*oneMinusCos);
    const double moment = 2.0*beta*gamma*elecmass*Vc*std::sin(Theta/2.0);
    const double momentSquared = moment*moment;
    
    
    double Rmott = 0.0;
    double alpha1;
    for(size_t j = 0; j < Constants::numMottjParam; j++) {
        alpha1 = 0.0;
        for(size_t k = 0; k < Constants::numMottkParam; k++) {
            alpha1 = alpha1 +
                mottScatteringParams[static_cast<size_t>(substrateCharge)][j][k]
                *std::pow(beta-betaAVE, k);
        }
        Rmott = Rmott+alpha1*std::pow(oneMinusCos,j/2.0);
    }
    
    double Feq = 0;
    for(size_t i = 0; i < Constants::numElecScreeningPotentialParams/2; i++) {
        Feq = Feq + elecScreeningParams[substrateCharge][i]*
            elecScreeningParams[substrateCharge][i + 3]
            *elecScreeningParams[substrateCharge][i + 3]
            /(elecScreeningParams[substrateCharge][i + 3]
            *elecScreeningParams[substrateCharge][i + 3]
            +momentSquared);
    }

    return Rmott*DSC*(1-Feq)*(1-Feq);
};

void Electron::generateRandomCrossSections() {
    for(size_t i = 0; i < input->getNumFlyingDistances(); i++) {
        randomCrossSections[i] = random(); 
    }
    std::sort(randomCrossSections.begin(), randomCrossSections.end());
};

double Electron::generateFlyingDistances(double totalMottCrossSection) {
    double distance;
    double flightGroupLength = 0;
    double crossSectionFlight = 0;
    double atomicSpacing = std::pow(substrateDensity*1e24, -1.0/3.0);

    for(size_t i = 0; i < input->getNumFlyingDistances(); i++) {
        crossSectionFlight = -std::log(random())
            /(substrateDensity*1e24*totalMottCrossSection);
        distance = crossSectionFlight+atomicSpacing;
        flightGroupLength += distance;
        flyingDistances[i] = distance; 
    }
    return flightGroupLength;
};

double Electron::getMottTotalCrossSection() {
    double totalCrossSection = 0.0;
    double differentialCrossSection =
        getMottDifferentialCrossSection(divisionAngles[0]);
    totalCrossSection += 
        differentialCrossSection*2.0*M_PI
        *std::sin(divisionAngles[0])
        *(divisionAngles[1]-divisionAngles[0])/2.0;
    partialCrossSections[0] = totalCrossSection;
    for(size_t i = 1; i < input->getNumAngleDivisors() - 1; i++) {
        differentialCrossSection =
            getMottDifferentialCrossSection(divisionAngles[i]);
        totalCrossSection += 
            differentialCrossSection*2.0
            *M_PI*std::sin(divisionAngles[i])
            *(divisionAngles[i+1] - divisionAngles[i-1])/2.0;
        partialCrossSections[i] = totalCrossSection;
    }
    totalCrossSection +=
        differentialCrossSection*2.0*M_PI*
        std::sin(divisionAngles[input->getNumAngleDivisors()-1])
        *(M_PI
            -(  divisionAngles[input->getNumAngleDivisors()-3]
                -divisionAngles[input->getNumAngleDivisors()-2]
            )/2.0);
    partialCrossSections[input->getNumAngleDivisors() - 1] = totalCrossSection;
    return totalCrossSection;
};

void Electron::generateIntegratedCrossSections(double totalMottCrossSection) {
    double prevFractionalCrossSection;
    double fractionalCrossSection;
    for(size_t i = 0; i < input->getNumFlyingDistances(); i++) {
        prevFractionalCrossSection = 0;
        for(size_t j = 0; j < input->getNumAngleDivisors(); j++) {
            fractionalCrossSection =
                partialCrossSections[j]
                /totalMottCrossSection;
            if(fractionalCrossSection > randomCrossSections[i]) {
                if(j > 0 && j < input->getNumAngleDivisors()-1) {
                    scatteringAngles[i] = 
                        (divisionAngles[j] + divisionAngles[j + 1])/2.0
                        -(fractionalCrossSection-randomCrossSections[i])
                        *(divisionAngles[j+1]-divisionAngles[j-1])/2.0
                        /(fractionalCrossSection - prevFractionalCrossSection);
                    prevFractionalCrossSection = fractionalCrossSection;
                } else if (j == 0) {
                    scatteringAngles[i] = 
                        (divisionAngles[j]+divisionAngles[j+1])/2.0
                        -(fractionalCrossSection-randomCrossSections[i])
                        *(divisionAngles[j] + divisionAngles[j+1])/2.0
                        /(fractionalCrossSection - prevFractionalCrossSection);
                    prevFractionalCrossSection = fractionalCrossSection;
                } else {
                    scatteringAngles[i] = M_PI
                        -(fractionalCrossSection-randomCrossSections[i])
                        *(M_PI - (divisionAngles[j-2]+divisionAngles[j-3])/2.0)
                        /(fractionalCrossSection - prevFractionalCrossSection);
                    prevFractionalCrossSection = fractionalCrossSection;
                }
                break;
            } else {
                continue;
            }
        }
    }
    std::shuffle(scatteringAngles.begin(), scatteringAngles.end(), randomGenerator);
};

double Electron::getElasticEnergyLoss() {
    double totalElasticEnergyLoss = 0;
    for(size_t i = 0; i < input->getNumFlyingDistances(); i++) {
        double ElecREup = 
            ((energy+511.0)*(std::sin(scatteringAngles[i]))
                *(std::sin(scatteringAngles[i]))+substrateMass*931.0*1000.0*(1
                -std::cos(scatteringAngles[i])))
            *energy*(energy+2.0*511.0);
        double ElecREdown = (energy+substrateMass*931.5*1000.0)
            *(energy+substrateMass*931.5*1000.0)-energy*(energy+2.0*511.0)
            *(std::cos(scatteringAngles[i]))*(std::cos(scatteringAngles[i]));
        elasticEnergyLoss[i] = ElecREup/ElecREdown;
        totalElasticEnergyLoss += elasticEnergyLoss[i];
    }
    return totalElasticEnergyLoss;
};

double Electron::getIonizationEnergyLoss() {
    double correctedEnergy = energy*correctionFactor;
    double beta = std::sqrt(1.0-1.0/((correctedEnergy/511.0+1.0)
        *(correctedEnergy/511.0+1.0)));
    double KforLoss = 0.734*std::pow(substrateCharge, 0.037);
    double JforLoss = (substrateCharge < 13) ? 0.0115*substrateCharge :
                        0.00976*substrateCharge+0.0585
                        *std::pow(substrateCharge, -0.19);

    JforLoss = JforLoss/(1.0+KforLoss*JforLoss/correctedEnergy);
    double ionizationEnergy = 153.55*1E24*substrateDensity*substrateCharge
        /(Constants::avogadrosNumber)/(beta*beta)*(std::log(511.0*beta*beta
        *correctedEnergy/(2.0*JforLoss*JforLoss*(1.0-beta*beta)))-std::log(2.0)
        *(2.0*std::sqrt(1.0-beta*beta)-1+beta*beta)+1-beta*beta+1.0/8.0
        *(1.0-std::sqrt(1.0-beta*beta))*(1.0-std::sqrt(1.0-beta*beta)));
    // Final expression ensures a positive value
    return (ionizationEnergy + std::abs(ionizationEnergy))/2.0;
};

inline double Electron::getBremsstrahlung() {
    double correctedEnergy = energy*correctionFactor;
    double bremsstralung = 1.4e-4*1E24*substrateDensity/(Constants::avogadrosNumber)
        *substrateCharge*(substrateCharge+1.0)*(correctedEnergy+511.0)*(4.0
        *std::log(2.0*(correctedEnergy+511.0)/511.0)-4.0/3.0);
    // Final expressio ensures a positive value
    return (bremsstralung + std::abs(bremsstralung))/2.0;
};

inline Coordinate Electron::calculateNewCoordinate(size_t i) {
    return {
        coordinate.x + flyingDistances[i]*10e8*std::sin(velocity.zAngle)
                        *std::cos(velocity.xAngle),
        coordinate.y + flyingDistances[i]*10e8*std::sin(velocity.zAngle)
                        *std::sin(velocity.xAngle),
        coordinate.z + flyingDistances[i]*10e8*std::cos(velocity.zAngle),
            coordinate.depth
    };
};

// InputFields functions
InputFields::InputFields() 
    : charge((Defaults::type == ION) ? Defaults::ionCharge : Constants::electronCharge),
    electronStoppingEnergy(Defaults::electronStoppingEnergy),
    energy((Defaults::type == ION) ? Defaults::ionEnergy : Defaults::electronEnergy),
    enableDamageCascade(Defaults::enableDamageCascade),
    electronScreeningParametersFilename(Defaults::electronScreeningParametersFilename),
    inputDirectoryName(Defaults::inputDirectory),
    ionDisplacementEnergy(Defaults::ionDisplacementEnergy),
    ionStoppingEnergy(Defaults::ionStoppingEnergy),
    logSingleDisplacement(Defaults::logSingleDisplacement),
    mass((Defaults::type == ION) ? Defaults::ionMass : Constants::electronMass),
    mottScatteringParametersFilename(Defaults::mottScatteringParametersFilename),
    numAngleDivisors(Defaults::numAngleDivisors),
    numElecScreeningPotentialElements(Defaults::numElecScreeningPotentialElements),
    numFlyingDistances(Defaults::numFlyingDistances),
    numMottScatteringPotentialElements(Defaults::numMottElements),
    outputFileEndMarker(Defaults::outputFileEnd),
    outputFileExtension(Defaults::outputFileExtension),
    outputCoordinateFilename(Defaults::outputCoordinateFilename),
    outputDirectory(Defaults::outputDirectory),
    range(Defaults::range),
    simulationCount(Defaults::simulationCount),
    substrateCharge(Defaults::substrateCharge),
    substrateDensity(Defaults::substrateDensity),
    substrateMass(Defaults::substrateMass),
    settingsFilename(Defaults::settingsFilename),
    type(Defaults::type),
    logEndOfFlyingDistanceOnly(Defaults::logEndOfFlyingDistanceOnly),
    logStoppingPointOnly(Defaults::logStoppingPointOnly),
    progressChecking(Defaults::progressChecking) {}


InputFields::InputFields(
    double charge,
    double electronStoppingEnergy,
    double energy,
    bool enableDamageCascade,
    const std::string& electronScreeningParametersFilename,
    const std::string& inputDirectoryName,
    double ionDisplacementEnergy,
    double ionStoppingEnergy,
    bool logSingleDisplacement, 
    double mass,
    const std::string& mottScatteringParametersFilename,
    size_t numAngleDivisors,
    size_t numElecScreeningPotentialElements,
    size_t numFlyingDistances,
    size_t numMottScatteringPotentialElements,
    const std::string& outputFileEndMarker,
    const std::string& outputFileExtension,
    const std::string& outputCoordinateFilename,
    const std::string& outputDirectory,
    double range,
    size_t simulationCount,
    double substrateCharge,
    double substrateDensity,
    double substrateMass,
    const std::string& settingsFilename,
    ParticleType type,
    bool logEndOfFlyingDistanceOnly,
    bool logStoppingPointOnly,
    bool progressChecking,
    size_t numThreads)
    : 
    charge(charge),
    electronStoppingEnergy(electronStoppingEnergy),
    energy(energy),
    enableDamageCascade(enableDamageCascade),
    electronScreeningParametersFilename(electronScreeningParametersFilename),
    inputDirectoryName(inputDirectoryName),
    ionDisplacementEnergy(ionDisplacementEnergy),
    ionStoppingEnergy(ionStoppingEnergy),
    logSingleDisplacement(logSingleDisplacement),
    mass(mass),
    mottScatteringParametersFilename(mottScatteringParametersFilename),
    numAngleDivisors(numAngleDivisors),
    numElecScreeningPotentialElements(numElecScreeningPotentialElements),
    numFlyingDistances(numFlyingDistances),
    numMottScatteringPotentialElements(numMottScatteringPotentialElements),
    outputFileEndMarker(outputFileEndMarker),
    outputFileExtension(outputFileExtension),
    outputCoordinateFilename(outputCoordinateFilename),
    outputDirectory(outputDirectory),
    range(range),
    simulationCount(simulationCount),
    substrateCharge(substrateCharge),
    substrateDensity(substrateDensity),
    substrateMass(substrateMass),
    settingsFilename(settingsFilename),
    type(type),
    logEndOfFlyingDistanceOnly(logEndOfFlyingDistanceOnly),
    logStoppingPointOnly(logStoppingPointOnly),
    progressChecking(progressChecking),
    numThreads(numThreads) {};

InputFields::InputFields(const std::string& settingsFilename)
        : InputFields() {
    this->settingsFilename = settingsFilename;
};

std::shared_ptr<InputFields> InputFields::instance;

std::shared_ptr<InputFields> InputFields::getInstance() {
    // Create the instance if it doesn't exist
    if (!instance) {
        instance = std::shared_ptr<InputFields>(new InputFields());
    }
    return instance;
}

std::shared_ptr<InputFields> InputFields::getInstance(
        const std::string& settingsFilename) {
    // Create the instance if it doesn't exist
    if (!instance) {
        instance = std::shared_ptr<InputFields>(new InputFields(settingsFilename));
    }
    return instance;
}

// Getters
bool InputFields::getEnableDamageCascade() const { return enableDamageCascade; }
double InputFields::getElectronStoppingEnergy() const { return electronStoppingEnergy; };
bool InputFields::getLogSingleDisplacement() const { return logSingleDisplacement; }
double InputFields::getCharge() const { return charge; }
double InputFields::getEnergy() const { return energy; }
double InputFields::getIonDisplacementEnergy() const { return ionDisplacementEnergy; };
double InputFields::getIonStoppingEnergy() const { return ionStoppingEnergy; };
double InputFields::getMass() const { return mass; }
double InputFields::getRange() const { return range; }
double InputFields::getSubstrateCharge() const { return substrateCharge; }
double InputFields::getSubstrateDensity() const { return substrateDensity; }
double InputFields::getSubstrateMass() const { return substrateMass; }
ParticleType InputFields::getType() const { return type; }
size_t InputFields::getSimulationCount() const { return simulationCount; }
size_t InputFields::getNumAngleDivisors() const { return numAngleDivisors; }
size_t InputFields::getNumElecScreeningPotentialElements() const { return numElecScreeningPotentialElements; }
size_t InputFields::getNumFlyingDistances() const { return numFlyingDistances; }
size_t InputFields::getNumMottScatteringPotentialElements() const { return numMottScatteringPotentialElements; }
const std::string& InputFields::getElectronScreeningParametersFilename() const { return electronScreeningParametersFilename; }
const std::string& InputFields::getInputDirectoryName() const { return inputDirectoryName; }
const std::string& InputFields::getMottScatteringParametersFilename() const { return mottScatteringParametersFilename; }
const std::string& InputFields::getOutputCoordinateFilename() const { return outputCoordinateFilename; }
const std::string& InputFields::getOutputDirectory() const { return outputDirectory; }
const std::string& InputFields::getOutputFileEndMarker() const { return outputFileEndMarker; }
const std::string& InputFields::getOutputFileExtension() const { return outputFileExtension; }
const std::string& InputFields::getSettingsFilename() const { return settingsFilename; }
bool InputFields::getLogEndOfFlyingDistanceOnly() const { return logEndOfFlyingDistanceOnly;};
bool InputFields::getLogStoppingPointOnly() const { return logStoppingPointOnly;};
bool InputFields::getProgressChecking() const { return progressChecking; };
size_t InputFields::getNumThreads() const { return numThreads; };

// Setters
void InputFields::setCharge(double charge) { this->charge = charge; }
void InputFields::setElectronStoppingEnergy(double energy) { this->electronStoppingEnergy = energy; };
void InputFields::setElectronScreeningParametersFilename(const std::string& filename) { this->electronScreeningParametersFilename = filename; }
void InputFields::setEnableDamageCascade(bool enable) { this->enableDamageCascade = enable; }
void InputFields::setEnergy(double energy) { this->energy = energy; }
void InputFields::setIonDisplacementEnergy(double energy) { this->ionDisplacementEnergy = energy; };
void InputFields::setIonStoppingEnergy(double energy) { this->ionStoppingEnergy = energy; };
void InputFields::setInputDirectoryName(const std::string& name) { this->inputDirectoryName = name; }
void InputFields::setLogSingleDisplacement(bool log) { this->logSingleDisplacement = log; }
void InputFields::setMass(double mass) { this->mass = mass; }
void InputFields::setMottScatteringParametersFilename(const std::string& filename) { this->mottScatteringParametersFilename = filename; }
void InputFields::setNumAngleDivisors(size_t divisors) { this->numAngleDivisors = divisors; }
void InputFields::setNumElecScreeningPotentialElements(size_t elements) { this->numElecScreeningPotentialElements = elements; }
void InputFields::setNumFlyingDistances(size_t distances) { this->numFlyingDistances = distances; }
void InputFields::setNumMottScatteringPotentialElements(size_t elements) { this->numMottScatteringPotentialElements = elements; }
void InputFields::setOutputCoordinateFilename(const std::string& filename) { this->outputCoordinateFilename = filename; }
void InputFields::setOutputDirectory(const std::string& directory) { this->outputDirectory = directory; }
void InputFields::setOutputFileEndMarker(const std::string& marker) { this->outputFileEndMarker = marker; }
void InputFields::setOutputFileExtension(const std::string& extension) { this->outputFileExtension = extension; }
void InputFields::setRange(size_t range) { this->range = range; }
void InputFields::setSettingsFilename(const std::string& filename) { this->settingsFilename = filename; }
void InputFields::setSimulationCount(size_t count) { this->simulationCount = count; }
void InputFields::setSubstrateCharge(double substrateCharge) { this->substrateCharge = substrateCharge; }
void InputFields::setSubstrateDensity(double substrateDensity) { this->substrateDensity = substrateDensity; }
void InputFields::setSubstrateMass(double substrateMass) { this->substrateMass = substrateMass; }
void InputFields::setType(ParticleType type) { this->type = type; }
void InputFields::setLogEndOfFlyingDistanceOnly(bool logEndOfFlyingDistanceOnly) {this->logEndOfFlyingDistanceOnly = logEndOfFlyingDistanceOnly; }
void InputFields::setLogStoppingPointOnly(bool logStoppingPointOnly) {this->logStoppingPointOnly = logStoppingPointOnly; }
void InputFields::setProgressChecking(bool progressChecking) {this->progressChecking = progressChecking; }
void InputFields::setNumThreads(size_t numThreads) {this->numThreads = numThreads; }

bool InputFields::parseSetting(const std::string& line, std::string& name, std::string& value) {
        std::istringstream iss(line);
    if (std::getline(iss, name, '=') && std::getline(iss, value)) {
        return true;
    }
    return false;
} 

bool InputFields::stringToBool(const std::string& str) {
    return str == "True" || str == "true" || str == "1";
}

void InputFields::readSettingsFromFile() {
    std::shared_ptr<InputFields> input;
    std::ifstream file(settingsFilename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file \""
            <<  settingsFilename
            << "\". Using default settings. Use -s for settings help."
            << std::endl;
        if(promptContinue()) {
            input = InputFields::getInstance();
        } else {
            std::exit(EXIT_FAILURE);
        }
    } else {
        std::string line;
        input = InputFields::getInstance(settingsFilename);
        while (std::getline(file, line)) {
            std::string name, value;
            if (parseSetting(line, name, value)) {
                if (name == "electronMode") {
                    input->setType(stringToBool(value) ? ELECTRON : ION);
                } else if (name == "electronEnergy(keV)") {
                    input->setEnergy(std::stod(value));
                } else if (name == "electronStoppingEnergy(keV)") {
                    input->setElectronStoppingEnergy(std::stod(value));
                } else if (name == "electronScreeningParametersFilename") {
                    value.erase(std::remove(value.begin(), value.end(), '\"'), value.end());
                    input->setElectronScreeningParametersFilename(value);
                } else if (name == "numElecScreeningPotentialElements") {
                    input->setNumElecScreeningPotentialElements(std::stoi(value));
                } else if (name == "mottScatteringParametersFilename") {
                    value.erase(std::remove(value.begin(), value.end(), '\"'), value.end());
                    input->setMottScatteringParametersFilename(value);
                } else if (name == "numMottScatteringPotentialElements") {
                    input->setNumMottScatteringPotentialElements(std::stoi(value));
                } else if (name == "numAngleDivisors") {
                    input->setNumAngleDivisors(std::stoi(value));
                } else if (name == "numFlyingDistances") {
                    input->setNumFlyingDistances(std::stoi(value));
                } else if (name == "enableDamageCascade") {
                    input->setEnableDamageCascade(stringToBool(value));
                } else if (name == "ionCharge(e)") {
                    input->setCharge(input->getType() == ELECTRON ? Constants::electronCharge : std::stod(value));
                } else if (name == "ionEnergy(keV)") {
                    input->setEnergy(input->getType() == ELECTRON ? Defaults::electronEnergy : std::stod(value));
                } else if (name == "ionMass(amu)") {
                    input->setMass(input->getType() == ELECTRON ? Constants::electronMass : std::stod(value));
                } else if (name == "substrateDisplacementEnergy(keV)") {
                    input->setIonDisplacementEnergy(std::stod(value));
                } else if (name == "ionStoppingEnergy(keV)") {
                    input->setIonStoppingEnergy(std::stod(value));
                } else if (name == "logSingleDisplacement") {
                    input->setLogSingleDisplacement(stringToBool(value));
                } else if (name == "inputDirectoryName") {
                    value.erase(std::remove(value.begin(), value.end(), '\"'), value.end());
                    input->setInputDirectoryName(value);
                } else if (name == "outputCoordinateFilename") {
                    value.erase(std::remove(value.begin(), value.end(), '\"'), value.end());
                    input->setOutputCoordinateFilename(value);
                } else if (name == "outputFileEndMarker") {
                    value.erase(std::remove(value.begin(), value.end(), '\"'), value.end());
                    input->setOutputFileEndMarker(value);
                } else if (name == "outputFileExtension") {
                    value.erase(std::remove(value.begin(), value.end(), '\"'), value.end());
                    input->setOutputFileExtension(value);
                } else if (name == "outputDirectory") {
                    value.erase(std::remove(value.begin(), value.end(), '\"'), value.end());
                    input->setOutputDirectory(value);
                } else if (name == "substrateCharge(e)") {
                    input->setSubstrateCharge(std::stod(value));
                } else if (name == "substrateMass(amu)") {
                    input->setSubstrateMass(std::stod(value));
                } else if (name == "settingsFilename") {
                    input->setSettingsFilename(value);
                } else if (name == "logEndOfFlyingDistanceOnly") {
                    input->setLogEndOfFlyingDistanceOnly(stringToBool(value));
                } else if (name == "logStoppingPointOnly") {
                    input->setLogStoppingPointOnly(stringToBool(value));
                } else if (name == "simulationCount") {
                    input->setSimulationCount(std::stoi(value));
                } else if (name == "numThreads") {
                    input->setNumThreads(std::stoi(value));
                } else
                    // Handle other settings here if needed
                    std::cerr << "Warning: Unknown setting name '" << name << "'" << std::endl;
                }
            }
        }
    file.close();
}

std::string InputFields::printInputFields() const {
    std::ostringstream oss; // Create a std::ostringstream to store the output

    // Redirect output to std::ostringstream
    oss << "Type:\t" << (type == ELECTRON ? "Electron" : "Ion") << std::endl;
    oss << "Energy(keV):\t" << energy << std::endl;
    oss << "Initial Particle Charge(e):\t" << charge << std::endl;
    oss << "Initial Particle Mass(amu for ion, kg for electron):\t" << mass << std::endl;
    oss << "Electron Stopping Energy(keV)\t" << charge << std::endl;
    oss << "Electron Screening Parameters Filename:\t" << electronScreeningParametersFilename << std::endl;
    oss << "Number of Electron Screening Potential Elements:\t" << numElecScreeningPotentialElements << std::endl;
    oss << "Mott Scattering Parameters Filename:\t" << mottScatteringParametersFilename << std::endl;
    oss << "Number of Mott Scattering Potential Elements:\t" << numMottScatteringPotentialElements << std::endl;
    oss << "Number of Angle Divisors:\t" << numAngleDivisors << std::endl;
    oss << "Number of Flying Distances:\t" << numFlyingDistances << std::endl;
    oss << "Enable Damage Cascade:\t" << (enableDamageCascade ? "true" : "false") << std::endl;
    oss << "Substrate Displacement Energy (keV):\t" << ionDisplacementEnergy << std::endl;
    oss << "Ion Stopping Energy (keV):\t" << ionStoppingEnergy << std::endl;
    oss << "Log Single Displacement:\t" << (logSingleDisplacement ? "true" : "false") << std::endl;
    oss << "Input Directory Name:\t" << inputDirectoryName << std::endl;
    oss << "Output Coordinate Filename:\t" << outputCoordinateFilename << std::endl;
    oss << "Output File Extension:\t" << outputFileExtension << std::endl;
    oss << "Output File End Marker:\t" << outputFileEndMarker << std::endl;
    oss << "Output Directory:\t" << outputDirectory << std::endl;
    oss << "Substrate Charge:\t" << substrateCharge << std::endl;
    oss << "Substrate Density:\t" << substrateDensity << std::endl;
    oss << "Substrate Mass:\t" << substrateMass << std::endl;
    oss << "Settings Filename:\t" << settingsFilename << std::endl;
    // oss << "Range(unused):\t" << range << std::endl;
    oss << "Log Stopping Point Only:\t" << (logStoppingPointOnly ? "true" : "false") << std::endl;
    oss << "Log End of Flying Distance Only:\t" << (logEndOfFlyingDistanceOnly ? "true" : "false") << std::endl;
    oss << "Number of Threads:\t" << numThreads << std::flush;
    oss << "Simulation Count:\t" << simulationCount << std::flush;

    // Return the contents of the std::ostringstream as a string
    return oss.str();
}

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

void Bombardment::initiate(std::shared_ptr<InputFields> input) {
    DEBUG_PRINT("Initiating bombardment " << id);
    std::unique_ptr<Particle> initialParticle;
    Coordinate initialCoordinate = {0.0, 0.0, 0.0, 0};
    Velocity initialVelocity = {0.0, 0.0, input->getEnergy()};

    if(input->getType() == ION) {
        initialParticle = std::make_unique<Ion>(initialCoordinate, initialVelocity,
            input, shared_from_this());
    } else {
        initialParticle = std::make_unique<Electron>(initialCoordinate, initialVelocity, 
            input, shared_from_this());
    }

    addParticle(std::move(initialParticle));

    if(auto locked = simulation.lock()) {
        #ifndef NO_OUTPUT
        locked->writeData(COORDINATE, getParticles(), id);
        #endif
    } else {
        std::filesystem::path cwd = std::filesystem::current_path();
        std::cerr << "Failed to initiate write to " << 
            (cwd / input->getOutputDirectory()
            / (input->getOutputCoordinateFilename()
            + input->getOutputFileExtension())) << std::endl;
    }
};

void Bombardment::addParticle(std::unique_ptr<Particle> particle) {
    particles.emplace_back(std::move(particle));
    particles.back()->fire();
}

Simulation::Simulation(std::shared_ptr<InputFields> input)
    : input(input),
    bombardments(input->getSimulationCount()),
    outputPath(fs::current_path() / input->getOutputDirectory()),
    progressCounter(0) {};

void Bombardment::popParticle() {
    particles.pop_back();
};

Simulation::~Simulation() {};

std::shared_ptr<InputFields> Simulation::getInput() { return input; };
size_t Simulation::getBombardmentSize() { return bombardments.capacity(); };
size_t Simulation::getBombardmentCapacity() { return bombardments.capacity(); };

void Simulation::initiate() {

    bombardmentCount = 0;

    #ifndef NO_OUTPUT
    checkOutputFiles();
    #endif

    for(size_t i = 0; i < input->getSimulationCount(); i++) {
        bombardments[i] = std::make_shared<Bombardment>(shared_from_this(), i);
    }

    const size_t simulation_count = input->getSimulationCount();
    threads.reserve(simulation_count);

    DEBUG_PRINT("Initiating simulation");
    size_t numThreadsPerBatch = input->getNumThreads();
    size_t numBatches = (simulation_count+numThreadsPerBatch-1)
        /numThreadsPerBatch;

    for (size_t batchIndex = 0; batchIndex < numBatches; ++batchIndex) {
        size_t startIndex = batchIndex * numThreadsPerBatch;
        size_t endIndex = std::min(startIndex + numThreadsPerBatch,
            simulation_count);

        std::vector<std::thread> batchThreads;

        for (size_t i = startIndex; i < endIndex; ++i) {
            batchThreads.emplace_back([this, i]() {
                bombardments[i]->initiate(input);
                bombardments[i].reset();
                if (input->getProgressChecking()) {
                    updateProgressCounter();
                }
            });
        }

        for (auto& thread : batchThreads) {
            if (thread.joinable()) {
                thread.join();
            }
        }
    }

    #ifndef NO_OUTPUT
    writeData(ENDOFFILE, {}, 0);
    #endif
}

void Simulation::writeData(
    OutputType outputType,
    const std::vector<std::unique_ptr<Particle>>& particles = {},
    size_t bombardmentID = 0) {
    
    // Acquire lock
    std::lock_guard<std::mutex> lock(fileLock);

    DEBUG_PRINT("Simulation " << bombardmentID << " Writing");

    // Check if output directory exists, create it if not
    if (!fs::exists(outputPath)) {
        const std::string ANSI_COLOR_RED = "\033[1;31m";
        const std::string ANSI_COLOR_RESET = "\033[0m";
        std::cerr << "Unable to access output directory " << outputPath.string()
            << ".\nCreating directory "
            << outputPath.string()
            << ".\n"
            << ANSI_COLOR_RED << "Warning: Cancellation at this point will lead to loss of simulation data." << ANSI_COLOR_RESET
            << std::endl;
        if(promptContinue()) {
            fs::create_directories(outputPath);
        } else {
            exit(EXIT_FAILURE);
        }
    }

    // Open output file
    std::string filename = std::string(input->getOutputCoordinateFilename());
    filename = filename + std::string(input->getOutputFileExtension());
    std::ofstream outputFile(outputPath / filename, std::ios_base::app);
    if (!outputFile.is_open() && outputType != ENDOFFILE) {
        std::cerr << "Failed to open output file "
                  << (outputPath / filename).string() 
                  << ".\nData will be lost for simulation "
                  << bombardmentID
                  << "."
                  << std::endl;
        return;
    } else if (!outputFile.is_open() && outputType == ENDOFFILE){
        std::cerr << "Failed to open output file "
                  << (outputPath / filename).string()
                  << ".\nFile end marker \""
                  << input->getOutputFileEndMarker() << "\""
                  << " could not be appended." 
                  << std::endl;
    }

    // Write headers if the file is empty
    if (outputFile.tellp() == 0) {
        outputFile << "bombardmentID,particleID,depth,x,y,z\n";
    }

    size_t particleID = 0;

    // Write data to file
    switch (outputType) {
        case OutputType::COORDINATE:
            for (const auto& particle : particles){
                if(particle->getCoordinates().empty()) {
                    if(input->getLogSingleDisplacement())
                    {
                        outputFile << bombardmentID << ','
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
                        outputFile << bombardmentID << ','
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
            outputFile << input->getOutputFileEndMarker();
            break; 
        default:
            std::cerr << "Unsupported output type." << std::endl;
            break;
    }

    // Close output file
    outputFile.close();
};

bool Simulation::fileIsWritten(const fs::path filename) {
    
    const std::streampos bufferSize =
        static_cast<std::streampos>(input->getOutputFileEndMarker().length());

    std::ifstream file(filename, std::ios::ate);
    if (!file.is_open()) {
        return false;
    }

    if (file.tellg() < bufferSize) {
        return false;
    }
    
    // Move file pointer to bufferSize characters before the end
    file.seekg(-bufferSize, std::ios::end);
    // Read the last bufferSize characters
    std::string lastChars(bufferSize, '\0');
    file.read(&lastChars[0], bufferSize);
    return lastChars == std::string(input->getOutputFileEndMarker());
};

void Simulation::renameFileWithTimestamp(const std::string& filename) {
    std::string timestamp = getCurrentDateTime();
    std::string newFilename = filename + "_" + timestamp
        + input->getOutputFileExtension();
    fs::path newFilePath = outputPath/newFilename;
    std::string oldFilename = filename + input->getOutputFileExtension();
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
}

void Simulation::checkOutputFiles() {
    fs::path filepath;
    std::string filename;
    for(int i = COORDINATE; i  < ENDOFFILE; i++) {
        switch(i) {
            case COORDINATE:
                filename = std::string(input->getOutputCoordinateFilename());
                filename = filename
                    + std::string(input->getOutputFileExtension());
                filepath = outputPath/filename;
                if(fileIsWritten(filepath)) {
                    filename = std::string(
                        input->getOutputCoordinateFilename());
                    renameFileWithTimestamp(filename);
                }
                break;
            default:
                break;
        }
    }
}

void Simulation::updateProgressCounter() {
    std::lock_guard<std::mutex> lock(counterLock);
    ++progressCounter;
        constexpr int progressBarWidth = 80;
        int progress = static_cast<int>((static_cast<double>(progressCounter)
            /static_cast<double>(input->getSimulationCount()))
            *progressBarWidth);

        std::cout << "Progress: [";
        for (int i = 0; i < progressBarWidth; ++i) {
            if (i < progress) {
                std::cout << "=";
            } else {
                std::cout << " ";
            }
        }
        std::cout << "] "
                    << std::setw(3)
                    << std::setfill(' ')
                    << std::right
                    << ((progressCounter * 100) / input->getSimulationCount())
                    << "%\r";
        std::cout.flush();

}

#endif