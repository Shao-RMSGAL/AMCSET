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

// Particle functions and static member initializers
std::mt19937 Particle::randomGenerator;
std::uniform_real_distribution<double> Particle::randomDistribution(
    0+std::numeric_limits<double>::epsilon(),
    1-std::numeric_limits<double>::epsilon());

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

inline double Particle::atomicSpacing() const {
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

inline Coordinate Particle::calculateNewCoordinate() {
    throw std::logic_error("Base class Particle::calculateNewCoordinate() called");
    Coordinate dummy;
    return dummy;
}

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

void Particle::createSubstrateKnockon(Coordinate& newCoordinate, Velocity& targetVelocity) {
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
};

void Ion::fire() {
    Velocity newVelocity;
    Velocity targetVelocity;
    Coordinate newCoordinate;
    if(velocity.energy <= Constants::ionStoppingEnergy) {
        // DEBUG_PRINT("Particle " << id  << " does not have enough energy");
        // addCoordinate(coordinate);
        return;
    } else {
        // DEBUG_PRINT("Firing particle " << coordinate.id);
        while(velocity.energy > Constants::ionStoppingEnergy) {
            // DEBUG_PRINT("Energy " << fireCount << ":\t" << velocity.energy);
            velocity.energy -= electronicStoppingEnergy();

            std::tie(newVelocity, targetVelocity) = recoilEnergyAndVelocity();

            newCoordinate = calculateNewCoordinate();

            if(newCoordinate.z < 0.0) {
                DEBUG_PRINT("Ion Sputter");
                return;
            }

            if(Defaults::enableDamageCascade) {
                createSubstrateKnockon(newCoordinate, targetVelocity);
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
Substrate::Substrate(Coordinate coordinate, Velocity velocity, double charge,
    double mass, double density, double range,
    std::weak_ptr<Bombardment> bombardment) :
    Ion(coordinate, velocity, charge, mass, mass, density, mass, SUBSTRATE, range,
        bombardment) {};

const size_t& Substrate::getDepth() const { return depth; };

std::array<std::array<std::array<double,
            Constants::numMottkParam>,Constants::numMottjParam>,
            Constants::numMottElements> Electron::mottScatteringParams;

std::array<std::array<double,
            Constants::numElecScreeningPotentialParams>,
            Defaults::numElecScreeningPotentialElements> Electron::elecScreeningParams;

std::vector<double> Electron::divisionAngles;

//  Electron functions
Electron::Electron(Coordinate initialCoordinate, Velocity initialVelocity, InputFields input,
        std::weak_ptr<Bombardment> bombardment)
        : Particle(initialCoordinate, initialVelocity,
        Constants::electronCharge, Constants::electronMass,
        input.getSubstrateCharge(), input.getSubstrateDensity(),
        input.getSubstrateMass(), ELECTRON, input.getRange(),
        bombardment), randomCrossSections(Defaults::numFlyingDistances),
        flyingDistances(Defaults::numFlyingDistances),
        scatteringAngles(Defaults::numFlyingDistances),
        partialCrossSections(Defaults::numAngleDivisors),
        elasticEnergyLoss(Defaults::numFlyingDistances),
        correctionFactor(1.0) {;
    DEBUG_PRINT("Electron made");
};

inline double Electron::findAngle(size_t i) {
    // Value is offset from 1 to avoid checking a zero-angle
    return -M_PI/(9.0)*(1-std::pow(10.0,(static_cast<double>(i+1)
        /Defaults::numAngleDivisors)));
};

void Electron::readParametersAndInitialize() {
    DEBUG_PRINT("Reading Mott scattering values");
    fs::path mottParametersPath(Defaults::inputDirectory);
    mottParametersPath = mottParametersPath / Defaults::mottScatteringParametersFilename ;

    std::string input;
    std::string parameter;
    std::istringstream stringStream;

    std::ifstream mottParametersFile(mottParametersPath);

    if(!mottParametersFile.is_open()) {
        throw std::ios_base::failure("Could not open file " + mottParametersPath.string() );
    }

    for(size_t i = 0; i < Constants::numMottElements; i++) {
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

    DEBUG_PRINT("Electron screening parameters");
    fs::path elecParametersPath(Defaults::inputDirectory);
    elecParametersPath = elecParametersPath / Defaults::electronScreeningParametersFilename;

    std::ifstream elecParametersFile(elecParametersPath);

    if(!elecParametersFile.is_open()) {
        throw std::ios_base::failure("Could not open file " + elecParametersPath.string() );
    }

    for(size_t i = 0; i < Defaults::numElecScreeningPotentialElements; i++) {
        std::getline(elecParametersFile, input);
        stringStream.clear();
        stringStream.str(input);
        for(size_t j = 0; j < Constants::numElecScreeningPotentialParams; j++) {
            std::getline(stringStream, parameter, ',');
            elecScreeningParams[i][j] = std::stod(parameter);
        }
    }

    elecParametersFile.close();

    DEBUG_PRINT("Calculating division angles");

    divisionAngles.reserve(Defaults::numAngleDivisors);

    for(size_t i = 0; i < Defaults::numAngleDivisors; i++) {
        divisionAngles.push_back(findAngle(i));
    }
    
};

void Electron::fire() {
    DEBUG_PRINT("Firing electron!");
    // size_t numParticles = getNumParticles();
    double totalMottCrossSection;
    double totalElasticEnergyLoss;
    double ionizationEnergyLoss;
    double bremsstrahlung;
    double newEnergy;
    bool firstIteration = true;
    double flightGroupLength;
    Velocity newVelocity;
    Velocity targetVelocity;
    Coordinate newCoordinate;
    while(velocity.energy > Constants::electronStoppingEnergy) {
        DEBUG_PRINT("Energy: " << energy << " keV");

        totalMottCrossSection = getMottTotalCrossSection();
        
        DEBUG_PRINT("Total Mott: " << totalMottCrossSection);

        flightGroupLength = generateFlyingDistances(totalMottCrossSection);

        generateRandomCrossSections();
        generateIntegratedCrossSections(totalMottCrossSection);

        

        // for(auto& i : scatteringAngles) DEBUG_PRINT("Random angle: " << i*180.0/M_PI << " degrees");

        DEBUG_PRINT("Factor: " << correctionFactor);

        totalElasticEnergyLoss = getElasticEnergyLoss();

        DEBUG_PRINT("Elastic: " << totalElasticEnergyLoss);

        DEBUG_PRINT("Length: " << flightGroupLength);

        ionizationEnergyLoss = flightGroupLength*getIonizationEnergyLoss();

        DEBUG_PRINT("Ionization: " << ionizationEnergyLoss);

        bremsstrahlung = flightGroupLength*getBremsstrahlung();

        DEBUG_PRINT("Bremsstralung: " << bremsstrahlung);

        newEnergy = energy - totalElasticEnergyLoss - ionizationEnergyLoss -
            bremsstrahlung;

        // TODO: Include logic for explicit and implicit methods
        if(firstIteration) {
            DEBUG_PRINT("Preassign: " << correctionFactor);
            DEBUG_PRINT("Energy: " << energy);
            DEBUG_PRINT("New Energy: " << newEnergy);
            correctionFactor = (energy + newEnergy)/2.0/energy;
            DEBUG_PRINT("Postcorrection: " << correctionFactor);
            firstIteration = false;
            continue;
        } 
        correctionFactor = (energy + newEnergy)/2.0/energy;
        energy = newEnergy;

        for(size_t i = 0; i < Defaults::numFlyingDistances; i++) {

            std::tie(newVelocity, targetVelocity) = relativeToAbsoluteVelocity(
                scatteringAngles[i], (M_PI-scatteringAngles[i])/2.0, elasticEnergyLoss[i]
            );

            newCoordinate = calculateNewCoordinate(i);

            if(newCoordinate.z < 0.0) {
                DEBUG_PRINT("Electron Sputter");
                return;
            }

            if(Defaults::enableDamageCascade) {
                createSubstrateKnockon(newCoordinate, targetVelocity);
            }

            addCoordinate(coordinate);
            coordinate = newCoordinate;
            velocity = newVelocity;
        }

        
    }   
};

double Electron::getMottDifferentialCrossSection(double Theta) {
    double correctedEnergy = energy*correctionFactor;
    double beta = std::sqrt(1.0-1.0/((correctedEnergy/511.0 + 1.0)
        *(correctedEnergy/511.0 + 1.0)));
    // DEBUG_PRINT("Beta: " << correctionFactor);
    constexpr double betaAVE = 0.7181287;
    double gamma = correctedEnergy/511.0+1.0;
    constexpr double elecmass = 1; // specific for a.u. unit,  only for Feq calcu.
    constexpr double Vc = 137;  //speed of light for a.u. unit, only for Feq calcu.
    constexpr double re1 = 2.817938E-13;  // in the unit of cm)
    double Rmott = 0.0;
    double alpha1;
    for(size_t j = 0; j < Constants::numMottjParam; j++) {
        alpha1 = 0.0;
        for(size_t k = 0; k < Constants::numMottkParam; k++) {
            alpha1 = alpha1 +
                mottScatteringParams[static_cast<size_t>(substrateCharge)][j][k]
                *std::pow(beta-betaAVE, k);
            // DEBUG_PRINT("Expr: " << std::pow(beta-betaAVE, k-1.0));
        }
        Rmott = Rmott+alpha1*std::pow(1.0-std::cos(Theta), (j)/2.0);
    }
    double DSC = (substrateCharge*re1)*(substrateCharge*re1)*(1.0-beta*beta)/
        std::pow(beta,4.0)*1.0/((1 - std::cos(Theta))*(1 - std::cos(Theta)));
    double moment = 2.0*beta*gamma*elecmass*Vc*std::sin(Theta/2.0);
    double Feq = 0;
    for(size_t i = 0; i < Constants::numElecScreeningPotentialParams/2; i++) {
        Feq = Feq + elecScreeningParams[substrateCharge][i]*
            elecScreeningParams[substrateCharge][i + 3]*elecScreeningParams[substrateCharge][i + 3]
            /(elecScreeningParams[substrateCharge][i + 3]*elecScreeningParams[substrateCharge][i + 3]
            +moment*moment);
    }
    return Rmott*DSC*(1-Feq)*(1-Feq);
};

void Electron::generateRandomCrossSections() {
    for(size_t i = 0; i < Defaults::numFlyingDistances; i++) {
        randomCrossSections[i] = random(); 
    }
    std::sort(randomCrossSections.begin(), randomCrossSections.end());
};

double Electron::generateFlyingDistances(double totalMottCrossSection) {
    double distance;
    double flightGroupLength = 0;
    double crossSectionFlight = 0;
    double atomicSpacing = std::pow(substrateDensity*1e24, -1.0/3.0);

    // DEBUG_PRINT("Cross section: " << totalMottCrossSection);

    for(size_t i = 0; i < Defaults::numFlyingDistances; i++) {
        crossSectionFlight = -std::log(random())/(substrateDensity*1e24*totalMottCrossSection);
        distance = crossSectionFlight+atomicSpacing;
        flightGroupLength += distance;
        flyingDistances[i] = distance; 
    }
    return flightGroupLength;
};

double Electron::getMottTotalCrossSection() {
    double totalCrossSection = 0.0;
    double differentialCrossSection = getMottDifferentialCrossSection(divisionAngles[0]);
    totalCrossSection += differentialCrossSection*2.0*M_PI
        *std::sin(divisionAngles[0])*(divisionAngles[1]-divisionAngles[0])/2.0;
    partialCrossSections[0] = totalCrossSection;
    for(size_t i = 1; i < Defaults::numAngleDivisors - 1; i++) {
        // DEBUG_PRINT("Diff: " << differentialCrossSection << "\tDiff (calc): " << differentialCrossSection*2.0*M_PI*std::sin(divisionAngles[i]) *(divisionAngles[i+1] - divisionAngles[i-1])/2.0 << "\tTotal " << i << ": " << totalCrossSection);
        differentialCrossSection = getMottDifferentialCrossSection(divisionAngles[i]);
        totalCrossSection += differentialCrossSection*2.0*M_PI*std::sin(divisionAngles[i])
            *(divisionAngles[i+1] - divisionAngles[i-1])/2.0;
        partialCrossSections[i] = totalCrossSection;
        // DEBUG_PRINT("Diff (calc): " << differentialCrossSection*2.0*M_PI*std::sin(divisionAngles[i])
        //     *(divisionAngles[i+1] - divisionAngles[i-1])/2.0);
        
        // DEBUG_PRINT("Partial cross " << i << ": " << totalCrossSection);
    }
    totalCrossSection += differentialCrossSection*2.0*M_PI*
        std::sin(divisionAngles[Defaults::numAngleDivisors-1])
        *(M_PI-(divisionAngles[Defaults::numAngleDivisors-3]
        -divisionAngles[Defaults::numAngleDivisors-2])/2.0);
    partialCrossSections[Defaults::numAngleDivisors-1] = totalCrossSection;
    return totalCrossSection;
};

void Electron::generateIntegratedCrossSections(double totalMottCrossSection) {
    // DEBUG_PRINT("Beginning integrated cross section");
    double prevFractionalCrossSection;
    double fractionalCrossSection;
    // TODO: Replace Defaults with parameter
    // DEBUG_PRINT("partialCrossSection Size: " << partialCrossSections.size());
    for(size_t i = 0; i < Defaults::numFlyingDistances; i++) {
        prevFractionalCrossSection = 0;
        // DEBUG_PRINT("Random: " << i << "\tCross: " << randomCrossSections[i]);
        for(size_t j = 0; j < partialCrossSections.size(); j++) {
            fractionalCrossSection = partialCrossSections[j]/totalMottCrossSection;
            if(fractionalCrossSection > randomCrossSections[i]) {
                // DEBUG_PRINT("Match! " << j);
                if(j > 0 && j < partialCrossSections.size()-1) {
                    scatteringAngles[i] = (divisionAngles[j] + divisionAngles[j + 1])/2.0
                        -(fractionalCrossSection-randomCrossSections[i])
                        *(divisionAngles[j+1]-divisionAngles[j-1])/2.0
                        /(fractionalCrossSection - prevFractionalCrossSection);
                    // DEBUG_PRINT("Between: " << scatteringAngles[i]);
                } else if (j == 0) {
                    scatteringAngles[i] = (divisionAngles[j]+divisionAngles[j+1])/2.0
                        -(fractionalCrossSection-randomCrossSections[i])
                        *(divisionAngles[j] + divisionAngles[j+1])/2.0
                        /(fractionalCrossSection - prevFractionalCrossSection);
                    // DEBUG_PRINT("Start: " << scatteringAngles[i]);
                } else {
                    scatteringAngles[i] = M_PI
                        -(fractionalCrossSection-randomCrossSections[i])
                        *(M_PI - (divisionAngles[j-2]-divisionAngles[j-3]))
                        /(fractionalCrossSection - prevFractionalCrossSection);
                    // DEBUG_PRINT("End: " << scatteringAngles[i]);
                }
                prevFractionalCrossSection = fractionalCrossSection;
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
    for(size_t i = 0; i < Defaults::numFlyingDistances; i++) {
        // DEBUG_PRINT("Angle: " << scatteringAngles[i]);
        double ElecREup = ((energy+511.0)*(std::sin(scatteringAngles[i]))
            *(std::sin(scatteringAngles[i]))+substrateMass*931.0*1000.0*(1
            -std::cos(scatteringAngles[i])))*energy*(energy+2.0*511.0);
        // DEBUG_PRINT("Top: " << ElecREup);
        double ElecREdown = (energy+substrateMass*931.5*1000.0)
            *(energy+substrateMass*931.5*1000.0)-energy*(energy+2.0*511.0)
            *(std::cos(scatteringAngles[i]))*(std::cos(scatteringAngles[i]));
        // DEBUG_PRINT("Bottom: " << ElecREdown);  
        elasticEnergyLoss[i] = ElecREup/ElecREdown;
        totalElasticEnergyLoss += elasticEnergyLoss[i];
    }
    // DEBUG_PRINT("Total loss: " << totalElasticEnergyLoss );
    return totalElasticEnergyLoss;
};

double Electron::getIonizationEnergyLoss() {
    // std::sqrt(1 - 1/((ElecE / 511 + 1)*(ElecE / 511 + 1)))
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
    DEBUG_PRINT("Ionization energy: " << ionizationEnergy << " keV/cm");
    // Final expression ensures a positive value
    return (ionizationEnergy + ::abs(ionizationEnergy))/2.0;
};

inline double Electron::getBremsstrahlung() {
    double correctedEnergy = energy*correctionFactor;
    double bremsstralung = 1.4e-4*1E24*substrateDensity/(Constants::avogadrosNumber)
        *substrateCharge*(substrateCharge+1.0)*(correctedEnergy+511.0)*(4.0
        *std::log(2.0*(correctedEnergy+511.0)/511.0)-4.0/3.0);
    DEBUG_PRINT("Bremsstralung: " << bremsstralung << " keV/cm");
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
        bombardments[i] = std::make_shared<Bombardment>(shared_from_this(), i);
    }

    const size_t simulation_count = inputs.getSimulationCount();
    threads.reserve(simulation_count); // Reserve space to avoid unnecessary reallocations

    DEBUG_PRINT("Initiating simulation");
    for(size_t i = 0; i < simulation_count; i++) {
        // DEBUG_PRINT("Initializing bombardment...");
        threads.emplace_back([this, i]() {
            bombardments[i]->initiate(inputs);
            // std::cout << "Use before: " << bombardments.at(i).use_count() << std::endl;
            bombardments[i].reset();
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