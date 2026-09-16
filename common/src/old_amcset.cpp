/******************************************************************************
 * Filename:        AMCSET.h
 * Project:         AMCSET
 * Author:          Nathaniel Thomas
 * Date Created:    March 15, 2024
 * Date Modified:   April 30, 2024
 * File Version:    1.1
 * Group:           Dr. Shao's RMSLCF Group
 *
 * Description:
 * This is the main header file for the AMCSET code developed by Dr. Shao and converted
 * to the C++ programming language by Nathaniel Thomaas
 *
 * Revision History
 * 1.0:
 * - File created and initial program structure made
 *****************************************************************************/
#ifndef AMCSET_H
#define AMCSET_H

#include <iostream>

// Main program  start
int startAMCSET(
    const int argc,
    const char* argv[],
    std::istream& inputStream = std::cin,
    std::ostream& outputStream = std::cout,
    std::ostream& errorStream = std::cerr
        );

#endif

/******************************************************************************
 * Filename:        AMCSET_classes.h
 * Project:         AMCSET
 * Author:          Nathaniel Thomas
 * Date Created:    March 8, 2024
 * Date Modified:   April 30, 2024
 * File Version:    1.0
 * Group:           Dr. Shao's RMSLCF Group
 *
 * Description:
 * This file contains the class declarations for the AMCSET program. 
 * implementations are in AMCSET_classes.cpp
 * 
 *
 * Revision History
 * 1.0:
 * - File created and class declarations added
 *
 *****************************************************************************/
#ifndef AMCSET_CLASSES_H
#define AMCSET_CLASSES_H

// Includes
#include <array>
#include <filesystem>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <random>
#include <thread>
#include <tuple>
#include <variant>
#include <vector>

// Local header files
#include "constants.h"


// TODO: Restore to only output on debug mode
// Debug compiler directives
#ifdef DEBUG_MODE
#define DEBUG_PRINT(x) std::cout << x << std::endl;
#else
#define DEBUG_PRINT(x)
#endif

// Namespaces
namespace fs = std::filesystem;

// TODO: Make this a dynamic type that can take three coordinates (type double),
// but can have an arbitrary number of additional fields of any data type. For example, 
// in addition to the size_t depth variable, it could have another entry called "energy"
// of type double. The size of the data structure should be determined at runtime depending on
// whether the user wants to include additional fields. The commented out code below is a
// good option. Consider enumerations instead of strings for the field names.
// // Define the types that can be stored in the additional fields
// using FieldValue = std::variant<size_t, double, std::string>;

// // Define a field as a pair of a string (the field name) and a FieldValue
// using Field = std::pair<std::string, FieldValue>;

// // Define a point as a tuple of three doubles and a vector of fields
// using Point = std::tuple<double, double, double, std::vector<Field>>;

struct Coordinate {
    double x;
    double y;
    double z;
    size_t depth;
};

struct Velocity{
    double zAngle;
    double xAngle;
    double energy;
};

class IOHandler;

class InputFields {
private:
    // Private constructors to prevent external instantiation 
    InputFields();
    InputFields(
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
        ParticleType type,
        bool logEndOfFlyingDistanceOnly,
        bool logStoppingPointOnly,
        bool progressChecking,
        size_t numThreads,
        std::shared_ptr<IOHandler> ioHandler);

    InputFields(std::shared_ptr<IOHandler> ioHandler);

    // Static instance of the class
    static std::shared_ptr<InputFields> instance;

    // Member variables

    // Varient type for storing an input option which can be a string, boolean, double, or integer
    using InputOptionVariant = std::variant<std::string, bool, double, size_t>;

    // InputOptions map that maps an InputOption to an INputOptionsVariant
    std::map<InputOptionType, InputOptionVariant> inputOptions;

    double charge;
    double electronStoppingEnergy;
    double energy;
    bool enableDamageCascade;
    std::string electronScreeningParametersFilename;
    std::string inputDirectoryName;
    double ionDisplacementEnergy;
    double ionStoppingEnergy;
    bool logSingleDisplacement;
    double mass;
    std::string mottScatteringParametersFilename;
    size_t numAngleDivisors;
    size_t numElecScreeningPotentialElements;
    size_t numFlyingDistances;
    size_t numMottScatteringPotentialElements;
    std::string outputFileEndMarker;
    std::string outputFileExtension;
    std::string outputCoordinateFilename;
    std::string outputDirectory;
    double range;
    size_t simulationCount;
    double substrateCharge;
    double substrateDensity;
    double substrateMass;
    ParticleType type;
    bool logEndOfFlyingDistanceOnly;
    bool logStoppingPointOnly;
    bool progressChecking;
    size_t numThreads;
    std::shared_ptr<IOHandler> ioHandler;

    // Functions
    bool parseSetting(
        const std::string& line,    
        std::string& name,
        std::string& value);
    bool stringToBool(const std::string& str);

public:
    // Static methods to get the instance of the class
    static std::shared_ptr<InputFields> getInstance();
    static std::shared_ptr<InputFields> getInstance(std::shared_ptr<IOHandler> ioHandler);  

    // Getters
    bool getEnableDamageCascade() const;
    double getElectronStoppingEnergy() const;
    bool getLogSingleDisplacement() const;
    double getCharge() const;
    double getEnergy() const;
    double getIonDisplacementEnergy() const;
    double getIonStoppingEnergy() const;
    double getMass() const;
    double getRange() const;
    double getSubstrateCharge() const;
    double getSubstrateDensity() const;
    double getSubstrateMass() const;
    ParticleType getType() const;
    size_t getSimulationCount() const;
    size_t getNumAngleDivisors() const;
    size_t getNumElecScreeningPotentialElements() const;
    size_t getNumFlyingDistances() const;
    size_t getNumMottScatteringPotentialElements() const;
    const std::string& getElectronScreeningParametersFilename() const;
    const std::string& getInputDirectoryName() const;
    const std::string& getMottScatteringParametersFilename() const;
    const std::string& getOutputCoordinateFilename() const;
    const std::string& getOutputDirectory() const;
    const std::string& getOutputFileEndMarker() const;
    const std::string& getOutputFileExtension() const;
    bool getLogEndOfFlyingDistanceOnly() const;
    bool getLogStoppingPointOnly() const;
    bool getProgressChecking() const;
    size_t getNumThreads() const;
    std::shared_ptr<IOHandler> getIOHandler() const;

    // Setters
    void setCharge(double charge);
    void setElectronScreeningParametersFilename(const std::string& filename);
    void setElectronStoppingEnergy(double energy);
    void setEnableDamageCascade(bool enable);
    void setEnergy(double energy);
    void setInputDirectoryName(const std::string& name);
    void setIonDisplacementEnergy(double energy);
    void setIonStoppingEnergy(double energy);
    void setLogSingleDisplacement(bool log);
    void setMass(double mass);
    void setMottScatteringParametersFilename(const std::string& filename);
    void setNumAngleDivisors(size_t divisors);
    void setNumElecScreeningPotentialElements(size_t elements);
    void setNumFlyingDistances(size_t distances);
    void setNumMottScatteringPotentialElements(size_t elements);
    void setOutputCoordinateFilename(const std::string& filename);
    void setOutputDirectory(const std::string& directory);
    void setOutputFileEndMarker(const std::string& marker);
    void setOutputFileExtension(const std::string& extension);
    void setRange(size_t range);
    void setSimulationCount(size_t count);
    void setSubstrateCharge(double substrateCharge);
    void setSubstrateDensity(double substrateDensity);
    void setSubstrateMass(double substrateMass);
    void setType(ParticleType type);
    void setLogEndOfFlyingDistanceOnly(bool logEndOfFlyingDistanceOnly);
    void setLogStoppingPointOnly(bool logStoppingPointOnly);
    void setProgressChecking(bool progressChecking);
    void setNumThreads(size_t numThreads);
    void setInputStream(std::istream& inputStream);
    void setIOHandler(std::shared_ptr<IOHandler> ioHandler);

    // Functions
    std::string printInputFields() const;
};


class Bombardment;

class Particle {
    protected:
        // Static member variables
        static std::mt19937 randomGenerator;
        static std::uniform_real_distribution<double> randomDistribution;
        static std::mutex randomGeneratorLock;

        // Member variables
        std::vector<Coordinate> coordinate_vector;
        Coordinate coordinate;
        Velocity velocity;
        double &energy = velocity.energy;
        const double charge;
        const double mass;
        const ParticleType type;
        std::shared_ptr<InputFields> input;
        std::weak_ptr<Bombardment> bombardment;

        // Static functions
        static double random();

    public:
        // Constructors
        Particle(
            const Coordinate& coordinate,
            const Velocity& velocity,
            std::shared_ptr<InputFields> input,
            std::weak_ptr<Bombardment> bombardment);

        Particle(
            const Coordinate& coordinate,
            const Velocity& velocity,
            double mass,
            double charge,
            ParticleType type,
            std::shared_ptr<InputFields> input,
            std::weak_ptr<Bombardment> bombardment); 

        // Destructor
        virtual ~Particle();

        // Modifiers
        virtual void addCoordinate(const Coordinate& coordinate);

        // Accessors
        virtual const Coordinate& getCoordinate() const;
        virtual const std::vector<Coordinate>& getCoordinates() const;
        virtual size_t getNumParticles() const;

        // Member functions
        virtual void fire();
        inline virtual double atomicSpacing() const;
        void relativeToAbsoluteVelocity(
            Velocity& newVelocity,
            Velocity& targetVelocity
            ) const noexcept;
        virtual double sign(const double x) const noexcept;
        virtual inline Coordinate calculateNewCoordinate();
        void createSubstrateKnockon(const Coordinate& newCoordinate,
            const Velocity&
            targetVelocity);

        // Public static functions
        static void seedRandomGenerator();
};  


class Ion : public Particle {
    public:
        // Constructor
        Ion(const Coordinate& coordinate,
            const Velocity& velocity,
            std::shared_ptr<InputFields> input,
            std::weak_ptr<Bombardment> bombardment);

        Ion(const Coordinate& coordinate,
            const Velocity& velocity,
            double charge,
            double mass,
            ParticleType type,
            std::shared_ptr<InputFields> input,
            std::weak_ptr<Bombardment> bombardment); 

        // Actions
        void fire() override;
        double F(double X, double COLUMBIAVK, double AU);
        double DF(double X, double COLUMBIAVK, double AU);
        void recoilEnergyAndVelocity(
            Velocity& newVelocity,
            Velocity& targetVelocity);
        double electronicStoppingEnergy();
        virtual const size_t& getDepth() const;
        inline Coordinate calculateNewCoordinate() override;
};


class Substrate : public Ion {
    private:
        size_t &depth = coordinate.depth;

    public:
        // Constructors
        Substrate(
            const Coordinate& coordinate,
            const Velocity& velocity,   
            std::shared_ptr<InputFields> input,
            std::weak_ptr<Bombardment> bombardment);
        // Functions
        virtual const size_t& getDepth() const override;
};


class Electron : public Particle {
    private:
        static std::vector<
            std::array<std::array<
            double,
            Constants::numMottkParam>,
            Constants::numMottjParam>>
            mottScatteringParams;
        
        static std::vector<
            std::array<
            double,
            Constants::numElecScreeningPotentialParams>> 
            elecScreeningParams;

        static std::vector<double> divisionAngles;

        std::vector<double> randomCrossSections;
        std::vector<double> flyingDistances;
        std::vector<double> scatteringAngles;
        std::vector<double> partialCrossSections;
        std::vector<double> elasticEnergyLoss;

        double correctionFactor;
        double electronStoppingEnergy;
        
    public:
        // Constructor
        Electron(
            const Coordinate& initialCoordinate,
            const Velocity& initialVelocity,
            std::shared_ptr<InputFields> input,
            std::weak_ptr<Bombardment> bombardment);
        
        // Initializer
        static inline double findAngle(
            size_t i,
            std::shared_ptr<InputFields> input);
        static void readParametersAndInitialize(
            std::shared_ptr<InputFields> input);

        // Functions
        void fire();
        double getMottDifferentialCrossSection(double Theta);
        void generateRandomCrossSections();
        double generateFlyingDistances(double totalMottCrossSection);
        double getMottTotalCrossSection();
        void generateIntegratedCrossSections(double totalMottCrossSection);
        double getElasticEnergyLoss();
        double getIonizationEnergyLoss();
        inline double getBremsstrahlung();
        inline Coordinate calculateNewElectronCoordinate(size_t i);
};

class Simulation;

class Bombardment : public std::enable_shared_from_this<Bombardment> {
    private:
        std::vector<std::unique_ptr<Particle>> particles;
        std::weak_ptr<Simulation> simulation;
        size_t id;
    
    public:
        // Constructors
        Bombardment(std::weak_ptr<Simulation> simulation, size_t id);

        // Destructor
        ~Bombardment();

        // Accessor
        const std::vector<std::unique_ptr<Particle>>& getParticles() const;

        void initiate(std::shared_ptr<InputFields> inputFields);
        void addParticle(std::unique_ptr<Particle> particle);
        void popParticle();
};

class Simulation : public std::enable_shared_from_this<Simulation> {
    private:
        size_t bombardmentCount;
        std::shared_ptr<InputFields> input;
        std::vector<std::shared_ptr<Bombardment>> bombardments;
        fs::path outputPath;
        std::mutex fileLock;
        std::mutex counterLock;
        std::vector<std::thread> threads;
        u_int64_t progressCounter;

    public:
        // Constructors
        Simulation();
        Simulation(std::shared_ptr<InputFields> input);

        // Destructors
        ~Simulation();

        // Accessors
        std::shared_ptr<InputFields> getInput();
        size_t getBombardmentSize();
        size_t getBombardmentCapacity();

        // Modifiers
        void addParticle();

        // Functions
        void initiate();

        void writeData(
            OutputType outputType,
            const std::vector<std::unique_ptr<Particle>>& particles = {},
            size_t simulationID = 0
            );
        bool fileIsWritten(const fs::path filename);
        void renameFileWithTimestamp(const std::string& filename);
        std::string getCurrentDateTime();
        void checkOutputFiles();
        void updateProgressCounter();
};

#endif

/******************************************************************************
 * Filename:        constants.cpp
 * Project:         AMCSET
 * Author:          Nathaniel Thomas
 * Date Created:    March 8, 2024
 * Date Modified:   April 30, 2024
 * File Version:    1.0
 * Group:           Dr. Shao's RMSLCF Group
 *
 * Description:
 * This file contains useful constants to be used in the AMCSET program.
 * 
 *
 * Revision History
 * 1.0:
 * - File created and constants added
 * See main.cpp for other revision history
 * 
 *****************************************************************************/
#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

enum OutputType{
    COORDINATE,
    ENDOFFILE
};

enum ParticleType{
    ION,
    ELECTRON,
    SUBSTRATE
};

// Options to add
// electronMode
// electronEnergy
// electronStoppingEnergy
// electronScreeningParametersFilename
// numElecScreeningPotentialElements
// mottScatteringParametersFilename
// numMottScatteringPotentialElements
// numAngleDivisors
// numFlyingDistances
// enableDamageCascade
// ionCharge
// ionEnergy
// ionMass
// substrateDisplacementEnergy
// ionStoppingEnergy
// logSingleDisplacement
// inputDirectoryName
// outputCoordinateFilename
// outputFileEndMarker
// outputFileExtension
// outputDirectory
// substrateCharge
// substrateMass
// logEndOfFlyingDistanceOnly
// logStoppingPointOnly
// simulationCount
// numThreads

enum InputOptionType{
    ELECTRON_MODE,
    ELECTRON_ENERGY,
    ELECTRON_STOPPING_ENERGY,
    ELECTRON_SCREENING_PARAMETERS_FILENAME,
    NUM_ELEC_SCREENING_POTENTIAL_ELEMENTS,
    MOTT_SCATTERING_PARAMETERS_FILENAME,
    NUM_MOTT_SCATTERING_POTENTIAL_ELEMENTS,
    NUM_ANGLE_DIVISORS,
    NUM_FLYING_DISTANCES,
    ENABLE_DAMAGE_CASCADE,
    ION_CHARGE,
    ION_ENERGY,
    ION_MASS,
    SUBSTRATE_DISPLACEMENT_ENERGY,
    ION_STOPPING_ENERGY,
    LOG_SINGLE_DISPLACEMENT,
    INPUT_DIRECTORY_NAME,
    OUTPUT_COORDINATE_FILENAME,
    OUTPUT_FILE_END_MARKER,
    OUTPUT_FILE_EXTENSION,
    OUTPUT_DIRECTORY,
    SUBSTRATE_CHARGE,
    SUBSTRATE_MASS,
    LOG_END_OF_FLYING_DISTANCE_ONLY,
    LOG_STOPPING_POINT_ONLY,
    SIMULATION_COUNT,
    NUM_THREADS
};

namespace Constants {
    constexpr double avogadrosNumber = 6.02214076e23; // atoms/mol
    constexpr double cmPerAngstrom = 1e-8; // cm/Angstrom
    constexpr double electronCharge = -1; // e
    constexpr double electronChargeSI = 1.602176634e-19; // C
    constexpr double electronMass = 9.1093837015E-31; // kg
    constexpr double keVToeV = 1000; // eV/keV
    constexpr size_t numElecScreeningPotentialParams = 6;
    constexpr size_t numMottjParam = 5; 
    constexpr size_t numMottkParam = 6;
    constexpr double plancksConstant = 6.62607015e-34; // J·s
    constexpr double speedOfLight = 299792458; // m/s

    // Dependent values
    constexpr double mec2 = electronMass
                            *speedOfLight*speedOfLight
                            /electronChargeSI
                            /keVToeV; // keV
}

namespace Defaults {
    constexpr double electronEnergy = 10000; // keV, default electron energy
    constexpr const char* electronScreeningParametersFilename = "electron_screeening_potentials.csv";
    constexpr double electronStoppingEnergy = 1; // keV
    constexpr bool enableDamageCascade = false;
    constexpr double ionCharge = 14; // e, silicon
    constexpr double ionDisplacementEnergy = 0.04; // keV
    constexpr double ionStoppingEnergy = 0.04; // keV
    constexpr double ionEnergy = 50; // keV, default ion energy
    constexpr double ionMass = 27.97692653442; // amu, silicon 29
    constexpr double range = 80000000; // angstrom
    constexpr bool logEndOfFlyingDistanceOnly = false;
    constexpr bool logSingleDisplacement = false; // Disable to shrink output file size
    constexpr bool logStoppingPointOnly = true;
    constexpr size_t numAngleDivisors = 1000;
    constexpr size_t numElecScreeningPotentialElements = 92;
    constexpr size_t numFlyingDistances = 1000;
    constexpr size_t numMottElements = 118; 
    constexpr const char* inputDirectory =  "input";
    constexpr const char* mottScatteringParametersFilename = "mott_scattering_parameters.csv";
    constexpr const char* outputCoordinateFilename =  "coordinateOutput";
    constexpr const char* outputFileEnd = "End of file";
    constexpr const char* outputFileExtension = ".csv";
    constexpr const char* outputDirectory =  "output";
    constexpr ParticleType type = ION;
    constexpr double substrateCharge = 26; // e, iron
    constexpr double substrateMass = 55.9349363; // amu, iron 56
    constexpr const char* settingsFilename = "settings.txt";
    constexpr size_t simulationCount = 1;
    constexpr bool progressChecking = false;
    constexpr size_t numThreads = 1;

    // Dependent values
    constexpr double windowRange = Defaults::range;
    constexpr double substrateDensity = 7.874
                                        /Defaults::substrateMass
                                        *Constants::avogadrosNumber
                                        *(   Constants::cmPerAngstrom
                                            *Constants::cmPerAngstrom
                                            *Constants::cmPerAngstrom
                                        ); // g/cm^3 -> atoms/angstrom^3
}

#endif

/******************************************************************************
 * Filename:        utilities.h
 * Project:         AMCSET
 * Author:          Nathaniel Thomas
 * Date Created:    March 12, 2024
 * Date Modified:   April 30, 2024
 * File Version:    1.0
 * Group:           Dr. Shao's RMSLCF Group
 *
 * Description:
 * This file contains useful utilities and constants used by the 
 * main function in the AMCSET program
 * 
 *
 * Revision History
 * 1.0:
 * - File created
 * See main.cpp for other revision history
 *****************************************************************************/
#ifndef UTILITIES_H
#define UTILITIES_H

// Library includes
#include <cstdlib>
#include <cstddef>
#include <iostream>

// Local includes
#include "AMCSET_classes.h"

// Define a struct to store settings
struct Arguments {
    std::string filename;
    bool time;
    bool displaySettings;
    bool progress;
};

// Define a custom exception that carries an exit status
class ExitException : public std::exception {
public:
    ExitException(int status);
    int status;
};

class InputFields;


// This class is for handler USER IO, output file data is handled by the Simulation::writeData() function
class IOHandler {
    private:
        // Private constructor to prevent external instantiation
        IOHandler(
            std::shared_ptr<InputFields> input,
            std::istream& inputStream,
            std::ostream& outputStream,
            std::ostream& errorStream);

        // Static instance pointer for a single IOHandler instance
        static std::shared_ptr<IOHandler> instance;

        // Constants

        const std::string_view inputNotInitializedMessage = "InputFields input has not been initialized.";

        const std::string_view helpMessage = 
        R"(AMCSET Version 1.0
        Author:         Nathaniel Thomas
        Contact:        nathaniel@swbell.net
        Release date:   March 12, 2024

        This is a program for simulating ion and electron bombardment in materials.

        Usage: ./AMCSET [options][paths...]

        Options
            -f --filename <filename>    Read settings for AMCSET from <filename>. [Default="settings.txt"]
            -t --time                   Record execution time and output to standard output.
            -s --settings               Display an example settings file.   
            -d --display                Output the active settings to the standard output.
            -h --help                   Display this help message.
            -p --progress               Display the progress of the simulation while it is running to the standard output.
        )";

        std::string settingsMessage;

        // Private Members
        std::shared_ptr<InputFields> input;
        std::istream& inputStream;
        std::ostream& outputStream;
        std::ostream& errorStream;
        Arguments arguments;
        bool argumentsSet = false;

    public:

        // Instantiators
        static std::shared_ptr<IOHandler> getInstance();  
        
        static std::shared_ptr<IOHandler> getInstance(
            std::shared_ptr<InputFields> input,
            std::istream& inputStream = std::cin,
            std::ostream& outputStream = std::cout,
            std::ostream& errorStream = std::cerr);

        // Getters
        std::shared_ptr<InputFields> getInput() const;
        Arguments& getArguments();
        std::istream& getInputStream() const;
        std::ostream& getOutputStream() const;
        std::ostream& getErrorStream() const;

        // Setters
        void setInput(std::shared_ptr<InputFields> input);
        // No setters for reference members, as they cannot be rebound.

        // Functions
        void parseCommandLine(const int argc, const char* argv[]); // Add scope resolution operator
        bool promptContinue();
        void clearLine();
        void checkHardwareThreads();
        void checkDisplayOption(); // Add scope resolution operator
        void writeSettingsToFile(const std::string& settingsFilename);
        bool parseSetting(const std::string& line, std::string& name, std::string& value);
        bool stringToBool(const std::string& str);
        void handleFileOpenError(const std::string& filename);
        void readSettingsFromFile();
        void handleNoOutputDirectory(const std::filesystem::path& outputPath);
        void handleOutputFileOpenError(
            OutputType outputType,
            std::filesystem::path outputPath,
            std::string filename,
            int bombardmentID);
        bool areArgumentsSet() const;

        // Templated function
        template<typename T>
        IOHandler& operator<<(const T& message) {
            errorStream << message;
            return *this;
        }

        IOHandler& operator<<(std::ostream& (*func)(std::ostream&)) {
            func(errorStream);
            return *this;
        }
};

#endif

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

/******************************************************************************
 * Filename:        AMCSET_classes.cpp
 * Project:         AMCSET
 * Author:          Nathaniel Thomas
 * Date Created:    March 8, 2024
 * Date Modified:   April 30, 2024
 * File Version:    1.0
 * Group:           Dr. Shao's RMSLCF Group
 *
 * Description:
 * This file contains the class declarations for the AMCSET program. 
 * implementations are in AMCSET_classes.cpp
 * 
 *
 * Revision History
 * 1.0:
 * - File created and class implementations added
 *
 *****************************************************************************/
#ifndef AMCSET_CLASSES_CXX
#define AMCSET_CLASSES_CXX

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
#include "AMCSET_classes.h"
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
    charge(input->getCharge()),
    mass(input->getMass()),
    type(input->getType()),
    input(input),
    bombardment(bombardment) {};

Particle::Particle(
    const Coordinate& coordinate,
    const Velocity& velocity,
    double charge,
    double mass,
    ParticleType type,
    std::shared_ptr<InputFields> input,
    std::weak_ptr<Bombardment> bombardment)
    : coordinate(coordinate),
    velocity(velocity),
    charge(charge),
    mass(mass),
    type(type),
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
    return std::pow(input->getSubstrateDensity(), -1.0/3.0);
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
    double charge,
    double mass,
    ParticleType type,
    std::shared_ptr<InputFields> input,
    std::weak_ptr<Bombardment> bombardment)
    : Particle(
        coordinate,
        velocity,
        charge,
        mass,
        type,
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
    while(energy > input->getIonStoppingEnergy()) {
        velocity.energy -= electronicStoppingEnergy();
        recoilEnergyAndVelocity(newVelocity, targetVelocity);
        velocity = newVelocity;
        coordinate = calculateNewCoordinate();
        if(coordinate.z < 0.0) {
            DEBUG_PRINT("Ion Sputter");
            return;
        }
        if(input->getEnableDamageCascade()
                && energy - targetVelocity.energy > input->getIonDisplacementEnergy()) {
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
    const double P = std::sqrt(random()/(M_PI*std::pow(input->getSubstrateDensity(), 2.0/3.0)));
    const double COLUMBIAVK = 0.0143992*charge*input->getSubstrateCharge();
    const double MC = mass*input->getSubstrateMass()/(mass+input->getSubstrateMass());
    const double INVLAB = std::sqrt(energy*2.0/mass);
    const double EC = 0.5*MC*INVLAB*INVLAB;
    const double AU = 0.8854*0.529/std::pow(std::sqrt(mass)
                +std::sqrt(input->getSubstrateCharge()), 2.0/3.0);
    const double ELINHARD = EC*AU/COLUMBIAVK;
    // Find RMIN for different energy
    double AA = P*P;
    if(AA == 0) {
        AA = .00001;
    }
    const double BB = COLUMBIAVK/EC;
    const double CC = -1;
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
    const double RMIN = 0.5*(RMINTRY2+RMINTRY1);
    // Calculate deflection angle
    const double RBIERSACK = 2.0*(EC-F(RMIN, COLUMBIAVK, AU))/RMIN*RMIN
        /DF(RMIN, COLUMBIAVK, AU);
    const double BBIERSACK = P/AU;
    const double ROBIERSACK = 1.0/(RMIN*AU);
    const double RCBIERSACK = RBIERSACK/AU;
    const double C1BIERSACK = 0.6743;
    const double C2BIERSACK = 0.009611;
    const double C3BIERSACK = 0.005175;
    const double C4BIERSACK = 10.0;
    const double C5BIERSACK = 6.314;
    const double ALTHABIERSACK = 1.0+1.0/std::sqrt(C1BIERSACK*ELINHARD);
    const double BELTABIERSACK = (C2BIERSACK+std::sqrt(ELINHARD))
        /(C3BIERSACK+std::sqrt(ELINHARD));
    const double GAMABIERSACK = (C4BIERSACK+ELINHARD)/(C5BIERSACK+ELINHARD);
    const double ABIERSACK = 2.0*ALTHABIERSACK*ELINHARD
        *std::pow(BBIERSACK, BELTABIERSACK);
    const double GBIERSACK = GAMABIERSACK*1.0
        /(std::sqrt((1.0+ABIERSACK*ABIERSACK))-ABIERSACK);
    const double DELTABIERSACK = ABIERSACK*(ROBIERSACK-BBIERSACK)/(1.0+GBIERSACK);
    double CALPHA1;
    if(P != 0.0) {
        CALPHA1 = 2.0*std::atan(std::sqrt((ROBIERSACK+RCBIERSACK)
            *(ROBIERSACK+RCBIERSACK)/((BBIERSACK+RCBIERSACK+DELTABIERSACK)
            *(BBIERSACK+RCBIERSACK+DELTABIERSACK))-1.0));
    } else {
        CALPHA1 = M_PI;
    }
    // Defelction angle in lab coordinator
    const double COSPLUSMASS = std::cos(CALPHA1)+mass/input->getSubstrateMass();
    double THETA1RELATIVE;
    if(COSPLUSMASS > 0.0) {
        THETA1RELATIVE = std::atan(std::sin(CALPHA1)/COSPLUSMASS);
    } else if (COSPLUSMASS < 0.0) {
        THETA1RELATIVE = M_PI+std::atan(std::sin(CALPHA1)/COSPLUSMASS);
    } else {
        THETA1RELATIVE = M_PI/2.0;
    }

    // Calculate target direction and energy
    const double RE = 4.0*EC*MC/input->getSubstrateMass()*std::sin(CALPHA1/2.0)*std::sin(CALPHA1
                /2.0);
    const double THETA2RELATIVE = (M_PI-CALPHA1)/2.0;

    newVelocity.zAngle = THETA1RELATIVE;
    newVelocity.energy = energy - RE;
    targetVelocity.zAngle = THETA2RELATIVE;
    targetVelocity.energy = RE;
    
    relativeToAbsoluteVelocity(newVelocity, targetVelocity);
};

double Ion::electronicStoppingEnergy() {
    double stoppingEnergy = 1.212*std::pow(charge, 7.0/6.0)*input->getSubstrateCharge()
                /(std::pow((std::pow(charge, 2.0/3.0) \
                + std::pow(input->getSubstrateCharge(), 2.0/3.0))
                , 3.0/2.0)*std::sqrt(mass))*std::sqrt(energy*1000);
    double energyLoss = 1.59*atomicSpacing()*input->getSubstrateDensity()*stoppingEnergy/1000;
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
        input->getSubstrateCharge(),
        input->getSubstrateMass(),   
        SUBSTRATE,
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
            coordinate = calculateNewElectronCoordinate(i);

            if(coordinate.z < 0.0) {
                DEBUG_PRINT("Electron Sputter");
                return;
            }

            if(input->getEnableDamageCascade()
                && targetVelocity.energy > input->getIonDisplacementEnergy()) {
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
    const double DSC = (input->getSubstrateCharge()*re1)*(input->getSubstrateCharge()*re1)*(1.0-betaSquared)/
        betaQuarted/(oneMinusCos*oneMinusCos);
    const double moment = 2.0*beta*gamma*elecmass*Vc*std::sin(Theta/2.0);
    const double momentSquared = moment*moment;
    
    double Rmott = 0.0;
    double alpha1;
    for(size_t j = 0; j < Constants::numMottjParam; j++) {
        alpha1 = 0.0;
        for(size_t k = 0; k < Constants::numMottkParam; k++) {
            alpha1 = alpha1 +
                mottScatteringParams[static_cast<size_t>(input->getSubstrateCharge())][j][k]
                *std::pow(beta-betaAVE, k);
        }
        Rmott = Rmott+alpha1*std::pow(std::sqrt(oneMinusCos), j);
    }
    
    double Feq = 0;
    for(size_t i = 0; i < Constants::numElecScreeningPotentialParams/2; i++) {
        Feq = Feq + elecScreeningParams[input->getSubstrateCharge()][i]*
            elecScreeningParams[input->getSubstrateCharge()][i + 3]
            *elecScreeningParams[input->getSubstrateCharge()][i + 3]
            /(elecScreeningParams[input->getSubstrateCharge()][i + 3]
            *elecScreeningParams[input->getSubstrateCharge()][i + 3]
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
    double atomicSpacing = std::pow(input->getSubstrateDensity()*1e24, -1.0/3.0);

    for(size_t i = 0; i < input->getNumFlyingDistances(); i++) {
        crossSectionFlight = -std::log(random())
            /(input->getSubstrateDensity()*1e24*totalMottCrossSection);
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
                *(std::sin(scatteringAngles[i]))+input->getSubstrateMass()*931.0*1000.0*(1
                -std::cos(scatteringAngles[i])))
            *energy*(energy+2.0*511.0);
        double ElecREdown = (energy+input->getSubstrateMass()*931.5*1000.0)
            *(energy+input->getSubstrateMass()*931.5*1000.0)-energy*(energy+2.0*511.0)
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
    double KforLoss = 0.734*std::pow(input->getSubstrateCharge(), 0.037);
    double JforLoss = (input->getSubstrateCharge() < 13) ? 0.0115*input->getSubstrateCharge() :
                        0.00976*input->getSubstrateCharge()+0.0585
                        *std::pow(input->getSubstrateCharge(), -0.19);

    JforLoss = JforLoss/(1.0+KforLoss*JforLoss/correctedEnergy);
    double ionizationEnergy = 153.55*1E24*input->getSubstrateDensity()*input->getSubstrateCharge()
        /(Constants::avogadrosNumber)/(beta*beta)*(std::log(511.0*beta*beta
        *correctedEnergy/(2.0*JforLoss*JforLoss*(1.0-beta*beta)))-std::log(2.0)
        *(2.0*std::sqrt(1.0-beta*beta)-1+beta*beta)+1-beta*beta+1.0/8.0
        *(1.0-std::sqrt(1.0-beta*beta))*(1.0-std::sqrt(1.0-beta*beta)));
    // Final expression ensures a positive value
    return (ionizationEnergy + std::abs(ionizationEnergy))/2.0;
};

inline double Electron::getBremsstrahlung() {
    double correctedEnergy = energy*correctionFactor;
    double bremsstralung = 1.4e-4*1E24*input->getSubstrateDensity()/(Constants::avogadrosNumber)
        *input->getSubstrateCharge()*(input->getSubstrateCharge()+1.0)*(correctedEnergy+511.0)*(4.0
        *std::log(2.0*(correctedEnergy+511.0)/511.0)-4.0/3.0);
    // Final expressio ensures a positive value
    return (bremsstralung + std::abs(bremsstralung))/2.0;
};

inline Coordinate Electron::calculateNewElectronCoordinate(size_t i) {
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
    type(Defaults::type),
    logEndOfFlyingDistanceOnly(Defaults::logEndOfFlyingDistanceOnly),
    logStoppingPointOnly(Defaults::logStoppingPointOnly),
    progressChecking(Defaults::progressChecking),
    numThreads(Defaults::numThreads),
    ioHandler(IOHandler::getInstance()) {}


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
    ParticleType type,
    bool logEndOfFlyingDistanceOnly,
    bool logStoppingPointOnly,
    bool progressChecking,
    size_t numThreads,
    std::shared_ptr<IOHandler> ioHandler)
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
    type(type),
    logEndOfFlyingDistanceOnly(logEndOfFlyingDistanceOnly),
    logStoppingPointOnly(logStoppingPointOnly),
    progressChecking(progressChecking),
    numThreads(numThreads),
    ioHandler(ioHandler) {};

// Using the delegated constructor allows default values to be used if they are not specified in the settings.
InputFields::InputFields(std::shared_ptr<IOHandler> ioHandler)
: InputFields()
{
    this->ioHandler = ioHandler;
};

std::shared_ptr<InputFields> InputFields::instance = nullptr;

std::shared_ptr<InputFields> InputFields::getInstance(std::shared_ptr<IOHandler> ioHandler) {
    // Create the instance if it doesn't exist
    if (!instance) {
        instance = std::shared_ptr<InputFields>(new InputFields(ioHandler));
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
bool InputFields::getLogEndOfFlyingDistanceOnly() const { return logEndOfFlyingDistanceOnly;};
bool InputFields::getLogStoppingPointOnly() const { return logStoppingPointOnly;};
bool InputFields::getProgressChecking() const { return progressChecking; };
size_t InputFields::getNumThreads() const { return numThreads; };
std::shared_ptr<IOHandler> InputFields::getIOHandler() const { return ioHandler; };

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
void InputFields::setSimulationCount(size_t count) { this->simulationCount = count; }
void InputFields::setSubstrateCharge(double substrateCharge) { this->substrateCharge = substrateCharge; }
void InputFields::setSubstrateDensity(double substrateDensity) { this->substrateDensity = substrateDensity; }
void InputFields::setSubstrateMass(double substrateMass) { this->substrateMass = substrateMass; }
void InputFields::setType(ParticleType type) { this->type = type; }
void InputFields::setLogEndOfFlyingDistanceOnly(bool logEndOfFlyingDistanceOnly) {this->logEndOfFlyingDistanceOnly = logEndOfFlyingDistanceOnly; }
void InputFields::setLogStoppingPointOnly(bool logStoppingPointOnly) {this->logStoppingPointOnly = logStoppingPointOnly; }
void InputFields::setProgressChecking(bool progressChecking) {this->progressChecking = progressChecking; }
void InputFields::setNumThreads(size_t numThreads) {this->numThreads = numThreads; }
void InputFields::setIOHandler(std::shared_ptr<IOHandler> ioHandler) {this->ioHandler = ioHandler; }

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
    // oss << "Range(unused):\t" << range << std::endl;
    oss << "Log Stopping Point Only:\t" << (logStoppingPointOnly ? "true" : "false") << std::endl;
    oss << "Log End of Flying Distance Only:\t" << (logEndOfFlyingDistanceOnly ? "true" : "false") << std::endl;
    oss << "Number of Threads:\t" << numThreads << std::endl;
    oss << "Simulation Count:\t" << simulationCount << std::endl;
    oss << "Progress Tracking:\t" << (progressChecking ? "true" : "false") << std::flush;

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
        *(input->getIOHandler()) << "Failed to initiate write to " << 
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
                if (input->getProgressChecking()) {
                    updateProgressCounter();
                }
                bombardments[i]->initiate(input);
                bombardments[i].reset();
            });
            bombardmentCount++;
        }

        for (auto& thread : batchThreads) {
            if (thread.joinable()) {
                thread.join();
            }
        }
    }

    #ifndef NO_OUTPUT
    writeData(ENDOFFILE);
    #endif
}

void Simulation::writeData(
    OutputType outputType,
    const std::vector<std::unique_ptr<Particle>>& particles,
    size_t bombardmentID
    ) {
    
    // Acquire lock
    std::lock_guard<std::mutex> lock(fileLock);

    DEBUG_PRINT("Simulation " << bombardmentID << " Writing");

    // Check if output directory exists, create it if not
    if (!fs::exists(outputPath)) {
        input->getIOHandler()->handleNoOutputDirectory(outputPath);
    }

    // Open output file
    std::string filename = std::string(input->getOutputCoordinateFilename());
    filename = filename + std::string(input->getOutputFileExtension());
    std::ofstream outputFile(outputPath / filename, std::ios_base::app);
    if (!outputFile.is_open()) {
        input->getIOHandler()->handleOutputFileOpenError(
            outputType,
            outputPath,
            filename,
            bombardmentID);;
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
            *(input->getIOHandler()) << "Unsupported output type." << std::endl;
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

        *(input->getIOHandler()) << "Progress: [";
        for (int i = 0; i < progressBarWidth; ++i) {
            if (i < progress) {
                *(input->getIOHandler()) << "=";
            } else {
                *(input->getIOHandler()) << " ";
            }
        }
        *(input->getIOHandler()) << "] "
                    << std::setw(3)
                    << std::setfill(' ')
                    << std::right
                    << ((progressCounter * 100) / input->getSimulationCount())
                    << "%\r"
                    << std::flush;

}

#endif

/******************************************************************************
 * Filename:        main.cpp
 * Project:         AMCSET
 * Author:          Nathaniel Thomas
 * Date Created:    March 15, 2024
 * Date Modified:   April 30, 2024
 * File Version:    1.0
 * Group:           Dr. Shao's RMSLCF Group
 *
 * Description:
 * This is the main file for the AMCSET code developed by Dr. Shao and converted
 * to the C++ programming language
 *
 * Revision History
 * 1.0:
 * - File created to account for testing requirements
 *****************************************************************************/
#include "AMCSET.h"

// Main program  start
int main(const int argc, const char* argv[]) {
    return startAMCSET(argc, argv);
}

/******************************************************************************
 * Filename:        utilities.cpp
 * Project:         AMCSET
 * Author:          Nathaniel Thomas
 * Date Created:    March 12, 2024
 * Date Modified:   April 30, 2024
 * File Version:    1.0
 * Group:           Dr. Shao's RMSLCF Group
 *
 * Description:
 * This file contains useful utilities and constants used by the 
 * main function in the AMCSET program
 * 
 *
 * Revision History
 * 1.0:
 * - File created
 * See main.cpp for other revision history
 *
 *****************************************************************************/
// Compiler directives
#ifndef UTILITIES_CXX
#define UTILITIES_CXX

// Library includes
#include <cstdlib>
#include <cstddef>
#include <iostream>
#include <fstream>

// Local includes
#include "constants.h"
#include "utilities.h"

std::shared_ptr<IOHandler> IOHandler::instance = nullptr;

IOHandler::IOHandler(
        std::shared_ptr<InputFields> input,
        std::istream& inputStream,
        std::ostream& outputStream,
        std::ostream& errorStream)
        : input(input),
        inputStream(inputStream),
        outputStream(outputStream),
        errorStream(errorStream),
        arguments({}), 
        argumentsSet(false) {
    std::stringstream settingsMessageStream;
    settingsMessageStream << "electronMode=" << ((Defaults::type == ELECTRON) ? "true" : "false") << std::endl;
    settingsMessageStream << "electronEnergy(keV)=" << Defaults::electronEnergy << std::endl;
    settingsMessageStream << "electronStoppingEnergy(keV)=" << Defaults::electronStoppingEnergy << std::endl;
    settingsMessageStream << "electronScreeningParametersFilename=\"" << Defaults::electronScreeningParametersFilename << "\"" << std::endl;
    settingsMessageStream << "numElecScreeningPotentialElements=" << Defaults::numElecScreeningPotentialElements << std::endl;
    settingsMessageStream << "mottScatteringParametersFilename=\"" << Defaults::mottScatteringParametersFilename << "\"" << std::endl;
    settingsMessageStream << "numMottScatteringPotentialElements=" << Defaults::numMottElements << std::endl;
    settingsMessageStream << "numAngleDivisors=" << Defaults::numAngleDivisors << std::endl;
    settingsMessageStream << "numFlyingDistances=" << Defaults::numFlyingDistances << std::endl;
    settingsMessageStream << "enableDamageCascade=" << ((Defaults::enableDamageCascade) ? "true" : "false") << std::endl;
    settingsMessageStream << "ionCharge(e)=" << Defaults::ionCharge << std::endl;
    settingsMessageStream << "ionEnergy(keV)=" << Defaults::ionEnergy << std::endl;
    settingsMessageStream << "ionMass(amu)=" << Defaults::ionMass << std::endl;
    settingsMessageStream << "substrateDisplacementEnergy(keV)=" << Defaults::ionDisplacementEnergy << std::endl;
    settingsMessageStream << "ionStoppingEnergy(keV)=" << Defaults::ionStoppingEnergy << std::endl;
    settingsMessageStream << "logSingleDisplacement=" << ((Defaults::logSingleDisplacement) ? "true" : "false") << std::endl;
    settingsMessageStream << "inputDirectoryName=\"" << Defaults::inputDirectory << "\"" << std::endl;
    settingsMessageStream << "outputCoordinateFilename=\"" << Defaults::outputCoordinateFilename << "\"" << std::endl;
    settingsMessageStream << "outputFileEndMarker=\"" << Defaults::outputFileEnd << "\"" << std::endl;
    settingsMessageStream << "outputFileExtension=\"" << Defaults::outputFileExtension << "\"" << std::endl;
    settingsMessageStream << "outputDirectory=\"" << Defaults::outputDirectory << "\"" << std::endl;
    settingsMessageStream << "substrateCharge(e)=" << Defaults::substrateCharge << std::endl;
    settingsMessageStream << "substrateMass(amu)=" << Defaults::substrateMass << std::endl;
    settingsMessageStream << "logEndOfFlyingDistanceOnly=" << ((Defaults::logEndOfFlyingDistanceOnly) ? "true" : "false") << std::endl;
    settingsMessageStream << "logStoppingPointOnly=" << ((Defaults::logStoppingPointOnly) ? "true" : "false") << std::endl;
    settingsMessageStream << "simulationCount=" << Defaults::simulationCount << std::endl;
    settingsMessageStream << "numThreads=" << Defaults::numThreads;
    settingsMessage = settingsMessageStream.str();
};

std::shared_ptr<IOHandler> IOHandler::getInstance()
{
    if(!instance) {
        throw std::logic_error("IOHandler instance has not been initialized.");
    }
    return instance;
}

std::shared_ptr<IOHandler> IOHandler::getInstance(
    std::shared_ptr<InputFields> input,
    std::istream& inputStream,
    std::ostream& outputStream,
    std::ostream& errorStream)
{
    if(!instance) {
        instance = std::shared_ptr<IOHandler>(new IOHandler(
            input,
            inputStream,
            outputStream,
            errorStream));
    }
    return instance;
}

// Getters
std::shared_ptr<InputFields> IOHandler::getInput() const { return input; }
Arguments& IOHandler::getArguments() { return arguments; }
std::istream& IOHandler::getInputStream() const { return inputStream; }
std::ostream& IOHandler::getOutputStream() const { return outputStream; }
std::ostream& IOHandler::getErrorStream() const { return errorStream; }

// Setters
void IOHandler::setInput(std::shared_ptr<InputFields> input) { this->input = input; }

// Function to parse command line arguments
void IOHandler::parseCommandLine(const int argc, const char* argv[]) {

    DEBUG_PRINT("-----------------------Progress Checking 116-----------------------");

    Arguments& settings = arguments;
    
    settings.filename = "settings.txt"; // Default filename
    settings.time = false; // Default time flag
    settings.displaySettings = false; // Default time flag
    settings.progress = false; // Default time flag

    // No arguments passed
    if(argc == 1) {
        DEBUG_PRINT("-----------------------Progress Checking 127-----------------------");
        return;
    }


    // Check for help option
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            outputStream << helpMessage << std::endl;
            throw ExitException(EXIT_SUCCESS);
        }
    }

    DEBUG_PRINT("-----------------------Progress Checking 141-----------------------");

    // Check for settings display option
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-s" || arg == "--settings") {
            outputStream << "Settings\n\n"
                    << "The settings file is a text file, usually named \"settings.txt\", "
                    << "but you can use a custom settings filename if you pass it via the command line."
                    << "Here is an example file of the default settings:\n\n"
                    << settingsMessage
                    << std::endl;
            throw ExitException(EXIT_SUCCESS);
        }
    }

    DEBUG_PRINT("-----------------------Progress Checking 157-----------------------");
    
    // Parse other options
    for (int i = 1; i < argc; ++i) {
        DEBUG_PRINT("Val: " << argc);
        std::string arg = argv[i];
        if (arg == "-f" || arg == "--filename") {
            if (i + 1 < argc) {
                settings.filename = argv[i + 1];
                ++i; // Skip next argument
            } else {
                errorStream << "Error: Missing argument for filename option" << std::endl;
                throw ExitException(EXIT_FAILURE);
            }
        } else if (arg == "-t" || arg == "--time") {
            settings.time = true;
        } else if (arg == "-d" || arg == "--display") {
            settings.displaySettings = true;
        } else if (arg == "-p" || arg == "--progress") {
            settings.progress = true;
        } else {
            errorStream << "Error: Unrecognized flag '" << arg << "'" << std::endl;
            errorStream << "Use -h or --help for usage information" << std::endl;
            throw ExitException(EXIT_FAILURE);
        }
    }

    argumentsSet = true;
}

bool IOHandler::promptContinue() {
    char answer;

    while(true) {
        outputStream << "Continue? (y/n) ";
        inputStream >> answer;

        // Clear input buffer to handle incorrect input
        inputStream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        if (answer == 'y' || answer == 'Y') {
            return true;
        } else if (answer == 'n' || answer == 'N') {
            return false;
        } else {
            outputStream << "Invalid input. Please enter 'y' or 'n'." << std::endl;
        }
    }
}

void IOHandler::clearLine() {
    constexpr int width = 120;
    outputStream << "\r";
    for (size_t i = 0; i < width; ++i) {
        outputStream << " ";
    }
    outputStream << "\r";
    outputStream.flush();
}

void IOHandler::checkHardwareThreads()
{
    if(input == nullptr) {
        throw std::logic_error(std::string(inputNotInitializedMessage));
    }

    DEBUG_PRINT("-----------------------Progress Checking 223-----------------------");


    size_t cores = std::thread::hardware_concurrency();
    DEBUG_PRINT("-----------------------Progress Checking 227-----------------------");
    size_t threads = input->getNumThreads();
    
    if (cores < threads)
    {
        DEBUG_PRINT("cores: " << cores << " threads: " << threads);
        DEBUG_PRINT("-----------------------Progress Checking 232-----------------------");

        const std::string ANSI_COLOR_GREEN = "\033[1;33m";
        const std::string ANSI_COLOR_RESET = "\033[0m";
        errorStream << ANSI_COLOR_GREEN
                << "Warning: Number of threads ("
                << threads
                << ") exceeds number of hardware threads ("
                << cores
                << "). Continuing may cause unexpected behavior on this system."
                << " Consider settings the \"numThreads\" setting to a value less"
                << " than or equal to "
                << cores
                << "."
                << ANSI_COLOR_RESET
                << "\nUse -s for settings help. Use -h for usage help."
                << std::endl;
        if (!promptContinue())
        {
            DEBUG_PRINT("-----------------------Progress Checking 251-----------------------");

            throw ExitException(EXIT_SUCCESS);
        }
    }
    DEBUG_PRINT("-----------------------Progress Checking 256-----------------------");

}

void IOHandler::checkDisplayOption()
{
    if(input == nullptr) {
        throw std::logic_error(std::string(inputNotInitializedMessage));
    }
    if (arguments.displaySettings)
    {
        const std::string ANSI_COLOR_GREEN = "\033[1;32m";
        const std::string ANSI_COLOR_RESET = "\033[0m";
        outputStream
            << ANSI_COLOR_GREEN
            << "Stimulation settings:\n\n"
            << ANSI_COLOR_RESET
            << input->printInputFields()
            << std::endl;
    }
}

void IOHandler::writeSettingsToFile(const std::string& settingsFilename) {
    // Open the file "settings.txt" for writing
    std::ofstream outputFile(settingsFilename);

    if (outputFile.is_open()) {
        outputFile << settingsMessage;
        outputFile.close();
        outputStream << "Created new \"" << Defaults::settingsFilename << "\" file.\n";
    } else {
        errorStream << "Error opening " 
            << Defaults::settingsFilename
            << " for writing.\n"
            << " Try again?";
        if (promptContinue()) {
            writeSettingsToFile(Defaults::settingsFilename);
        } else {
            throw ExitException(EXIT_SUCCESS);
        }
    }
}

bool IOHandler::parseSetting(const std::string& line, std::string& name, std::string& value) {
        std::istringstream iss(line);
    if (std::getline(iss, name, '=') && std::getline(iss, value)) {
        return true;
    }
    return false;
} 

bool IOHandler::stringToBool(const std::string& str) {
    return str == "True" || str == "true" || str == "1";
}

void IOHandler::handleFileOpenError(const std::string& settingsFilename) {

    DEBUG_PRINT("-----------------------Progress Checking 302-----------------------");

    if(input == nullptr) {
                DEBUG_PRINT("-----------------------Progress Checking 305-----------------------");
        throw std::logic_error(std::string(inputNotInitializedMessage));
    }

    DEBUG_PRINT("-----------------------Progress Checking 309-----------------------");

    // errorStream = std::cerr;

    DEBUG_PRINT("File that could not be opened: " << settingsFilename);
    errorStream << "Error: Unable to open file \""
                <<  settingsFilename << "\". ";
            // Check if the default settings file exists

    DEBUG_PRINT("-----------------------Progress Checking 315-----------------------");

    DEBUG_PRINT("Tryng to open file:" << Defaults::settingsFilename);
    std::ifstream defaultSettingsFile(Defaults::settingsFilename);
    DEBUG_PRINT("File opened:" << Defaults::settingsFilename);


    if (defaultSettingsFile.is_open()) {
                DEBUG_PRINT("-----------------------Progress Checking 317-----------------------");
        outputStream << "The file \""
                    << Defaults::settingsFilename
                    << "\" has been found. Continue with this file? Use -s for settings help" 
                    << std::endl;
    } else {
                DEBUG_PRINT("-----------------------Progress Checking 323-----------------------");
        errorStream << "Creating default settings file \""
                        << Defaults::settingsFilename
                        << "\". Use -s for settings help."
                        << std::endl;
    }

    if(promptContinue()) {
                DEBUG_PRINT("-----------------------Progress Checking 331-----------------------");
        arguments.filename = Defaults::settingsFilename;
        if(!defaultSettingsFile.is_open()) { 
                    DEBUG_PRINT("-----------------------Progress Checking 334-----------------------");
            writeSettingsToFile(arguments.filename);
        } else {
                    DEBUG_PRINT("-----------------------Progress Checking 337-----------------------");
            defaultSettingsFile.close();
        }
        return;
    } else {
        throw ExitException(EXIT_SUCCESS);
    }
}

void IOHandler::readSettingsFromFile() {

    DEBUG_PRINT("-----------------------Progress Checking 339-----------------------");

    if(input == nullptr) {
                DEBUG_PRINT("-----------------------Progress Checking 342-----------------------");
        throw std::logic_error(std::string(inputNotInitializedMessage));
    }

    DEBUG_PRINT("-----------------------Progress Checking 346-----------------------");

    std::ifstream file(arguments.filename);

    while (!file.is_open()) {
        DEBUG_PRINT("-----------------------Progress Checking 351-----------------------");
        handleFileOpenError(arguments.filename);
        file.open(arguments.filename);
        return;
    }
    
    std::string line;
    while (std::getline(file, line)) {
        DEBUG_PRINT("-----------------------Progress Checking 558-----------------------");
        DEBUG_PRINT("Line: " << line);
        std::string name, value;
        if (parseSetting(line, name, value)) {
            DEBUG_PRINT("Name: " << name << " Value: " << value);
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
            } else if (name == "logEndOfFlyingDistanceOnly") {
                input->setLogEndOfFlyingDistanceOnly(stringToBool(value));
            } else if (name == "logStoppingPointOnly") {
                input->setLogStoppingPointOnly(stringToBool(value));
            } else if (name == "simulationCount") {
                input->setSimulationCount(std::stoi(value));
            } else if (name == "numThreads") {
                input->setNumThreads(std::stoi(value));
            } else
                // Handle other `tings here if needed
                errorStream << "Warning: Unknown setting name '" << name << "'" << std::endl;
            }
        }
    file.close();
}

void IOHandler::handleNoOutputDirectory(const std::filesystem::path& outputPath) {
    const std::string ANSI_COLOR_RED = "\033[1;31m";
    const std::string ANSI_COLOR_RESET = "\033[0m";
    errorStream << "Unable to access output directory " << outputPath.string()
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

void IOHandler::handleOutputFileOpenError(
        OutputType outputType,
        std::filesystem::path outputPath,
        std::string filename,
        int bombardmentID) {

    if(input == nullptr) {
        throw std::logic_error(std::string(inputNotInitializedMessage));
    }

    if(outputType != ENDOFFILE) {
        errorStream << "Failed to open output file "
                    << (outputPath / filename).string() 
                    << ".\nData will be lost for simulation "
                    << bombardmentID
                    << "."
                    << std::endl;
        return;
    } else {
        errorStream << "Failed to open output file "
                    << (outputPath / filename).string()
                    << ".\nFile end marker \""
                    << input->getOutputFileEndMarker() << "\""
                    << " could not be appended." 
                    << std::endl;
    }
}

bool IOHandler::areArgumentsSet() const { return argumentsSet; };

#endif

// Electron screening potentials:

-184.39,185.39,0,2.0027,1.9973,0
-0.2259,1.2259,0,5.5272,2.3992,0
0.6045,0.3955,0,2.8174,0.6625,0
0.3278,0.6722,0,4.543,0.9852,0
0.2327,0.7673,0,5.99,1.2135,0
0.1537,0.8463,0,8.0404,1.4913,0
0.0996,0.9004,0,10.812,1.7687,0
0.0625,0.9375,0,14.823,2.0403,0
0.0368,0.9632,0,21.4,2.306,0
0.0188,0.9812,0,34.999,2.5662,0
0.7444,0.2556,0,4.1205,0.8718,0
0.6423,0.3577,0,4.7266,1.0025,0
0.6002,0.3998,0,5.1405,1.0153,0
0.516,0.484,0,5.8492,1.1732,0
0.4387,0.5613,0,6.6707,1.341,0
0.5459,-0.5333,0.9874,6.3703,2.5517,1.6753
0.7249,-0.7548,1.0299,6.2118,3.3883,1.8596
2.1912,-2.2852,1.094,5.547,4.5687,2.0446
0.0486,0.7759,0.1755,30.26,3.1243,0.7326
0.58,0.42,0,6.3218,1.0094,0
0.5543,0.4457,0,6.6328,1.1023,0
0.0112,0.6832,0.3056,99.757,4.1286,1.009
0.0318,0.6753,0.2929,42.533,3.9404,1.0533
0.1075,0.7162,0.1763,18.959,3.0638,1.0014
0.0498,0.6866,0.2636,31.864,3.7811,1.1279
0.0512,0.6995,0.2493,31.825,3.7716,1.1606
0.05,0.7142,0.2358,32.915,3.7908,1.1915
0.0474,0.7294,0.2232,34.758,3.8299,1.2209
0.0771,0.7951,0.1278,25.326,3.3928,1.1426
0.04,0.759,0.201,40.343,3.9465,1.2759
0.1083,0.7489,0.1428,20.192,3.4733,1.0064
0.061,0.7157,0.2233,29.2,4.1252,1.1845
0.0212,0.6709,0.3079,62.487,4.9502,1.3582
0.4836,0.5164,0,8.7824,1.6967,0
0.4504,0.5496,0,9.3348,1.79,0
0.419,0.581,0,9.9142,1.8835,0
0.1734,0.7253,0.1013,17.166,3.1103,0.7177
0.0336,0.7816,0.1848,55.208,4.2842,0.8578
0.0689,0.7202,0.2109,31.366,4.2412,0.9472
0.1176,0.6581,0.2243,22.054,4.0325,1.0181
0.2257,0.5821,0.1922,14.24,2.9702,1.017
0.2693,0.5763,0.1544,14.044,2.8611,1.0591
0.2201,0.5618,0.2181,15.918,3.3672,1.1548
0.2751,0.5943,0.1306,14.314,2.737,1.1092
0.2711,0.6119,0.117,14.654,2.7183,1.1234
0.2784,0.6067,0.1149,14.645,2.6155,1.4318
0.2562,0.6505,0.0933,15.588,2.7412,1.1408
0.2271,0.6155,0.1574,16.914,3.0841,1.2619
0.2492,0.644,0.1068,16.155,2.8819,0.9942
0.2153,0.6115,0.1732,17.793,3.2937,1.1478
0.1806,0.5767,0.2427,19.875,3.8092,1.2829
0.1308,0.5504,0.3188,24.154,4.6119,1.4195
0.0588,0.5482,0.393,39.996,5.9132,1.5471
0.4451,0.5549,0,11.805,1.7967,0
0.2708,0.6524,0.0768,16.591,2.6964,0.6814
0.1728,0.6845,0.1427,22.397,3.4595,0.8073
0.1947,0.6384,0.1669,20.764,3.4657,0.8911
0.1913,0.6467,0.162,21.235,3.4819,0.9011
0.1868,0.6558,0.1574,21.803,3.5098,0.9106
0.1665,0.7057,0.1278,23.949,3.5199,0.8486
0.1624,0.7133,0.1243,24.598,3.556,0.8569
0.158,0.721,0.121,25.297,3.5963,0.865
0.1538,0.7284,0.1178,26.017,3.6383,0.8731
0.1587,0.7024,0.1389,25.497,3.7364,0.955
0.1453,0.7426,0.1121,27.547,3.7288,0.889
0.1413,0.7494,0.1093,28.346,3.7763,0.8969
0.1374,0.7558,0.1068,29.16,3.8244,0.9048
0.1336,0.7619,0.1045,29.99,3.8734,0.9128
0.1299,0.768,0.1021,30.835,3.9233,0.9203
0.1267,0.7734,0.0999,31.681,3.9727,0.9288
0.1288,0.7528,0.1184,31.353,4.0904,1.0072
0.1303,0.7324,0.1373,31.217,4.2049,1.0946
0.1384,0.7096,0.152,30.077,4.2492,1.1697
0.15,0.6871,0.1629,28.63,4.2426,1.234
0.1608,0.6659,0.1733,27.568,4.2341,1.297
0.1722,0.6468,0.181,26.586,4.1999,1.3535
0.1834,0.6306,0.186,25.734,4.1462,1.4037
0.223,0.6176,0.1594,22.994,3.7346,1.4428
0.2289,0.6114,0.1597,22.864,3.6914,1.4886
0.2098,0.6004,0.1898,24.408,3.9643,1.5343
0.2708,0.6428,0.0864,20.941,3.2456,1.1121
0.238,0.6308,0.1312,22.987,3.6217,1.2373
0.2288,0.622,0.1492,23.792,3.7796,1.2534
0.1941,0.6105,0.1954,26.695,4.2582,1.3577
0.15,0.6031,0.2469,31.84,4.9285,1.4683
0.0955,0.606,0.2985,43.489,5.852,1.5736
0.3192,0.6233,0.0575,20.015,2.9091,0.7207
0.2404,0.6567,0.1029,24.501,3.5524,0.8376
0.2266,0.6422,0.1312,25.684,3.7922,0.9335
0.2176,0.624,0.1584,26.554,4.0044,1.0238
0.2413,0.6304,0.1283,25.193,3.678,0.9699
0.2448,0.6298,0.1254,25.252,3.6397,0.9825

# Mott scattering parameters:

1.00000E+00,2.67363E-08,7.11530E-08,-9.77030E-08,-6.69132E-07,-3.09263E-07
1.17182E-02,1.62222E-02,-5.90397E-05,-1.05585E-04,4.17873E-04,9.13843E-04
-2.65955E-01,-7.29531E-01,-4.99796E-01,2.83507E-04,-9.09042E-04,-2.20244E-03
-1.82348E-04,-8.86355E-05,-1.90554E-04,-2.49708E-04,6.35004E-04,1.73523E-03
4.70966E-05,-4.09705E-06,3.75218E-05,8.05645E-05,-1.90534E-04,-5.42847E-04
1.00000E+00,3.76476E-08,-3.05313E-07,-3.27422E-07,2.44235E-06,4.08754E-06
2.35767E-02,3.24642E-02,-6.37269E-04,-7.69160E-04,5.28004E-03,9.45642E-03
-2.73743E-01,-7.40767E-01,-4.98195E-01,1.74337E-03,-1.25798E-02,-2.24046E-02
-7.79128E-04,-4.14495E-04,-1.62657E-03,-1.37286E-03,1.04319E-02,1.83488E-02
2.02855E-04,1.94598E-06,4.30102E-04,4.32180E-04,-3.31526E-03,-5.81788E-03
1.00000E+00,7.00357E-08,-3.15076E-07,-4.24915E-07,2.45516E-06,4.90187E-06
3.55657E-02,4.87956E-02,-1.95525E-03,-2.78660E-03,1.65490E-02,3.11496E-02
-2.81171E-01,-7.52015E-01,-4.95329E-01,5.83548E-03,-3.39830E-02,-6.55379E-02
-1.83452E-03,-8.12746E-04,-3.84675E-03,-4.44467E-03,2.55871E-02,4.99483E-02
4.79031E-04,-3.89615E-05,1.01022E-03,1.39133E-03,-7.99398E-03,-1.56366E-02
1.00000E+00,7.58881E-08,4.70500E-08,2.48041E-07,-2.06053E-06,-1.97319E-06
4.76788E-02,6.52200E-02,-4.54331E-03,-6.50318E-03,3.76564E-02,7.17176E-02
-2.88203E-01,-7.63217E-01,-4.90337E-01,1.22839E-02,-6.86398E-02,-1.35769E-01
-3.37733E-03,-1.36514E-03,-7.51614E-03,-8.78592E-03,4.78572E-02,9.69021E-02
8.81822E-04,-1.02577E-04,1.99797E-03,2.72661E-03,-1.48296E-02,-3.01060E-02
9.99999E-01,7.91498E-08,1.84164E-06,2.68534E-06,-1.81630E-05,-2.69021E-05
5.98818E-02,8.17654E-02,-7.70811E-03,-1.12378E-02,6.38329E-02,1.25339E-01
-2.94716E-01,-7.74405E-01,-4.86220E-01,1.77367E-02,-9.46825E-02,-2.01789E-01
-5.52375E-03,-2.05348E-03,-9.44915E-03,-1.08135E-02,5.41024E-02,1.25257E-01
1.44555E-03,-1.99404E-04,2.36742E-03,3.29655E-03,-1.64122E-02,-3.83750E-02
9.99999E-01,7.68158E-08,5.18185E-06,7.34245E-06,-4.94780E-05,-7.71923E-05
7.21461E-02,9.84618E-02,-1.06535E-02,-1.62358E-02,8.59238E-02,1.78727E-01
-3.00622E-01,-7.85616E-01,-4.85735E-01,1.91563E-02,-8.10204E-02,-2.15386E-01
-8.34809E-03,-2.85241E-03,-7.03252E-03,-7.56786E-03,1.44975E-02,8.79093E-02
2.18964E-03,-3.42022E-04,1.32930E-03,2.20108E-03,-3.57927E-03,-2.59280E-02
9.99999E-01,8.36312E-08,1.09116E-05,1.47812E-05,-1.02733E-04,-1.62724E-04
8.44142E-02,1.15310E-01,-1.17230E-02,-1.94732E-02,8.92604E-02,2.09303E-01
-3.05743E-01,-7.96809E-01,-4.93957E-01,1.01607E-02,1.67761E-02,-1.05909E-01
-1.20090E-02,-3.80678E-03,4.51195E-03,6.93472E-03,-1.12405E-01,-8.15484E-02
3.16048E-03,-5.22237E-04,-2.58261E-03,-2.38303E-03,3.63393E-02,2.75127E-02
9.99998E-01,1.57323E-08,1.77595E-05,2.56082E-05,-1.67537E-04,-2.73755E-04
9.66438E-02,1.32264E-01,-9.53841E-03,-1.83707E-02,6.01664E-02,1.93357E-01
-3.09969E-01,-8.07790E-01,-5.14392E-01,-1.67153E-02,2.33870E-01,1.91600E-01
-1.65906E-02,-5.11585E-03,2.80424E-02,3.94663E-02,-3.55720E-01,-4.39251E-01
4.37866E-03,-6.81795E-04,-1.01520E-02,-1.24875E-02,1.11484E-01,1.38105E-01
9.99997E-01,-8.06132E-08,2.49797E-05,3.85120E-05,-2.37451E-04,-3.99607E-04
1.08782E-01,1.49306E-01,-2.50975E-03,-1.05471E-02,-1.64831E-02,1.05733E-01
-3.13165E-01,-8.18489E-01,-5.50832E-01,-6.74447E-02,6.06357E-01,7.39717E-01
-2.21976E-02,-6.84023E-03,6.65411E-02,9.48702E-02,-7.43989E-01,-1.03582E+00
5.87088E-03,-8.10480E-04,-2.21731E-02,-2.94422E-02,2.29540E-01,3.19669E-01
9.99997E-01,-1.87404E-07,3.10276E-05,5.20007E-05,-2.98132E-04,-5.19259E-04
1.20783E-01,1.66407E-01,1.06608E-02,6.48772E-03,-1.53031E-01,-7.59354E-02
-3.15222E-01,-8.28793E-01,-6.05740E-01,-1.47812E-01,1.15760E+00,1.58565E+00
-2.89055E-02,-9.08096E-03,1.21467E-01,1.77575E-01,-1.29110E+00,-1.90333E+00
7.65342E-03,-8.85417E-04,-3.89092E-02,-5.44040E-02,3.93087E-01,5.79439E-01
9.99996E-01,-2.44548E-07,3.31019E-05,6.29483E-05,-3.24667E-04,-5.95527E-04
1.32615E-01,1.83566E-01,3.04158E-02,3.40925E-02,-3.54681E-01,-3.63044E-01
-3.16092E-01,-8.38704E-01,-6.78558E-01,-2.59346E-01,1.88547E+00,2.73632E+00
-3.67233E-02,-1.18139E-02,1.91089E-01,2.87408E-01,-1.98397E+00,-3.03075E+00
9.72033E-03,-9.26380E-04,-5.95654E-02,-8.69829E-02,5.95744E-01,9.10242E-01
9.99995E-01,-2.12227E-07,2.95645E-05,6.92848E-05,-3.02153E-04,-6.05145E-04
1.44258E-01,2.00775E-01,5.67845E-02,7.35166E-02,-6.22861E-01,-7.62213E-01
-3.15754E-01,-8.48196E-01,-7.67318E-01,-4.02984E-01,2.77477E+00,4.18114E+00
-4.56307E-02,-1.50425E-02,2.72232E-01,4.23528E-01,-2.79606E+00,-4.38863E+00
1.20560E-02,-9.44637E-04,-8.28738E-02,-1.26564E-01,8.26726E-01,1.29882E+00
9.99995E-01,-4.03407E-08,1.86047E-05,6.85201E-05,-2.14503E-04,-5.22528E-04
1.55704E-01,2.18048E-01,8.88994E-02,1.24878E-01,-9.51331E-01,-1.26824E+00
-3.14244E-01,-8.57322E-01,-8.67190E-01,-5.75787E-01,3.78571E+00,5.87052E+00
-5.55526E-02,-1.86861E-02,3.58860E-01,5.81094E-01,-3.67623E+00,-5.90800E+00
1.46269E-02,-9.79742E-04,-1.06652E-01,-1.71226E-01,1.06737E+00,1.71918E+00
9.99994E-01,3.00267E-07,-1.11840E-06,5.88256E-05,-4.78456E-05,-3.25731E-04
1.66960E-01,2.35405E-01,1.25215E-01,1.87646E-01,-1.32685E+00,-1.86549E+00
-3.11630E-01,-8.66152E-01,-9.71254E-01,-7.72715E-01,4.85654E+00,7.72150E+00
-6.63778E-02,-2.26481E-02,4.42898E-01,7.53182E-01,-4.55172E+00,-7.48670E+00
1.73883E-02,-1.07669E-03,-1.28075E-01,-2.18389E-01,1.29217E+00,2.13475E+00
9.99994E-01,8.94829E-07,-2.98434E-05,3.82193E-05,2.00584E-04,-6.40482E-06
1.78039E-01,2.52912E-01,1.63761E-01,2.60132E-01,-1.73287E+00,-2.53185E+00
-3.08007E-01,-8.74905E-01,-1.07146E+00,-9.85062E-01,5.91697E+00,9.63265E+00
-7.79747E-02,-2.66797E-02,5.15288E-01,9.29261E-01,-5.34252E+00,-9.00574E+00
2.02892E-02,-1.33011E-03,-1.44039E-01,-2.64330E-01,1.47360E+00,2.50398E+00
9.99994E-01,1.75397E-06,-6.73310E-05,6.29524E-06,5.29623E-04,4.35288E-04
1.88968E-01,2.70612E-01,2.01975E-01,3.40574E-01,-2.14737E+00,-3.23836E+00
-3.03499E-01,-8.83717E-01,-1.15816E+00,-1.20414E+00,6.88176E+00,1.14841E+01
-9.01806E-02,-3.06202E-02,5.65581E-01,1.09902E+00,-5.95552E+00,-1.03302E+01
2.32694E-02,-1.80614E-03,-1.51041E-01,-3.05449E-01,1.58037E+00,2.78083E+00
9.99994E-01,3.07931E-06,-1.11876E-04,-4.10164E-05,9.17095E-04,9.80145E-04
1.99765E-01,2.88611E-01,2.37501E-01,4.25803E-01,-2.55105E+00,-3.95585E+00
-2.98206E-01,-8.92930E-01,-1.22279E+00,-1.41690E+00,7.67836E+00,1.31601E+01
-1.02865E-01,-3.40967E-02,5.84677E-01,1.24786E+00,-6.31301E+00,-1.13328E+01
2.62800E-02,-2.63995E-03,-1.46076E-01,-3.36795E-01,1.58677E+00,2.92251E+00
9.99993E-01,4.49776E-06,-1.65136E-04,-9.76754E-05,1.39664E-03,1.66293E-03
2.10469E-01,3.06924E-01,2.66793E-01,5.13797E-01,-2.90958E+00,-4.63816E+00
-2.92294E-01,-9.02560E-01,-1.25307E+00,-1.61470E+00,8.18574E+00,1.44912E+01
-1.15831E-01,-3.70891E-02,5.59807E-01,1.36619E+00,-6.28824E+00,-1.18327E+01
2.92513E-02,-3.84903E-03,-1.24976E-01,-3.55149E-01,1.45127E+00,2.86925E+00
9.99993E-01,6.01488E-06,-2.22125E-04,-1.61774E-04,1.92058E-03,2.41975E-03
2.21091E-01,3.25648E-01,2.87732E-01,6.01632E-01,-3.20436E+00,-5.25724E+00
-2.85814E-01,-9.12889E-01,-1.24213E+00,-1.78635E+00,8.34196E+00,1.53776E+01
-1.29005E-01,-3.92986E-02,4.84255E-01,1.44206E+00,-5.81999E+00,-1.17275E+01
3.21555E-02,-5.54272E-03,-8.56572E-02,-3.56589E-01,1.15480E+00,2.58829E+00
9.99993E-01,8.01467E-06,-2.79242E-04,-2.36820E-04,2.45459E-03,3.21683E-03
2.31651E-01,3.44948E-01,2.97820E-01,6.85187E-01,-3.41294E+00,-5.77715E+00
-2.78858E-01,-9.24428E-01,-1.18215E+00,-1.91691E+00,8.07489E+00,1.56969E+01
-1.42276E-01,-4.01888E-02,3.50466E-01,1.45983E+00,-4.83806E+00,-1.08936E+01
3.49529E-02,-7.90933E-03,-2.58002E-02,-3.36028E-01,6.75740E-01,2.04052E+00
9.99992E-01,1.04277E-05,-3.35126E-04,-3.21042E-04,2.98507E-03,4.03325E-03
2.42172E-01,3.64954E-01,2.94606E-01,7.60693E-01,-3.51409E+00,-6.16460E+00
-2.71512E-01,-9.37543E-01,-1.06570E+00,-1.99328E+00,7.31863E+00,1.53396E+01
-1.55540E-01,-3.93862E-02,1.51477E-01,1.40614E+00,-3.28024E+00,-9.22338E+00
3.76066E-02,-1.10812E-02,5.66831E-02,-2.89210E-01,-4.60274E-03,1.19273E+00
9.99992E-01,1.30838E-05,-3.84070E-04,-4.09294E-04,3.47025E-03,4.81071E-03
2.52646E-01,3.85718E-01,2.76330E-01,8.25665E-01,-3.48888E+00,-6.39017E+00
-2.63758E-01,-9.52326E-01,-8.88179E-01,-2.00720E+00,6.01960E+00,1.42157E+01
-1.68806E-01,-3.68095E-02,-1.16429E-01,1.27309E+00,-1.09991E+00,-6.63409E+00
4.01184E-02,-1.50938E-02,1.62740E-01,-2.13760E-01,-8.99142E-01,2.08129E-02
9.99991E-01,1.59363E-05,-4.27315E-04,-5.01461E-04,3.91365E-03,5.55096E-03
2.63096E-01,4.07357E-01,2.40649E-01,8.76450E-01,-3.31641E+00,-6.42069E+00
-2.55682E-01,-9.69090E-01,-6.43149E-01,-1.94678E+00,4.11915E+00,1.22255E+01
-1.81974E-01,-3.21462E-02,-4.58817E-01,1.04913E+00,1.75388E+00,-3.03434E+00
4.24541E-02,-2.00595E-02,2.93916E-01,-1.06138E-01,-2.02203E+00,-1.50181E+00
9.99990E-01,1.88950E-05,-4.59994E-04,-5.93663E-04,4.27684E-03,6.20090E-03
2.73504E-01,4.29990E-01,1.86473E-01,9.09921E-01,-2.98441E+00,-6.23344E+00
-2.47225E-01,-9.88118E-01,-3.28759E-01,-1.80252E+00,1.59090E+00,9.30968E+00
-1.95100E-01,-2.51242E-02,-8.76244E-01,7.25592E-01,5.29650E+00,1.62177E+00
4.46307E-02,-2.60754E-02,4.50054E-01,3.61638E-02,-3.37524E+00,-3.38601E+00
9.99989E-01,2.18906E-05,-4.79199E-04,-6.83164E-04,4.53467E-03,6.72485E-03
2.83854E-01,4.53722E-01,1.12877E-01,9.23203E-01,-2.48247E+00,-5.80855E+00
-2.38338E-01,-1.00965E+00,5.61453E-02,-1.56605E+00,-1.58353E+00,5.42138E+00
-2.08226E-01,-1.55138E-02,-1.36846E+00,2.95168E-01,9.53392E+00,7.36550E+00
4.66619E-02,-3.32255E-02,6.30711E-01,2.15167E-01,-4.95748E+00,-5.63744E+00
9.99987E-01,2.48200E-05,-4.82895E-04,-7.67488E-04,4.66900E-03,7.09581E-03
2.94123E-01,4.78653E-01,1.93256E-02,9.13569E-01,-1.80354E+00,-5.13040E+00
-2.28945E-01,-1.03390E+00,5.11210E-01,-1.22996E+00,-5.40939E+00,5.31490E-01
-2.21426E-01,-3.12679E-03,-1.93353E+00,-2.48179E-01,1.44576E+01,1.42077E+01
4.85718E-02,-4.15796E-02,8.34879E-01,4.32418E-01,-6.76244E+00,-8.25464E+00
9.99985E-01,2.76168E-05,-4.67522E-04,-8.43978E-04,4.65008E-03,7.27476E-03
3.04280E-01,5.04856E-01,-9.42256E-02,8.78905E-01,-9.44256E-01,-4.18847E+00
-2.18939E-01,-1.06098E+00,1.03431E+00,-7.89123E-01,-9.87815E+00,-5.36997E+00
-2.34805E-01,1.21281E-02,-2.56765E+00,-9.08015E-01,2.00436E+01,2.21379E+01
5.03942E-02,-5.11768E-02,1.06090E+00,6.88631E-01,-8.77867E+00,-1.12289E+01
9.99982E-01,3.00792E-05,-4.33447E-04,-9.09366E-04,4.47786E-03,7.25072E-03
3.14283E-01,5.32444E-01,-2.27528E-01,8.16951E-01,9.53704E-02,-2.97730E+00
-2.08186E-01,-1.09112E+00,1.62237E+00,-2.38394E-01,-1.49699E+01,-1.22743E+01
-2.48493E-01,3.04543E-02,-3.26601E+00,-1.68760E+00,2.62568E+01,3.11259E+01
5.21709E-02,-6.20908E-02,1.30683E+00,9.84353E-01,-1.09909E+01,-1.45449E+01
9.99979E-01,3.24569E-05,-3.76717E-04,-9.64000E-04,4.11997E-03,6.98761E-03
3.24082E-01,5.61534E-01,-3.79993E-01,7.25313E-01,1.31342E+00,-1.49163E+00
-1.96521E-01,-1.12457E+00,2.27091E+00,4.27678E-01,-2.06558E+01,-2.01690E+01
-2.62655E-01,5.20812E-02,-4.02224E+00,-2.59048E+00,3.30508E+01,4.11350E+01
5.39556E-02,-7.44077E-02,1.57015E+00,1.32022E+00,-1.33802E+01,-1.81848E+01
9.99976E-01,3.39628E-05,-2.98845E-04,-9.99651E-04,3.58523E-03,6.47782E-03
3.33647E-01,5.92068E-01,-5.51020E-01,6.03630E-01,2.70712E+00,2.67853E-01
-1.83849E-01,-1.16095E+00,2.97559E+00,1.20728E+00,-2.69053E+01,-2.90214E+01
-2.77375E-01,7.65838E-02,-4.83023E+00,-3.61238E+00,4.03786E+01,5.21084E+01
5.57745E-02,-8.79952E-02,1.84848E+00,1.69425E+00,-1.59274E+01,-2.21242E+01
9.99972E-01,3.48473E-05,-1.97828E-04,-1.01659E-03,2.85509E-03,5.69725E-03
3.42904E-01,6.24164E-01,-7.39023E-01,4.50306E-01,4.26553E+00,2.29299E+00
-1.69940E-01,-1.20051E+00,3.72868E+00,2.10276E+00,-3.36594E+01,-3.87726E+01
-2.92882E-01,1.04194E-01,-5.68033E+00,-4.75339E+00,4.81633E+01,6.39614E+01
5.77008E-02,-1.02941E-01,2.13826E+00,2.10589E+00,-1.86034E+01,-2.63294E+01
9.99968E-01,3.47804E-05,-7.11898E-05,-1.01028E-03,1.90919E-03,4.61426E-03
3.51793E-01,6.57815E-01,-9.42175E-01,2.64940E-01,5.97606E+00,4.57222E+00
-1.54595E-01,-1.24305E+00,4.52170E+00,3.11205E+00,-4.08539E+01,-4.93518E+01
-3.09367E-01,1.34661E-01,-6.56214E+00,-6.00882E+00,5.63228E+01,7.65965E+01
5.97948E-02,-1.19174E-01,2.43570E+00,2.55300E+00,-2.13779E+01,-3.07628E+01
9.99963E-01,3.37519E-05,8.13037E-05,-9.78638E-04,7.41412E-04,3.21498E-03
3.60246E-01,6.93065E-01,-1.15880E+00,4.68519E-02,7.82662E+00,7.09456E+00
-1.37615E-01,-1.28855E+00,5.34673E+00,4.23402E+00,-4.84263E+01,-6.06897E+01
-3.27017E-01,1.67941E-01,-7.46593E+00,-7.37492E+00,6.47773E+01,8.99181E+01
6.21158E-02,-1.36692E-01,2.73726E+00,3.03374E+00,-2.42209E+01,-3.53873E+01
9.99958E-01,3.14983E-05,2.59741E-04,-9.19008E-04,-6.49202E-04,1.49008E-03
3.68196E-01,7.29888E-01,-1.38664E+00,-2.03807E-01,9.80062E+00,9.84155E+00
-1.18788E-01,-1.33674E+00,6.19400E+00,5.46446E+00,-5.62994E+01,-7.26915E+01
-3.46033E-01,2.03729E-01,-8.38012E+00,-8.84466E+00,7.34319E+01,1.03805E+02
6.47255E-02,-1.55409E-01,3.03879E+00,3.54516E+00,-2.70980E+01,-4.01572E+01
9.99952E-01,2.79961E-05,4.60479E-04,-8.33486E-04,-2.23214E-03,-5.16285E-04
3.75580E-01,7.68292E-01,-1.62392E+00,-4.87674E-01,1.18854E+01,1.28025E+01
-9.79355E-02,-1.38749E+00,7.05560E+00,6.80198E+00,-6.44113E+01,-8.52915E+01
-3.66572E-01,2.41863E-01,-9.29524E+00,-1.04141E+01,8.22093E+01,1.18166E+02
6.76714E-02,-1.75289E-01,3.33691E+00,4.08538E+00,-2.99809E+01,-4.50381E+01
9.99947E-01,2.29520E-05,6.82639E-04,-7.17139E-04,-4.00522E-03,-2.81601E-03
3.82332E-01,8.08194E-01,-1.86848E+00,-8.03415E-01,1.40640E+01,1.59560E+01
-7.48735E-02,-1.44031E+00,7.92233E+00,8.23820E+00,-7.26858E+01,-9.83865E+01
-3.88797E-01,2.81837E-01,-1.02005E+01,-1.20716E+01,9.10175E+01,1.32873E+02
7.10017E-02,-1.96182E-01,3.62780E+00,4.65000E+00,-3.28365E+01,-4.99823E+01
9.99941E-01,1.63607E-05,9.28242E-04,-5.69257E-04,-5.98245E-03,-5.42476E-03
3.88363E-01,8.49548E-01,-2.11701E+00,-1.15003E+00,1.63119E+01,1.92737E+01
-4.93364E-02,-1.49491E+00,8.78129E+00,9.76596E+00,-8.10212E+01,-1.11851E+02
-4.12953E-01,3.23322E-01,-1.10814E+01,-1.38073E+01,9.97386E+01,1.47775E+02
7.47911E-02,-2.18001E-01,3.90642E+00,5.23509E+00,-3.56231E+01,-5.49352E+01
9.99935E-01,8.31520E-06,1.19244E-03,-3.95258E-04,-8.12525E-03,-8.28836E-03
3.93603E-01,8.92358E-01,-2.36688E+00,-1.52786E+00,1.86071E+01,2.27314E+01
-2.11345E-02,-1.55112E+00,9.62205E+00,1.13827E+01,-8.93266E+01,-1.25574E+02
-4.39199E-01,3.66161E-01,-1.19262E+01,-1.56159E+01,1.08267E+02,1.62736E+02
7.90853E-02,-2.40719E-01,4.16872E+00,5.83839E+00,-3.83040E+01,-5.98484E+01
9.99929E-01,-1.65608E-06,1.47476E-03,-1.86250E-04,-1.04363E-02,-1.14318E-02
3.97989E-01,9.36430E-01,-2.61569E+00,-1.93387E+00,2.09305E+01,2.63021E+01
9.88744E-03,-1.60813E+00,1.04351E+01,1.30740E+01,-9.75207E+01,-1.39437E+02
-4.67657E-01,4.09494E-01,-1.27240E+01,-1.74797E+01,1.16508E+02,1.77614E+02
8.39169E-02,-2.64077E-01,4.41092E+00,6.45344E+00,-4.08455E+01,-6.46698E+01
9.99922E-01,-1.37624E-05,1.77335E-03,5.99400E-05,-1.29013E-02,-1.48443E-02
4.01450E-01,9.81600E-01,-2.86053E+00,-2.36562E+00,2.32591E+01,2.99551E+01
4.39150E-02,-1.66522E+00,1.12093E+01,1.48279E+01,-1.05512E+02,-1.53310E+02
-4.98475E-01,4.52607E-01,-1.34628E+01,-1.93836E+01,1.24357E+02,1.92256E+02
8.93258E-02,-2.87869E-01,4.62891E+00,7.07473E+00,-4.32116E+01,-6.93461E+01
9.99916E-01,-2.78975E-05,2.08753E-03,3.41607E-04,-1.55141E-02,-1.85156E-02
4.03889E-01,1.02783E+00,-3.09762E+00,-2.82183E+00,2.55629E+01,3.36548E+01
8.12122E-02,-1.72210E+00,1.19310E+01,1.66362E+01,-1.13185E+02,-1.67048E+02
-5.31880E-01,4.95236E-01,-1.41275E+01,-2.13166E+01,1.31687E+02,2.06495E+02
9.53781E-02,-3.12041E-01,4.81765E+00,7.69813E+00,-4.53587E+01,-7.38182E+01
9.99910E-01,-4.39880E-05,2.41439E-03,6.54636E-04,-1.82489E-02,-2.24062E-02
4.05249E-01,1.07495E+00,-3.32445E+00,-3.30073E+00,2.78221E+01,3.73760E+01
1.21896E-01,-1.77813E+00,1.25907E+01,1.84893E+01,-1.20461E+02,-1.80542E+02
-5.67942E-01,5.36741E-01,-1.47080E+01,-2.32664E+01,1.38408E+02,2.20203E+02
1.02086E-01,-3.36419E-01,4.97373E+00,8.31905E+00,-4.72562E+01,-7.80409E+01
9.99904E-01,-6.22073E-05,2.75066E-03,1.00063E-03,-2.10817E-02,-2.64933E-02
4.05472E-01,1.12279E+00,-3.53825E+00,-3.79964E+00,3.00132E+01,4.10871E+01
1.66089E-01,-1.83255E+00,1.31785E+01,2.03743E+01,-1.27252E+02,-1.93664E+02
-6.06729E-01,5.76392E-01,-1.51936E+01,-2.52171E+01,1.44423E+02,2.33232E+02
1.09463E-01,-3.60803E-01,5.09358E+00,8.93178E+00,-4.88706E+01,-8.19632E+01
9.99898E-01,-8.26232E-05,3.09353E-03,1.37924E-03,-2.39902E-02,-3.07502E-02
4.04498E-01,1.17115E+00,-3.73628E+00,-4.31580E+00,3.21142E+01,4.47580E+01
2.13904E-01,-1.88458E+00,1.36848E+01,2.22782E+01,-1.33471E+02,-2.06293E+02
-6.48301E-01,6.13439E-01,-1.55742E+01,-2.71531E+01,1.49638E+02,2.45444E+02
1.17516E-01,-3.84990E-01,5.17383E+00,9.53070E+00,-5.01709E+01,-8.55368E+01
9.99893E-01,-1.05293E-04,3.43932E-03,1.78981E-03,-2.69462E-02,-3.51436E-02
4.02272E-01,1.21982E+00,-3.91575E+00,-4.84621E+00,3.41013E+01,4.83563E+01
2.65438E-01,-1.93341E+00,1.40998E+01,2.41876E+01,-1.39034E+02,-2.18300E+02
-6.92691E-01,6.47135E-01,-1.58398E+01,-2.90581E+01,1.53960E+02,2.56697E+02
1.26243E-01,-4.08783E-01,5.21125E+00,1.01100E+01,-5.11262E+01,-8.87130E+01
9.99888E-01,-1.30188E-04,3.78436E-03,2.22971E-03,-2.99185E-02,-3.96315E-02
3.98750E-01,1.26855E+00,-4.07427E+00,-5.38796E+00,3.59550E+01,5.18546E+01
3.20748E-01,-1.97817E+00,1.44153E+01,2.60890E+01,-1.43866E+02,-2.29578E+02
-7.39891E-01,6.76679E-01,-1.59819E+01,-3.09162E+01,1.57312E+02,2.66867E+02
1.35630E-01,-4.31972E-01,5.20310E+00,1.06642E+01,-5.17104E+01,-9.14493E+01
9.99884E-01,-1.57459E-04,4.12052E-03,2.69984E-03,-3.28468E-02,-4.41512E-02
3.93891E-01,1.31710E+00,-4.20944E+00,-5.93691E+00,3.76535E+01,5.52188E+01
3.79857E-01,-2.01791E+00,1.46238E+01,2.79653E+01,-1.47893E+02,-2.39997E+02
-7.89852E-01,7.01196E-01,-1.59929E+01,-3.27077E+01,1.59613E+02,2.75813E+02
1.45646E-01,-4.54328E-01,5.14701E+00,1.11863E+01,-5.18977E+01,-9.36987E+01
9.99881E-01,-1.86816E-04,4.45491E-03,3.19574E-03,-3.57812E-02,-4.87457E-02
3.87612E-01,1.36524E+00,-4.31686E+00,-6.49005E+00,3.91612E+01,5.84032E+01
4.42940E-01,-2.05189E+00,1.47103E+01,2.98033E+01,-1.50988E+02,-2.49389E+02
-8.42683E-01,7.20045E-01,-1.58579E+01,-3.44169E+01,1.60733E+02,2.83354E+02
1.56312E-01,-4.75705E-01,5.03810E+00,1.16708E+01,-5.16454E+01,-9.53999E+01
9.99878E-01,-2.18480E-04,4.78376E-03,3.71897E-03,-3.86926E-02,-5.33866E-02
3.79879E-01,1.41265E+00,-4.39391E+00,-7.04282E+00,4.04550E+01,6.13721E+01
5.10001E-01,-2.07896E+00,1.46666E+01,3.15843E+01,-1.53072E+02,-2.57624E+02
-8.98306E-01,7.32213E-01,-1.55686E+01,-3.60229E+01,1.60591E+02,2.89348E+02
1.67591E-01,-4.95837E-01,4.87389E+00,1.21106E+01,-5.09274E+01,-9.65060E+01
9.99876E-01,-2.52173E-04,5.09558E-03,4.26518E-03,-4.14939E-02,-5.79712E-02
3.70689E-01,1.45908E+00,-4.43983E+00,-7.59192E+00,4.15261E+01,6.41080E+01
5.80924E-01,-2.09831E+00,1.44911E+01,3.32945E+01,-1.54118E+02,-2.64634E+02
-9.56518E-01,7.37009E-01,-1.51243E+01,-3.75098E+01,1.59159E+02,2.93723E+02
1.79397E-01,-5.14573E-01,4.65426E+00,1.25003E+01,-4.97362E+01,-9.69932E+01
9.99874E-01,-2.87917E-04,5.39066E-03,4.83313E-03,-4.41852E-02,-6.24953E-02
3.60001E-01,1.50420E+00,-4.45150E+00,-8.13270E+00,4.23476E+01,6.65709E+01
6.55716E-01,-2.10882E+00,1.41741E+01,3.49157E+01,-1.54034E+02,-2.70276E+02
-1.01724E+00,7.33501E-01,-1.45156E+01,-3.88572E+01,1.56348E+02,2.96328E+02
1.91690E-01,-5.31690E-01,4.37639E+00,1.28328E+01,-4.80435E+01,-9.68125E+01
9.99874E-01,-3.25847E-04,5.67083E-03,5.42945E-03,-4.67877E-02,-6.70100E-02
3.47783E-01,1.54767E+00,-4.42606E+00,-8.65995E+00,4.28934E+01,6.87183E+01
7.34353E-01,-2.10941E+00,1.37070E+01,3.64273E+01,-1.52736E+02,-2.74401E+02
-1.08036E+00,7.20771E-01,-1.37344E+01,-4.00427E+01,1.52073E+02,2.97007E+02
2.04417E-01,-5.46972E-01,4.03784E+00,1.31007E+01,-4.58224E+01,-9.59123E+01
9.99875E-01,-3.65406E-04,5.92115E-03,6.03332E-03,-4.91674E-02,-7.12951E-02
3.34062E-01,1.58924E+00,-4.36361E+00,-9.17241E+00,4.31638E+01,7.05522E+01
8.16604E-01,-2.09933E+00,1.30915E+01,3.78235E+01,-1.50231E+02,-2.77015E+02
-1.14554E+00,6.98298E-01,-1.27840E+01,-4.10591E+01,1.46348E+02,2.95770E+02
2.17451E-01,-5.60344E-01,3.63994E+00,1.33015E+01,-4.30802E+01,-9.42981E+01
9.99877E-01,-4.06637E-04,6.14278E-03,6.64937E-03,-5.13391E-02,-7.53881E-02
3.18822E-01,1.62854E+00,-4.26169E+00,-9.66459E+00,4.31357E+01,7.20331E+01
9.02373E-01,-2.07743E+00,1.23209E+01,3.90832E+01,-1.46445E+02,-2.77984E+02
-1.21259E+00,6.65174E-01,-1.16583E+01,-4.18839E+01,1.39104E+02,2.92479E+02
2.30710E-01,-5.71608E-01,3.18104E+00,1.34277E+01,-3.97963E+01,-9.19255E+01
9.99881E-01,-4.49175E-04,6.32631E-03,7.27184E-03,-5.32297E-02,-7.91972E-02
3.02086E-01,1.66529E+00,-4.11996E+00,-1.01330E+01,4.28031E+01,7.31464E+01
9.91429E-01,-2.04292E+00,1.13959E+01,4.01931E+01,-1.41369E+02,-2.77260E+02
-1.28117E+00,6.20870E-01,-1.03596E+01,-4.25024E+01,1.30338E+02,2.87090E+02
2.44065E-01,-5.80702E-01,2.66217E+00,1.34745E+01,-3.59721E+01,-8.87824E+01
9.99886E-01,-4.93136E-04,6.47914E-03,7.90029E-03,-5.48909E-02,-8.27736E-02
2.83840E-01,1.69905E+00,-3.93534E+00,-1.05724E+01,4.21392E+01,7.38521E+01
1.08366E+00,-1.99451E+00,1.03076E+01,4.11340E+01,-1.34917E+02,-2.74710E+02
-1.35107E+00,5.64397E-01,-8.87993E+00,-4.28945E+01,1.19972E+02,2.79470E+02
2.57431E-01,-5.87417E-01,2.08111E+00,1.34353E+01,-3.15844E+01,-8.48268E+01
9.99892E-01,-5.38122E-04,6.59206E-03,8.52848E-03,-5.62508E-02,-8.60251E-02
2.64130E-01,1.72952E+00,-3.70812E+00,-1.09794E+01,4.11430E+01,7.41410E+01
1.17876E+00,-1.93141E+00,9.05928E+00,4.18928E+01,-1.27099E+02,-2.70309E+02
-1.42185E+00,4.95292E-01,-7.22427E+00,-4.30461E+01,1.08024E+02,2.69601E+02
2.70647E-01,-5.91727E-01,1.43982E+00,1.33056E+01,-2.66423E+01,-8.00557E+01
9.99900E-01,-5.83896E-04,6.66373E-03,9.15389E-03,-5.72988E-02,-8.89368E-02
2.42961E-01,1.75635E+00,-3.43619E+00,-1.13482E+01,3.97928E+01,7.39747E+01
1.27654E+00,-1.85263E+00,7.64589E+00,4.24492E+01,-1.17853E+02,-2.63933E+02
-1.49322E+00,4.12926E-01,-5.38909E+00,-4.29363E+01,9.44419E+01,2.57363E+02
2.83597E-01,-5.93557E-01,7.37540E-01,1.30785E+01,-2.11316E+01,-7.44329E+01
9.99909E-01,-6.30100E-04,6.69295E-03,9.77077E-03,-5.80220E-02,-9.14799E-02
2.20373E-01,1.77919E+00,-3.11907E+00,-1.16754E+01,3.80818E+01,7.33384E+01
1.37668E+00,-1.75734E+00,6.06819E+00,4.27905E+01,-1.07169E+02,-2.55543E+02
-1.56475E+00,3.16885E-01,-3.37687E+00,-4.25518E+01,7.92269E+01,2.42726E+02
2.96124E-01,-5.92915E-01,-2.45072E-02,1.27501E+01,-1.50562E+01,-6.79514E+01
9.99920E-01,-6.77211E-04,6.67917E-03,1.03866E-02,-5.84208E-02,-9.36775E-02
1.96406E-01,1.79753E+00,-2.75530E+00,-1.19533E+01,3.59929E+01,7.21943E+01
1.47887E+00,-1.64416E+00,4.32381E+00,4.28906E+01,-9.50034E+01,-2.45021E+02
-1.63600E+00,2.06169E-01,-1.18704E+00,-4.18661E+01,6.23476E+01,2.25579E+02
3.08067E-01,-5.89622E-01,-8.46080E-01,1.23116E+01,-8.40939E+00,-6.05788E+01
9.99932E-01,-7.23543E-04,6.61026E-03,1.09756E-02,-5.83892E-02,-9.53462E-02
1.71135E-01,1.81127E+00,-2.34665E+00,-1.21818E+01,3.35368E+01,7.05547E+01
1.58269E+00,-1.51310E+00,2.42140E+00,4.27487E+01,-8.14091E+01,-2.32422E+02
-1.70642E+00,8.13302E-02,1.16945E+00,-4.08790E+01,4.38728E+01,2.05995E+02
3.19226E-01,-5.84023E-01,-1.72313E+00,1.17636E+01,-1.21776E+00,-5.23434E+01
9.99946E-01,-7.69368E-04,6.48923E-03,1.15425E-02,-5.79523E-02,-9.65284E-02
1.44613E-01,1.81998E+00,-1.89196E+00,-1.23548E+01,3.06998E+01,6.83894E+01
1.68776E+00,-1.36316E+00,3.59801E-01,4.23443E+01,-6.63559E+01,-2.17658E+02
-1.77552E+00,-5.81387E-02,3.69184E+00,-3.95701E+01,2.37855E+01,1.83897E+02
3.29420E-01,-5.76119E-01,-2.65496E+00,1.10995E+01,6.52027E+00,-4.32241E+01
9.99962E-01,-8.13654E-04,6.31217E-03,1.20763E-02,-5.70782E-02,-9.71683E-02
1.16915E-01,1.82343E+00,-1.39208E+00,-1.24689E+01,2.74833E+01,6.56914E+01
1.79364E+00,-1.19400E+00,-1.85513E+00,4.16665E+01,-4.98680E+01,-2.00724E+02
-1.84272E+00,-2.12012E-01,6.37189E+00,-3.79293E+01,2.12742E+00,1.59295E+02
3.38449E-01,-5.66167E-01,-3.63838E+00,1.03165E+01,1.47866E+01,-3.32296E+01
9.99980E-01,-8.57331E-04,6.07736E-03,1.25819E-02,-5.57471E-02,-9.72475E-02
8.81540E-02,1.82105E+00,-8.47718E-01,-1.25186E+01,2.38891E+01,6.24510E+01
1.89976E+00,-1.00428E+00,-4.21799E+00,4.06970E+01,-3.19703E+01,-1.81604E+02
-1.90731E+00,-3.81080E-01,9.20190E+00,-3.59386E+01,-2.10594E+01,1.32192E+02
3.46059E-01,-5.54093E-01,-4.67039E+00,9.40924E+00,2.35631E+01,-2.23653E+01
9.99999E-01,-8.99343E-04,5.78680E-03,1.30512E-02,-5.39764E-02,-9.67723E-02
5.84024E-02,1.81266E+00,-2.59193E-01,-1.25002E+01,1.99139E+01,5.86545E+01
2.00567E+00,-7.93885E-01,-6.72461E+00,3.94244E+01,-1.26714E+01,-1.60273E+02
-1.96873E+00,-5.64782E-01,1.21753E+01,-3.35880E+01,-4.57487E+01,1.02578E+02
3.52050E-01,-5.40282E-01,-5.74837E+00,8.37487E+00,3.28364E+01,-1.06341E+01
1.00002E+00,-9.39311E-04,5.43607E-03,1.34763E-02,-5.17249E-02,-9.56765E-02
2.78005E-02,1.79795E+00,3.70737E-01,-1.24114E+01,1.55759E+01,5.43184E+01
2.11070E+00,-5.62473E-01,-9.36239E+00,3.78424E+01,7.94516E+00,-1.36806E+02
-2.02617E+00,-7.62769E-01,1.52768E+01,-3.08728E+01,-7.18354E+01,7.05554E+01
3.56140E-01,-5.25061E-01,-6.86679E+00,7.21247E+00,4.25672E+01,1.92463E+00
1.00004E+00,-9.76747E-04,5.02266E-03,1.38530E-02,-4.89737E-02,-9.39340E-02
-3.54744E-03,1.77664E+00,1.04132E+00,-1.22473E+01,1.08741E+01,4.94269E+01
2.21429E+00,-3.09743E-01,-1.21254E+01,3.59359E+01,2.98600E+01,-1.11177E+02
-2.07893E+00,-9.74603E-01,1.84979E+01,-2.77795E+01,-9.92800E+01,3.61173E+01
3.58090E-01,-5.08796E-01,-8.02233E+00,5.91833E+00,5.27374E+01,1.53069E+01
1.00007E+00,-1.01117E-03,4.54349E-03,1.41720E-02,-4.56907E-02,-9.14834E-02
-3.54805E-02,1.74848E+00,1.74921E+00,-1.20063E+01,5.83106E+00,4.40029E+01
2.31570E+00,-3.56472E-02,-1.49988E+01,3.37013E+01,5.29727E+01,-8.34879E+01
-2.12613E+00,-1.19953E+00,2.18208E+01,-2.43067E+01,-1.27958E+02,-6.05797E-01
3.57593E-01,-4.91966E-01,-9.20870E+00,4.49274E+00,6.33009E+01,2.94631E+01
1.00009E+00,-1.04219E-03,3.99979E-03,1.44301E-02,-4.18849E-02,-8.83308E-02
-6.78832E-02,1.71323E+00,2.49363E+00,-1.16834E+01,4.45978E-01,3.80313E+01
2.41433E+00,2.59886E-01,-1.79765E+01,3.11249E+01,7.72617E+01,-5.37152E+01
-2.16704E+00,-1.43680E+00,2.52369E+01,-2.04429E+01,-1.57828E+02,-3.96140E+01
3.54400E-01,-4.75067E-01,-1.04225E+01,2.93269E+00,7.42385E+01,4.43866E+01
1.00012E+00,-1.06927E-03,3.39091E-03,1.46208E-02,-3.75471E-02,-8.44506E-02
-1.00611E-01,1.67068E+00,3.27246E+00,-1.12756E+01,-5.27058E+00,3.15152E+01
2.50948E+00,5.76648E-01,-2.10477E+01,2.81988E+01,1.02665E+02,-2.18989E+01
-2.20084E+00,-1.68528E+00,2.87325E+01,-1.61828E+01,-1.88802E+02,-8.08394E+01
3.48226E-01,-4.58722E-01,-1.16589E+01,1.23735E+00,8.55166E+01,6.00477E+01
1.00015E+00,-1.09216E-03,2.71745E-03,1.47404E-02,-3.26779E-02,-7.98343E-02
-1.33495E-01,1.62059E+00,4.08307E+00,-1.07793E+01,-1.13042E+01,2.44607E+01
2.60038E+00,9.14519E-01,-2.42000E+01,2.49142E+01,1.29108E+02,1.19101E+01
-2.22661E+00,-1.94391E+00,3.22924E+01,-1.15206E+01,-2.20784E+02,-1.24202E+02
3.38764E-01,-4.43546E-01,-1.29122E+01,-5.94280E-01,9.70972E+01,7.64128E+01
1.00018E+00,-1.11005E-03,1.98019E-03,1.47814E-02,-2.72814E-02,-7.44728E-02
-1.66392E-01,1.56288E+00,4.92359E+00,-1.01915E+01,-1.76472E+01,1.68656E+01
2.68631E+00,1.27279E+00,-2.74231E+01,2.12642E+01,1.56537E+02,4.76851E+01
-2.24352E+00,-2.21096E+00,3.59035E+01,-6.45262E+00,-2.53694E+02,-1.69646E+02
3.25734E-01,-4.30391E-01,-1.41779E+01,-2.56236E+00,1.08949E+02,9.34555E+01
1.00021E+00,-1.12270E-03,1.18750E-03,1.47467E-02,-2.14258E-02,-6.84575E-02
-1.99145E-01,1.49739E+00,5.79226E+00,-9.50812E+00,-2.42938E+01,8.72312E+00
2.76656E+00,1.65088E+00,-3.07077E+01,1.72391E+01,1.84905E+02,8.54178E+01
-2.25073E+00,-2.48478E+00,3.95541E+01,-9.72452E-01,-2.87462E+02,-2.17133E+02
3.08850E-01,-4.20088E-01,-1.54515E+01,-4.66797E+00,1.21043E+02,1.11156E+02
1.00025E+00,-1.12974E-03,3.25083E-04,1.46196E-02,-1.49818E-02,-6.15916E-02
-2.31499E-01,1.42403E+00,6.68301E+00,-8.72894E+00,-3.12005E+01,8.13039E-02
2.84004E+00,2.04799E+00,-3.40292E+01,1.28426E+01,2.14035E+02,1.24907E+02
-2.24701E+00,-2.76343E+00,4.32160E+01,4.91177E+00,-3.21880E+02,-2.66423E+02
2.87708E-01,-4.13579E-01,-1.67235E+01,-6.90724E+00,1.33305E+02,1.29426E+02
1.00028E+00,-1.13018E-03,-5.85208E-04,1.44008E-02,-8.12721E-03,-5.40793E-02
-2.63328E-01,1.34285E+00,7.59541E+00,-7.85118E+00,-3.83723E+01,-9.07534E+00
2.90612E+00,2.46267E+00,-3.73829E+01,8.06978E+00,2.43918E+02,1.66172E+02
-2.23164E+00,-3.04426E+00,4.68824E+01,1.12007E+01,-3.56918E+02,-3.17506E+02
2.62068E-01,-4.12056E-01,-1.79911E+01,-9.27922E+00,1.45719E+02,1.48255E+02
1.00032E+00,-1.12442E-03,-1.55891E-03,1.40881E-02,-7.32260E-04,-4.57754E-02
-2.94323E-01,1.25374E+00,8.52167E+00,-6.87354E+00,-4.57537E+01,-1.86910E+01
2.96357E+00,2.89403E+00,-4.07385E+01,2.92128E+00,2.74335E+02,2.08986E+02
-2.20324E+00,-3.32515E+00,5.05197E+01,1.78893E+01,-3.92323E+02,-3.70113E+02
2.31474E-01,-4.16536E-01,-1.92429E+01,-1.17810E+01,1.58197E+02,1.67545E+02
1.00036E+00,-1.11107E-03,-2.58076E-03,1.36692E-02,7.08911E-03,-3.67730E-02
-3.24316E-01,1.15695E+00,9.45947E+00,-5.79675E+00,-5.33328E+01,-2.87528E+01
3.01159E+00,3.33981E+00,-4.40850E+01,-2.59468E+00,3.05217E+02,2.53266E+02
-2.16093E+00,-3.60244E+00,5.41141E+01,2.49638E+01,-4.28002E+02,-4.24126E+02
1.95645E-01,-4.28567E-01,-2.04742E+01,-1.44066E+01,1.70702E+02,1.87249E+02
1.00039E+00,-1.09011E-03,-3.65347E-03,1.31413E-02,1.53649E-02,-2.70335E-02
-3.53028E-01,1.05257E+00,1.04028E+01,-4.62034E+00,-6.10681E+01,-3.92186E+01
3.04907E+00,3.79831E+00,-4.73979E+01,-8.47391E+00,3.36394E+02,2.98829E+02
-2.10347E+00,-3.87306E+00,5.76385E+01,3.24147E+01,-4.63751E+02,-4.79318E+02
1.54174E-01,-4.49506E-01,-2.16756E+01,-1.71516E+01,1.83161E+02,2.07283E+02
1.00043E+00,-1.06114E-03,-4.77340E-03,1.25001E-02,2.40713E-02,-1.65728E-02
-3.80207E-01,9.40778E-01,1.13464E+01,-3.34423E+00,-6.89259E+01,-5.00549E+01
3.07493E+00,4.26732E+00,-5.06567E+01,-1.47102E+01,3.67721E+02,3.45521E+02
-2.02970E+00,-4.13334E+00,6.10696E+01,4.02298E+01,-4.99399E+02,-5.35495E+02
1.06694E-01,-4.80932E-01,-2.28393E+01,-2.00106E+01,1.95512E+02,2.27574E+02
1.00047E+00,-1.02377E-03,-5.93877E-03,1.17393E-02,3.32021E-02,-5.38025E-03
-4.05577E-01,8.21834E-01,1.22846E+01,-1.96909E+00,-7.68675E+01,-6.12207E+01
3.08808E+00,4.74434E+00,-5.38383E+01,-2.12945E+01,3.99035E+02,3.93158E+02
-1.93842E+00,-4.37923E+00,6.43816E+01,4.83938E+01,-5.34750E+02,-5.92428E+02
5.28167E-02,-5.24544E-01,-2.39565E+01,-2.29768E+01,2.07685E+02,2.48037E+02
1.00051E+00,-9.77691E-04,-7.14416E-03,1.08563E-02,4.27196E-02,6.50699E-03
-4.28882E-01,6.96049E-01,1.32126E+01,-4.94649E-01,-8.48638E+01,-7.26907E+01
3.08745E+00,5.22659E+00,-5.69235E+01,-2.82201E+01,4.30208E+02,4.41611E+02
-1.82849E+00,-4.60634E+00,6.75530E+01,5.68935E+01,-5.69648E+02,-6.49946E+02
-7.81506E-03,-5.82177E-01,-2.50201E+01,-2.60444E+01,2.19621E+02,2.68607E+02
1.00055E+00,-9.22411E-04,-8.38967E-03,9.83942E-03,5.26394E-02,1.91419E-02
-4.49833E-01,5.63865E-01,1.41246E+01,1.07680E+00,-9.28737E+01,-8.44177E+01
3.07189E+00,5.71077E+00,-5.98888E+01,-3.54715E+01,4.61066E+02,4.90670E+02
-1.69869E+00,-4.80968E+00,7.05579E+01,6.57062E+01,-6.03887E+02,-7.07793E+02
-7.55851E-02,-6.55871E-01,-2.60213E+01,-2.92042E+01,2.31249E+02,2.89188E+02
1.00059E+00,-8.58152E-04,-9.66020E-03,8.70229E-03,6.28328E-02,3.23376E-02
-4.68192E-01,4.25668E-01,1.50174E+01,2.74763E+00,-1.00883E+02,-9.64012E+01
3.04043E+00,6.19366E+00,-6.27207E+01,-4.30487E+01,4.91532E+02,5.40291E+02
-1.54799E+00,-4.98431E+00,7.33810E+01,7.48252E+01,-6.37365E+02,-7.65885E+02
-1.50810E-01,-7.47659E-01,-2.69551E+01,-3.24524E+01,2.42528E+02,3.09744E+02
1.00062E+00,-7.84582E-04,-1.09553E-02,7.43533E-03,7.33106E-02,4.61358E-02
-4.83646E-01,2.82050E-01,1.58842E+01,4.51493E+00,-1.08846E+02,-1.08589E+02
2.99185E+00,6.67128E+00,-6.53938E+01,-5.09329E+01,5.21419E+02,5.90246E+02
-1.37512E+00,-5.12439E+00,7.59949E+01,8.42238E+01,-6.69862E+02,-8.23945E+02
-2.33878E-01,-8.59879E-01,-2.78125E+01,-3.57783E+01,2.53381E+02,3.30173E+02
1.00066E+00,-7.01484E-04,-1.22684E-02,6.03454E-03,8.40298E-02,6.05018E-02
-4.95899E-01,1.33690E-01,1.67195E+01,6.37574E+00,-1.16723E+02,-1.20936E+02
2.92499E+00,7.13921E+00,-6.78861E+01,-5.91051E+01,5.50562E+02,6.40337E+02
-1.17889E+00,-5.22361E+00,7.83760E+01,9.38748E+01,-7.01185E+02,-8.81728E+02
-3.25146E-01,-9.95032E-01,-2.85858E+01,-3.91713E+01,2.63738E+02,3.50385E+02
1.00070E+00,-6.08729E-04,-1.35873E-02,4.50061E-03,9.48991E-02,7.53277E-02
-5.04647E-01,-1.86445E-02,1.75175E+01,8.32729E+00,-1.24477E+02,-1.33405E+02
2.83871E+00,7.59269E+00,-7.01766E+01,-6.75463E+01,5.78809E+02,6.90388E+02
-9.58143E-01,-5.27523E+00,8.05020E+01,1.03751E+02,-7.31155E+02,-9.39013E+02
-4.24964E-01,-1.15578E+00,-2.92679E+01,-4.26205E+01,2.73537E+02,3.70298E+02
1.00073E+00,-5.06594E-04,-1.49195E-02,2.82752E-03,1.05984E-01,9.07179E-02
-5.09498E-01,-1.74114E-01,1.82691E+01,1.03640E+01,-1.32042E+02,-1.45912E+02
2.73154E+00,8.02658E+00,-7.22315E+01,-7.62278E+01,6.05902E+02,7.40063E+02
-7.11399E-01,-5.27209E+00,8.23375E+01,1.13815E+02,-7.59479E+02,-9.95409E+02
-5.33773E-01,-1.34491E+00,-2.98474E+01,-4.61115E+01,2.82675E+02,3.89772E+02
1.00076E+00,-3.94938E-04,-1.62461E-02,1.01616E-03,1.17156E-01,1.06529E-01
-5.10179E-01,-3.31690E-01,1.89693E+01,1.24811E+01,-1.39382E+02,-1.58412E+02
2.60247E+00,8.43503E+00,-7.40330E+01,-8.51228E+01,6.31692E+02,7.89166E+02
-4.37673E-01,-5.20620E+00,8.38642E+01,1.24030E+02,-7.85984E+02,-1.05067E+03
-6.51850E-01,-1.56551E+00,-3.03186E+01,-4.96307E+01,2.91091E+02,4.08719E+02
1.00079E+00,-2.74014E-04,-1.75614E-02,-9.33046E-04,1.28369E-01,1.22714E-01
-5.06317E-01,-4.90282E-01,1.96102E+01,1.46725E+01,-1.46440E+02,-1.70838E+02
2.45014E+00,8.81187E+00,-7.55530E+01,-9.42005E+01,6.55963E+02,8.37419E+02
-1.35635E-01,-5.06921E+00,8.50532E+01,1.34356E+02,-8.10425E+02,-1.10448E+03
-7.79579E-01,-1.82077E+00,-3.06725E+01,-5.31628E+01,2.98700E+02,4.27020E+02
1.00082E+00,-1.43895E-04,-1.88590E-02,-3.02404E-03,1.39581E-01,1.39236E-01
-4.97537E-01,-6.48585E-01,2.01838E+01,1.69296E+01,-1.53157E+02,-1.83108E+02
2.27321E+00,9.15011E+00,-7.67630E+01,-1.03420E+02,6.78484E+02,8.84496E+02
1.95990E-01,-4.85189E+00,8.58760E+01,1.44742E+02,-8.32543E+02,-1.15645E+03
-9.17319E-01,-2.11420E+00,-3.09003E+01,-5.66889E+01,3.05413E+02,4.44542E+02
1.00085E+00,-5.30411E-06,-2.01239E-02,-5.24780E-03,1.50683E-01,1.55963E-01
-4.83526E-01,-8.05262E-01,2.06846E+01,1.92461E+01,-1.59496E+02,-1.95177E+02
2.07059E+00,9.44256E+00,-7.76442E+01,-1.12749E+02,6.99103E+02,9.30193E+02
5.58224E-01,-4.54471E+00,8.63141E+01,1.55146E+02,-8.52163E+02,-1.20633E+03
-1.06534E+00,-2.44938E+00,-3.09967E+01,-6.01932E+01,3.11168E+02,4.61192E+02
1.00087E+00,1.42049E-04,-2.13455E-02,-7.61031E-03,1.61591E-01,1.72814E-01
-4.63920E-01,-9.58725E-01,2.11056E+01,2.16134E+01,-1.65402E+02,-2.06974E+02
1.84100E+00,9.68106E+00,-7.81727E+01,-1.22145E+02,7.17614E+02,9.74218E+02
9.52222E-01,-4.13715E+00,8.63446E+01,1.65515E+02,-8.69054E+02,-1.25378E+03
-1.22395E+00,-2.83024E+00,-3.09548E+01,-6.36564E+01,3.15887E+02,4.76849E+02
1.00089E+00,2.96049E-04,-2.25199E-02,-1.00947E-02,1.72304E-01,1.89770E-01
-4.38356E-01,-1.10758E+00,2.14388E+01,2.40250E+01,-1.70813E+02,-2.18421E+02
1.58322E+00,9.85807E+00,-7.83221E+01,-1.31575E+02,7.33783E+02,1.01626E+03
1.37908E+00,-3.61927E+00,8.59420E+01,1.75804E+02,-8.82957E+02,-1.29844E+03
-1.39342E+00,-3.26049E+00,-3.07672E+01,-6.70618E+01,3.19480E+02,4.91383E+02
1.00090E+00,4.56346E-04,-2.36386E-02,-1.27003E-02,1.82771E-01,2.06789E-01
-4.06465E-01,-1.24997E+00,2.16774E+01,2.64711E+01,-1.75677E+02,-2.29447E+02
1.29603E+00,9.96437E+00,-7.80698E+01,-1.40992E+02,7.47404E+02,1.05603E+03
1.83985E+00,-2.97940E+00,8.50852E+01,1.85958E+02,-8.93641E+02,-1.33995E+03
-1.57401E+00,-3.74444E+00,-3.04281E+01,-7.03888E+01,3.21870E+02,5.04671E+02
1.00091E+00,6.22359E-04,-2.46850E-02,-1.54269E-02,1.92829E-01,2.23669E-01
-3.67940E-01,-1.38391E+00,2.18175E+01,2.89415E+01,-1.79969E+02,-2.40018E+02
9.78442E-01,9.99015E+00,-7.74041E+01,-1.50349E+02,7.58380E+02,1.09335E+03
2.33533E+00,-2.20518E+00,8.37655E+01,1.95915E+02,-9.00998E+02,-1.37812E+03
-1.76589E+00,-4.28658E+00,-2.99358E+01,-7.36156E+01,3.23018E+02,5.16638E+02
1.00091E+00,7.92875E-04,-2.56477E-02,-1.82659E-02,2.02416E-01,2.40347E-01
-3.22419E-01,-1.50732E+00,2.18520E+01,3.14257E+01,-1.83632E+02,-2.50054E+02
6.29314E-01,9.92529E+00,-7.63035E+01,-1.59597E+02,7.66494E+02,1.12791E+03
2.86640E+00,-1.28398E+00,8.19642E+01,2.05617E+02,-9.04788E+02,-1.41257E+03
-1.96925E+00,-4.89150E+00,-2.92855E+01,-7.67206E+01,3.22845E+02,5.27151E+02
1.00091E+00,9.66485E-04,-2.65182E-02,-2.12089E-02,2.11471E-01,2.56746E-01
-2.69552E-01,-1.61793E+00,2.17751E+01,3.39120E+01,-1.86616E+02,-2.59484E+02
2.47591E-01,9.75896E+00,-7.47505E+01,-1.68681E+02,7.71562E+02,1.15941E+03
3.43387E+00,-2.02402E-01,7.96675E+01,2.14997E+02,-9.04813E+02,-1.44297E+03
-2.18422E+00,-5.56401E+00,-2.84744E+01,-7.96806E+01,3.21284E+02,5.36093E+02
1.00090E+00,1.14163E-03,-2.72849E-02,-2.42457E-02,2.19912E-01,2.72759E-01
-2.08991E-01,-1.71333E+00,2.15806E+01,3.63880E+01,-1.88871E+02,-2.68235E+02
-1.67743E-01,9.47970E+00,-7.27286E+01,-1.77547E+02,7.73399E+02,1.18758E+03
4.03843E+00,1.05357E+00,7.68633E+01,2.23990E+02,-9.00876E+02,-1.46897E+03
-2.41092E+00,-6.30912E+00,-2.75000E+01,-8.24716E+01,3.18270E+02,5.43346E+02
1.00088E+00,1.31654E-03,-2.79362E-02,-2.73643E-02,2.27648E-01,2.88269E-01
-1.40385E-01,-1.79089E+00,2.12627E+01,3.88400E+01,-1.90344E+02,-2.76230E+02
-6.17665E-01,9.07535E+00,-7.02219E+01,-1.86133E+02,7.71818E+02,1.21209E+03
4.68070E+00,2.49864E+00,7.35409E+01,2.32521E+02,-8.92777E+02,-1.49023E+03
-2.64940E+00,-7.13202E+00,-2.63610E+01,-8.50681E+01,3.13739E+02,5.48788E+02
1.00086E+00,1.48934E-03,-2.84596E-02,-3.05545E-02,2.34591E-01,3.03168E-01
-6.33861E-02,-1.84775E+00,2.08152E+01,4.12512E+01,-1.90978E+02,-2.83376E+02
-1.10309E+00,8.53286E+00,-6.72151E+01,-1.94369E+02,7.66617E+02,1.23259E+03
5.36119E+00,4.14839E+00,6.96913E+01,2.40511E+02,-8.80309E+02,-1.50636E+03
-2.89970E+00,-8.03816E+00,-2.50563E+01,-8.74417E+01,3.07622E+02,5.52280E+02
1.00083E+00,1.65771E-03,-2.88470E-02,-3.38011E-02,2.40678E-01,3.17362E-01
2.23629E-02,-1.88091E+00,2.02316E+01,4.36050E+01,-1.90711E+02,-2.89575E+02
-1.62493E+00,7.83872E+00,-6.36928E+01,-2.02185E+02,7.57583E+02,1.24872E+03
6.08035E+00,6.01881E+00,6.53055E+01,2.47876E+02,-8.63252E+02,-1.51693E+03
-3.16179E+00,-9.03311E+00,-2.35856E+01,-8.95638E+01,2.99851E+02,5.53685E+02
1.00079E+00,1.81924E-03,-2.90832E-02,-3.70879E-02,2.45791E-01,3.30698E-01
1.17187E-01,-1.88709E+00,1.95065E+01,4.58826E+01,-1.89486E+02,-2.94734E+02
-2.18391E+00,6.97844E+00,-5.96443E+01,-2.09503E+02,7.44530E+02,1.26012E+03
6.83839E+00,8.12682E+00,6.03812E+01,2.54528E+02,-8.41421E+02,-1.52156E+03
-3.43558E+00,-1.01226E+01,-2.19506E+01,-9.14044E+01,2.90368E+02,5.52871E+02
1.00075E+00,1.97121E-03,-2.91599E-02,-4.03994E-02,2.49868E-01,3.43083E-01
2.21414E-01,-1.86285E+00,1.86345E+01,4.80646E+01,-1.87245E+02,-2.98750E+02
-2.78070E+00,5.93700E+00,-5.50600E+01,-2.16243E+02,7.27260E+02,1.26641E+03
7.63543E+00,1.04897E+01,5.49176E+01,2.60377E+02,-8.14622E+02,-1.51985E+03
-3.72092E+00,-1.13126E+01,-2.01538E+01,-9.29324E+01,2.79114E+02,5.49699E+02
1.00069E+00,2.11063E-03,-2.90649E-02,-4.37166E-02,2.52811E-01,3.54381E-01
3.35343E-01,-1.80450E+00,1.76110E+01,5.01296E+01,-1.83930E+02,-3.01520E+02
-3.41584E+00,4.69858E+00,-4.99338E+01,-2.22321E+02,7.05590E+02,1.26721E+03
8.47134E+00,1.31256E+01,4.89188E+01,2.65326E+02,-7.82678E+02,-1.51137E+03
-4.01757E+00,-1.26092E+01,-1.81997E+01,-9.41154E+01,2.66036E+02,5.44032E+02
1.00062E+00,2.23438E-03,-2.87912E-02,-4.70251E-02,2.54570E-01,3.64523E-01
4.59241E-01,-1.70813E+00,1.64324E+01,5.20553E+01,-1.79492E+02,-3.02947E+02
-4.08966E+00,3.24663E+00,-4.42655E+01,-2.27646E+02,6.79364E+02,1.26217E+03
9.34577E+00,1.60530E+01,4.23955E+01,2.69276E+02,-7.45446E+02,-1.49575E+03
-4.32520E+00,-1.40184E+01,-1.60952E+01,-9.49201E+01,2.51093E+02,5.35741E+02
1.00054E+00,2.33856E-03,-2.83172E-02,-5.02884E-02,2.54953E-01,3.73216E-01
5.93303E-01,-1.56969E+00,1.50988E+01,5.38223E+01,-1.73905E+02,-3.02969E+02
-4.80218E+00,1.56424E+00,-3.80667E+01,-2.32140E+02,6.48516E+02,1.25106E+03
1.02579E+01,1.92906E+01,3.53725E+01,2.72140E+02,-7.02882E+02,-1.47274E+03
-4.64333E+00,-1.55463E+01,-1.38520E+01,-9.53176E+01,2.34279E+02,5.24749E+02
1.00045E+00,2.41947E-03,-2.76379E-02,-5.34904E-02,2.53921E-01,3.80403E-01
7.37729E-01,-1.38486E+00,1.36088E+01,5.54064E+01,-1.67129E+02,-3.01495E+02
-5.55335E+00,-3.66307E-01,-3.13473E+01,-2.35709E+02,6.12931E+02,1.23355E+03
1.12069E+01,2.28579E+01,2.78734E+01,2.73815E+02,-6.54893E+02,-1.44199E+03
-4.97143E+00,-1.71992E+01,-1.14816E+01,-9.52743E+01,2.15573E+02,5.10941E+02
1.00035E+00,2.47268E-03,-2.67441E-02,-5.66022E-02,2.51390E-01,3.85940E-01
8.92647E-01,-1.14912E+00,1.19639E+01,5.67848E+01,-1.59141E+02,-2.98457E+02
-6.34279E+00,-2.56322E+00,-2.41272E+01,-2.38265E+02,5.72562E+02,1.20938E+03
1.21914E+01,2.67741E+01,1.99334E+01,2.74204E+02,-6.01463E+02,-1.40324E+03
-5.30880E+00,-1.89832E+01,-8.99983E+00,-9.47588E+01,1.94979E+02,4.94237E+02
1.00024E+00,2.49375E-03,-2.56255E-02,-5.95976E-02,2.47266E-01,3.89678E-01
1.05812E+00,-8.57709E-01,1.01681E+01,5.79340E+01,-1.49935E+02,-2.93803E+02
-7.16980E+00,-5.04529E+00,-1.64365E+01,-2.39716E+02,5.27420E+02,1.17838E+03
1.32096E+01,3.10591E+01,1.15998E+01,2.73212E+02,-5.42641E+02,-1.35631E+03
-5.65461E+00,-2.09042E+01,-6.42647E+00,-9.37406E+01,1.72527E+02,4.74581E+02
1.00011E+00,2.47775E-03,-2.42766E-02,-6.24493E-02,2.41487E-01,3.91511E-01
1.23421E+00,-5.05710E-01,8.22428E+00,5.88279E+01,-1.39487E+02,-2.87452E+02
-8.03353E+00,-7.83160E+00,-8.30429E+00,-2.39965E+02,4.77471E+02,1.14028E+03
1.42597E+01,3.57324E+01,2.91990E+00,2.70736E+02,-4.78440E+02,-1.30094E+03
-6.00791E+00,-2.29682E+01,-3.78184E+00,-9.21874E+01,1.48234E+02,4.51893E+02
9.99975E-01,2.41925E-03,-2.26882E-02,-6.51219E-02,2.33958E-01,3.91263E-01
1.42085E+00,-8.80256E-02,6.13985E+00,5.94426E+01,-1.27806E+02,-2.79367E+02
-8.93265E+00,-1.09415E+01,2.24068E-01,-2.38922E+02,4.22794E+02,1.09495E+03
1.53388E+01,4.08131E+01,-6.04035E+00,2.66685E+02,-4.08996E+02,-1.23702E+03
-6.36755E+00,-2.51806E+01,-1.09260E+00,-9.00713E+01,1.22161E+02,4.26151E+02
9.99822E-01,2.31254E-03,-2.08550E-02,-6.75815E-02,2.24607E-01,3.88796E-01
1.61799E+00,4.00698E-01,3.92095E+00,5.97493E+01,-1.14882E+02,-2.69475E+02
-9.86566E+00,-1.43946E+01,9.10466E+00,-2.36487E+02,3.63422E+02,1.04216E+03
1.64443E+01,4.63204E+01,-1.52151E+01,2.60954E+02,-3.34406E+02,-1.16435E+03
-6.73227E+00,-2.75466E+01,1.61430E+00,-8.73612E+01,9.43594E+01,3.97300E+02
9.99655E-01,2.15142E-03,-1.87711E-02,-6.97883E-02,2.13350E-01,3.83936E-01
1.82545E+00,9.65942E-01,1.57686E+00,5.97207E+01,-1.00725E+02,-2.57723E+02
-1.08305E+01,-1.82107E+01,1.82808E+01,-2.32560E+02,2.99460E+02,9.81786E+02
1.75726E+01,5.22723E+01,-2.45237E+01,2.53451E+02,-2.54850E+02,-1.08283E+03
-7.10059E+00,-3.00710E+01,4.30673E+00,-8.40296E+01,6.49087E+01,3.65326E+02
9.99473E-01,1.92924E-03,-1.64318E-02,-7.17008E-02,2.00109E-01,3.76517E-01
2.04302E+00,1.61334E+00,-8.82149E-01,5.93264E+01,-8.53384E+01,-2.44043E+02
-1.18250E+01,-2.24092E+01,2.76900E+01,-2.27039E+02,2.31010E+02,9.13629E+02
1.87198E+01,5.86859E+01,-3.38779E+01,2.44074E+02,-1.70518E+02,-9.92335E+02
-7.47086E+00,-3.27576E+01,6.94978E+00,-8.00487E+01,3.38953E+01,3.30205E+02
9.99276E-01,1.63895E-03,-1.38361E-02,-7.32738E-02,1.84822E-01,3.66376E-01
2.27037E+00,2.34860E+00,-3.44417E+00,5.85364E+01,-6.87344E+01,-2.28373E+02
-1.28460E+01,-2.70092E+01,3.72616E+01,-2.19825E+02,1.58211E+02,8.37547E+02
1.98811E+01,6.55768E+01,-4.31793E+01,2.32735E+02,-8.16486E+01,-8.92777E+02
-7.84119E+00,-3.56097E+01,9.50496E+00,-7.53943E+01,1.42399E+00,2.91936E+02
9.99064E-01,1.27324E-03,-1.09916E-02,-7.44696E-02,1.67497E-01,3.53455E-01
2.50713E+00,3.17754E+00,-6.09661E+00,5.73172E+01,-5.09157E+01,-2.10626E+02
-1.38905E+01,-3.20289E+01,4.69203E+01,-2.10810E+02,8.11857E+01,7.53321E+02
2.10515E+01,7.29586E+01,-5.23231E+01,2.19337E+02,1.15276E+01,-7.84022E+02
-8.20951E+00,-3.86295E+01,1.19315E+01,-7.00414E+01,-3.23996E+01,2.50503E+02
9.98836E-01,8.23950E-04,-7.89491E-03,-7.52274E-02,1.48028E-01,3.37491E-01
2.75270E+00,4.10593E+00,-8.82051E+00,5.56392E+01,-3.19293E+01,-1.90772E+02
-1.49543E+01,-3.74855E+01,5.65661E+01,-1.99903E+02,2.22903E-01,6.60953E+02
2.22247E+01,8.08420E+01,-6.11776E+01,2.03805E+02,1.08591E+02,-6.66179E+02
-8.57337E+00,-4.18175E+01,1.41794E+01,-6.39741E+01,-6.74045E+01,2.05974E+02
9.98592E-01,2.82966E-04,-4.55740E-03,-7.55030E-02,1.26429E-01,3.18414E-01
3.00639E+00,5.13948E+00,-1.15960E+01,5.34702E+01,-1.18151E+01,-1.68761E+02
-1.60324E+01,-4.33944E+01,6.60933E+01,-1.87010E+02,-8.43973E+01,5.60386E+02
2.33936E+01,8.92345E+01,-6.96034E+01,1.86068E+02,2.09119E+02,-5.39313E+02
-8.93007E+00,-4.51728E+01,1.61962E+01,-5.71780E+01,-1.03415E+02,1.58410E+02