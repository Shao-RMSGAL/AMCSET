# eSRIM (electron Stopping and Range of Ions in Matter)


<div style="text-align:center;">
  <img src="assets/images/1000electrons.svg" alt="1000 Electron Simulation", style="trim-path">
</div>

## Table of Contents  
- [Introduction](#introduction)
- [Dependencies](#dependencies)
- [Features](#features)
- [Building](#building)
- [Running](#running)
- [Planned Features](#planned-features)
- [Author](#author)
- [Acknowledgements](#Acknowledgements)
- [License](#license)

## Introduction

This is eSRIM, a version of Stopping and Range of Ions in Matter (SRIM) with added capabilities for electron bombardment simulation. It is written in C++ using the GNU standard libraries. Originally conceived of and written by Dr. Lin Shao, this is a program written to implement the methods used by Dr. Shao to efficiently perform electron bombardment calculations efficiently, a traditionally computationally expensive and slow task. 

Having the ability to simulate electron and ion bombardment is an invaluable tool in the field of materials science, as it enables the prediction of irradiation penetration depth, allowing researchers to correlate experimental data with simulated data.

## Features

This is a program for simulating ion bombardment in materials. It can be used for simulating electron bombardment.

Written in C++, it is a highly optimized program, which takes full advantage of multi-threading, combined with performance-oriented programming to achieve extremely high performance calculations of both ion bombardment and electron bombardment. 

Thanks to multi-threading support, thousands of simulations can be run in parallel and the performance can be tuned by setting the number of allowed threads in the settings.

### Ion Bombardment Simulation

Just like the classic SRIM program, this simulation can simulate the Stopping and Range of Ions in Matter. Using the [settings.txt](build/settings.txt) file, set any bombardment ion mass and Z-number, as well as any substrate mass and Z-number, and the program will simulate ion bombardment, outputting to a file with a name specified in the settings.

### Electron Bombardment Simulation

In addition to ion bombardment, this program can also simulate electron bombardment at high throughput. Through highly optimized Mott scattering calculations, a single electron bombardment can be completed in a fraction of a second, while this calculation would take several seconds or even minutes using traditional techniques.

This works by calculating the total Mott scattering cross-section for a group of many local substrate atom interactions, with the assumption that elastic energy loss is negligible over small distances. This assumption is demonstrated to be valid for groupings of up to local 1000 interactions. 

### Multithreading

This program employs multi-threading, which allows for multiple CPU cores to be employed generating bombardment data. Since each bombardment is independent, a separate thread can simulate an entire bombardment and write the simulation data to a file. This can be done with many threads simultaneously.

### Protected Output Files

The program has a built-in protection which ensures that data from previous runs is not accidentally overridden. When the program encounters a file that has the same name as the name specified in the settings file, it checks to see if it has a file marker at the end (this marker can also be customized). If a marker is present, the old file is renamed \<filename>\_\<Day>\_\<Month>\_\<Date>\_\<Hour>\_\<Minute>\_\<Second>\_.\<file extension>. For example, ```coordinateOutput.csv``` will become ```coordinateOutput_Thu_Mar_14_19-40-03_2024.csv``` if it is already present in the output directory when a new run is initiated.

### Data Toggles

There are several toggles available in the settings file that enables a user to select what data is reported by eSRIM. This is a valuable feature, as eSRIM can generate gigabytes of data if all the generated data is reported. Three toggles allow for reduction of output file size:
- _logStoppingPointOnly_ - This option only logs the final resting point of each particle, and does not record the path that the particle took to get to the stopping point. 
- _logEndOfFlyingDistance_ (__Electron only__) This option reduces the number of coordinates logged by not logging every coordinate within a local substrate interaction group, but only logging the last interaction within a group. This can reduce the data size by a factor of 

## Dependencies

### C++ Compiler

Building the program requires the C++ standard libraries, as well as a C++ compiler. This program was developed using the g++ and clang compiler. The compiler can be chosen by setting the ```-DCOMPILER``` flag when using cmake. On Linux, a C++ compiler and the associated standard libraries generally included with most distributions. This compiler is called g++, and it is all you will need to build and run this project.

### Google gtest (Developers only)

This project uses Google [gtest](https://github.com/google/googletest/) to handle test cases. The [cmake configuration file](CMakeFile.txt) automatically handles the installation of this dependency. If you encounter any gtest related issues on your system, ensure that it is installed using

```bash
sudo apt-get install libgtest-dev
```

The testing environment has been configured, but no tests have been run so far. To run the ["Trivial test"](tests/mainTest.cpp#TrivialTest)

### Supported Compilers

Three compilers are supported, g++, clang++, and icpx (Intel compiler).

- There is currently a linker error generated when trying to compile with clang. Please build with ipcx or g++.

## Building


**IMPORTANT: Please follow the build directions exactly to ensure successful compilation of the project. It is recoomended to copy and paste the instructions exactly**

Please build this project using a Linux. If you are on Windows, install [WSL2](https://learn.microsoft.com/en-us/windows/wsl/install) before continuing. Support for building on Windows is planned for future releases.

This project can be built using cmake. Once you clone the repository to your local machine, ensure you have cmake installed. If you are on a Debian Linux distribution, you can do this using
```bash
sudo apt install cmake
```

If you are using another Linux distribution or Windows, ensure you have cmake, a C++ compiler, and the C++ standard libraries installed on your system before continuing. 

Note that Windows is not officially supported, but it is not expected that this code would not compile and run on Windows.

To build the project, first nagivate to the build directory:

```bash
cd build
```

Then configure cmake

```bash    
cmake ..
```

You should see an output like this

```bash
-- The CXX compiler identification is GNU XX.X.X
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- Check for working CXX compiler: /usr/bin/c++ - skipped
-- Detecting CXX compile features
-- Detecting CXX compile features - done
-- The C compiler identification is GNU XX.X.X
-- Detecting C compiler ABI info
-- Detecting C compiler ABI info - done
-- Check for working C compiler: /usr/bin/gcc - skipped
-- Detecting C compile features
-- Detecting C compile features - done
-- Found Python3: /usr/bin/python3.XX (found version "3.XX.XX") found components: Interpreter
-- Performing Test CMAKE_HAVE_LIBC_PTHREAD
-- Performing Test CMAKE_HAVE_LIBC_PTHREAD - Success
-- Found Threads: TRUE
-- Using /usr/bin/gcc C compiler
-- Using /usr/bin/c++ C++ compiler
-- Production mode is enabled.
-- Configuring done (2.5s)
-- Generating done (0.0s)
-- Build files have been written to: /path/to/project/eSRIM/build
```

You can also set flags to configure cmake to use a different compiler. Read the [CMakeLists.txt](CMakeLists.txt) file to learn more.

Once cmake is configured, build the project using

```bash
cmake --build .
```

This will build the project and place the binary in the [/bin] folder.

### Testing (For Developers only)

Test cases can be found in the [tests](tests) directory. To run the all the test suites, run 

```bash
ctest
```

once you've [built](#building) the project.


## Running

Run the program using the following command

```bash
./bin/eSRIM
```

By default, the program does not produce any output, it simply reads the [settings.txt](build/settings.txt) file and completes a simulation based on those settings before exiting. 

eSRIM supports a number of command-line arguments. You can see all the available arguments using the _-h_ flag

```bash
./bin/eSRIM -h
```

## Cleaning the Build Directory

If you wish to clear the build directory and associated files, run 
```bash
make clean
```

This will automatically delete any executables and intermediate files generated by the build process.

### Output

By default, the program outputs a .csv file. The name and file extension can be changed in the [settings.txt](build/settings.txt) file. 

The output has 6 columns. Here is an example

#### coordinateOutput.csv

```
bombardmentID,  particleID, depth,  x,            y,            z
3,              0,          0,      1.70705e+08,  -6.06866e+07, 4.03509e+08
4,              0,          0,      -4.90489e+08, -2.02786e+08, 2.654e+08
1,              0,          0,      -1.17807e+08,  8.74504e+06, 2.74198e+08
```

- The _bombardmentID_ is the simulation number. Each bombardment is independent from the other.
- The _particleID_ is an identifier for each particle generated in a simulation. _particleID_ 0 is the original particle, while all other particles have an ID of 1 or more. 
- The depth is the cascade generation number of the particle. When a particle displaces a substrate atom, that substrate particle is given a _particleID_ and a depth of 1 more than the parent particle. 
- _x_,_y_,and _z_ are the respective coordinates.

#### Graphing the Output

A python script is located in the /python directory called [graphing.py](python/graphing.py). To use it, ensure you have [python](https://www.python.org/) installed. In addition, you will need to install the [pandas](https://pandas.pydata.org/) and [matplotlib](https://matplotlib.org/) libraries.

Install these using

```bash
pip install pandas
pip install matplotlib
```

After you've run eSRIM, you can generate a plot using this script like so (run this from the /build directory).

```bash
../python/graphing.py
```


## Planned Features

Below is a list of planned features that are not yet implemented in the program. There are three different types of planned features
- (Required) Sections marked as Required are features that should be an inherent part of the program -- they cannot be controlled by the user. They are inherent program behavior that should be the only behavior of the program, with no option to disable it.
- (Option) Sections marked as Option are features that can be toggled by the user. These are optional features that are not critical to the program's function.
- (Development) These are  codebase features that should be implemented to make development go by more smoothly. They do not affect the final program's behavior, only the development process.

### Building on Windows (Development)

Currently this project is only tested to build and run on Linux. Windows support should be added.

### Add Planned Features and Issues to GitHub (Development)

While this file tracks planned features, they should also be added to the [GitHub](https://github.com/Shao-RMSGAL/eSRIM/issues) page as issues. Similarly, issues should also be added to GitHub.

### Signal Handling (Required)

The program should be able to interpret and respond to signals such as SIGINT (Usually sent my Ctrl+C) by stopping any file I/O and writing a file end marker to the output file.

### Auxillary Output Information (Option)

In addition to the primary output file, an auxillary information file should be created which contains information about the simulation

### Graphical User Interface (Option)

A graphical interface should be created upon running the program. This should be the default behavior, unless the program is specifically invoked through the command line with a flag disabling it. 

The user should be able to specify settings, such as the substrate isotope/element, as well as the ion isotope/element (if an ion). They should also be able to run the program and view the data, all inside the interface. Ideally this should be implemented using a modern UI interface. Some promising candidates are:

### High and Low Energy Electron Adjustment (Required)

Currently, the code does not change the number of grouped substrate interactions for electrons with very small or very large energy. This can prove to be a problem for very high energy electrons, where total flying distance become exceedingly large, and for low energies, where the lack of granular elastic energy calculations could result in a negative electron energy. The number of flying distances should be reduced in both cases.

### Output File Compression (Option)

The program should have an option to enable/disable output file compression. It should be able to compress output to .7z or .zip files.

### Alternative Output Filetypes (Option)

Other output filetypes should be supported. While ```.csv``` files are ubiquitous, they are not good for extremely large files (like those that can be generated by this program) as they are very simple files (text separated by commas, as indiacted by the file extension, _.CommaSeparatedList_, or ```.csv```).

### Energy Bounds Checking (Required)

Currently the program does not check the bounds on the energy of particles before performing calculations. This is especially important for electron bombardment calculations, where several of the calculations made are only valid for certain ranges of energy. The program should check that the user-provided energy is within the energy bounds of the calculation.

### Test Case Writing (Development)

The codebase should be tested using Google gtest to ensure that the program works between updates. This should be integrated with GitHub to ensure that code that does not function correctly is not accidentally merged into the main branch.

### GPU Optimization (Option)

The code should be able to leverage the GPU to perform parallelized calculations on a large scale. The most promising candidate for this option is the [function for generating Mott Differential Cross Sections](src/eSRIM_classes.cpp#getMottDifferentialCrossSection). This function is a major bottleneck for electron bombardment simulation (determined by performance analysis). Parallelizing this operation using a GPU may result in drastic speedups. 

Due to the non-standardized nature of GPU-Accelerated programming (NVidea uses [CUDA](https://developer.nvidia.com/cuda-toolkit), AMD uses [ROCm](https://www.amd.com/en/products/software/rocm.html), Intel uses [OpenCL](https://www.intel.com/content/www/us/en/developer/articles/tool/tools-for-opencl-applications.html)) multiple GPU-Acceleration algorithms may need to be implemented depending on the system configuration. 

[SYCL](https://www.khronos.org/sycl/) is a promising library which attempts to unify the many different GPU vendors under a single development library. Using this can avoid multiple implementations, but may be less performant.

### Mott Scattering Differential Cross Section Calculation Caching (Option)

The Mott scattering differential cross section calculation is a major bottleneck in the code. It involves exponentiation and many array accesses to calculate a single differential cross section for a particular energy.

However, this is a deterministic function, albeit a complicated one. The results for certain energies can be cached and used by multiple threads to reduce the computational workload for large simulations. Bins should be sufficiently small to ensure that accuracy is not lost. A promising data structure for this purpose is std::unordered_map.

### Unified Default Settings (Development)

Currently, the default settings (in [constants.h](include/constants.h)) and default settings file image (in [utilities.h](include/utilities.h)). These should be unified so that a change in the [default settings](include/constants.h#electronEnergy) in constants.h is reflected in the [settingsMessage](include/utilities.h#settingsMessage)<!--First variable in the Defaults namespace--> variable.

### Integrated Mott Scattering Parameters and Electron Screening Parameters (Option)

The Mott Scattering Parameters and Electron Screening Parameters are currently read in from a separate file. While this functionality is useful, an option should also be provided to regenerate the files from the executable in case the files are lost. This is similar to the existing settings preservation settings, where if the ["settings.txt"](/build/settings.txt) cannot be located, a new one is generated.

### Clean Up Main Function Wrapping and Input Stream Handling (Development)

In order to get tests to run properly, the main function was wrapped in a function called eSRIM(). and a [separate header file](include/eSRIM.h) was made for it. However, the input stream (which is passed as an argument to the wrapper) is not properly passed through the InputFields class, and is instead passed as a parameter to several functions which read the standard output, either directly or indirectly. This behavior is undesirable as it makes function signatures unnecessarily complicated and should be corrected.

Add the input stream reference (std::istream& inputStream) as a member in the InputFields class so that it can be better handled.

### Implement Object-Oriented Programming for I/O functions (Development)

Currently several user-interface functions (those in [utilities.cpp](src/utilities.cpp)) are not wrapped with objects, and are in the global scope. This is bad practice that can lead to disorganized code. The I/O functions should be wrapped into an IOHandler class which manages I/O operations through member variables and functions. 

### Proper Licensing (Development)

While the project has a license file, the project is not yet properly licensed. It needs a copyright disclaimer, as well as proper header files with the copyright notice. Notice should be aquired from Dr. Shao and any other invested parties before proper licensing is enforced. 

It should adhere to the [GNU Licensing Guidelines](https://www.gnu.org/licenses/gpl-howto.html).

## Known issues

### Command line progress bar
- The command line progress bar does not properly clear at the end of a run. 

![Progress Bar Bug](assets/images/progressBarBug.png)

### Markdown file image bug

- Image at the top of the markdown file as axis multipliers of "10e8", but are only partially visible. The graph was generated using matplotlib.pylot.

## Author

<div style="text-align:center;">
  <img src="assets/images/portrait.jpg" alt="Nathaniel Thomas" width="300">
</div>

I am a chemical engineer working on my PhD in nuclear engineering. I have a passion for molten salt corrosion research, with hopes of applying my computational knowledge and skills towards solving corrosion problems in molten salt reactors. This particular project is important to me, as it is relevant to predicting irradiation effects in nuclear materials, a particularly important problem in my field.

In my free time, I love to play piano, swim, play video games, and program. 
- Name:     Nathaniel Thomas
- Email:    nathaniel@swbell.net
- GitHub:   [github.com/At11011](https://github.com/At11011)

## License

This projoect will be licensed under the [GNU General Public License v2.0](LICENSE) - see the [LICENSE](LICENSE) file for details. 

Currently the licensing is not yet properly implemented. See proper [Proper Licensing (Development)](#proper-licensing-development) for details.


## Acknowledgements

All credit goes to Dr. Lin Shao at Texas A&M University for providing the original implementation of this software using QuickBASIC. His program included ion and electron bombardment support, as well as the novel flying-distance grouping method and the midpoint energy approximation method that drastically improved the performance of electron bombardment calculations while still maintaining high accuracy.

This program builds upon his original work by improving the ease of use (settings file instead of CLI input), as well as performance improvements.

Thanks to Kenneth Cooper and Benjamin Mejia Diaz for their help in reviewing this code.

---
### Last modified March 16, 2024