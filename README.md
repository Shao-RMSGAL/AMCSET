# eSRIM (electron Stopping and Range of Ions in Matter)

## Table of Contents
- [Introduction](#introduction)
- [Dependencies](#dependencies)
- [Features](#features)
- [Building](#building)
- [Running](#running)
- [Planned Features](#planned_features)
- [Author](#author)
- [Acknowledgements](#Acknowledgements)
- [License](#license)

## Introduction

This is eSRIM, a version of Stopping and Range of Ions in Matter (SRIM) with added capabilities for electron bombardment simulation. It is implemented in C++.

## Features

This is a program for simulating ion bombardment in materials. It can be used for simulating electron bombardment.

Written in C++, it is a highly optimized program, which takes full advantage of multi-threading, combined with performance-oriented programming to achieve extremely high performance calculations of both ion bombardment and electron bombardment. 

Thanks to multi-threading support, thousands of simulations can be run in parallel and the performance can be tuned by setting the number of allowed threads in the settings.

### Ion Bombardment Simulation

### Electron Bombardment Simulation

### Multithreading

### Protected Output Files



## Dependencies
List any dependencies required to build and run the program.

- Boost Libraries:

### Installation:

#### Debian (Ubuntu, Mint, etc.)

```bash
sudo apt-get install libboost-all-dev
```

## Building

The 
    
Configure cmake using the following:

```bash    
cmake -DDEBUG_MODE=OFF -DPERF=OFF -DPROD=ON .
```

Options:
    -DDEBUG_MODE=<option>   Enable or disable debug messages. Also adds -g flag for debugging using gdb or valgrind.
    -DPERF=<option>        Enable performance analysis mode. For use with valgrind callgrind. This is ifnored if -DDEBUG_MODE is enabled.
    -DPROD=<option>         For production-mode compilation. Enables optimizations. This is ignored if -DDEBUG_MODE or -DPERF is enabled.

Make the project (compiles program based on cmake configuration above)

```bash
make .
```

## Running

Run the project

```bash
./bin/eSRIM [options] [arguments]
```

Use -h or --help option to see usage information.  

### Output

By default, the program outputs a .csv file. The name and file extension can be changed in the settings.txt file. 

## Planned Features

### Signal Handling

The program should be able to interpret and respond to signals such as SIGINT (Usually sent my Ctrl+C) by stopping any file I/O and writing a file end marker to the output file.

### Auxillary Output Information

In addition to the primary output file, an auxillary information file should be created which contains information about the simulation

### Graphical User Interface

A graphical interface should be created upon running the program. This should be the default behavior, unless the program is specifically invoked through the command line with a flag disabling it. 

The user should be able to specify settings, such as the substrate isotope/element, as well as the ion isotope/element (if an ion).

### Reorganized codebaes

Currently the project is part of the "DR.SHAO-RMSLCF" project on GitHub, under the subdirectory /src. This project should be ported to a separate project called eSRIM and have it's own project as part of the Dr. Shao RMSLCF GitHub Repository ([github.com/Shao-RMSGAL](https://github.com/Shao-RMSGAL)).

### Output File Compression

The program should have an option to enable/disable output file compression. It should be able to compress output to .7z or .zip files.

### Alternative Output File Types



## Author

<div style="text-align:center;">
  <img src="photo.jpg" alt="Nathaniel THomas" width="300">
</div>

I am a chemical engineer working on my PhD in nuclear engineering. I have a passion for molten salt corrosion research, with hopes of applying my computational knowledge and skills towards solving corrosion problems in molten salt reactors. In my free time, I love to play piano, swim, play video games, and program. 
- Name:     Nathaniel Thomas
- Email:    nathaniel@swbell.net
- GitHub:   [github.com/At11011](https://github.com/At11011)

# License

This project is licensed under the [GNU General Public License v2.0](LICENSE) - see the [LICENSE](LICENSE) file for details.


# Acknowledgements

All credit goes to Dr. Lin Shao at Texas A&M University for providing the original implementation of this software using QuickBASIC. His program included ion and electron bombardment support, as well as the novel flying-distance grouping method and the midpoint energy approximation method that drastically improved the performance of electron bombardment calculations while still maintaining high accuracy.

Thanks to Kenneth Cooper and Benjamin Mejia Diaz for their help in reviewing this code.

---
Last modified: March 14, 2024