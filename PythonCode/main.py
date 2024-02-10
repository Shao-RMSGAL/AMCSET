"""
File: main.py
Author:         Nathaniel Thomas
Date created:   February 9th, 2024
Last modified:  February 9th, 2024
Group:          Shao Group, Texas A&M Nuclear Engineering

Description:
- This is the primary file for Dr. Shao's particle collision simulation software. TODO: Provide a more detailed description
"""

import numpy as np
import shao_utils as ut

# Depth of interest in Angstroms
# e_range = 80000000

# Interval number for the region of interest
e_interval = 30

# Range to averageing neighboring data
e_ave_range = 1

e_deltaLat = ut.SUBDENSITY_DEFAULT/e_interval
e_deltaDep = ut.SUBDENSITY_DEFAULT/e_interval

# Electrons are assumed to stop at an energy of 00 eV
e_stopping = 1

# Ions are assumed to stop at an energy of 20 eV
ion_stopping = 0.04

# Threshold displacement energy
ion_Ed = 0.04

# Number of angle intervals from 0 to 180 degrees for electron bombardment
Divisor = 1000

fly0 = 1000
fly1 = 1000
flyC = fly1
flyJudge = 1
sizeScan = 100
Resolution = 10
scandelta = 20
EMAX = 10000

# Switch value. '0' is for an explicit method, '1' is for the middle-point method, and '2' is for back (implicit) method
switch = 1

# The code has several type declaration statements at this point. Python does not require this.

Ltotal = 0.0

# Matrix for storing Mott potentials (TODO: Include all elements up to 118)
mott = np.zeros((96,5,6))

# TODO: Write code for reading mott potential values in from a text file.

# Matrix for storing screening potentials up to element 92. TODO: Change to read up to element 118
scr = np.zeros((92,6))

# TODO: Write code for reading screening potential values

if __name__ == "__main__":
    # Various output files
    out_1 = open("./output/PRANGEE40.dat", 'w')
    out_2 = open("./output/simulEd40.dat", 'w')
    out_3 = open("./output/simulvimageEd40.dat", 'w')
    out_4 = open("./output/#DINTLimageEd40.dat", 'w')
    out_5 = open("./output/eRangeimageEd40.dat", 'w')
    out_6 = open("./output/ProjectedeRangeimageEd40.dat", 'w')
    out_7 = open("./output/e-outimageEd40.dat", 'w')
    out_9 = open("./output/Projected_e_rangeImageEd40.txt", 'w')
    out_10 = open("./output/20imageEd40.txt", 'w')
    out_11 = open("./output/crossimageEd40.txt", 'w')
    out_12 = open("./output/imageEd40.txt", 'w')
    out_13 = open("./output/VacancyEd40.txt", 'w')
    out_14 = open("./output/dedx.txt", 'w')

    # iontype determines whether it is an electron bombardment or ion bombardment. The default is 0
    is_electron = ut.type_check_input("Is this an electron bombardment?", bool, ut.IS_ELECTRON_DEFAULT)

    if ut.DEBUG:
        print("Is electron:", is_electron)
    # In the case of electron bombardment, ask for energy:
    if is_electron:
        # Ask for the electron input if electron bombardment is being simulated
        ElecE0 = ut.type_check_input("Input electron energy in keV", float, ut.ELECE0_DEFAULT)
        if ut.DEBUG:
            print("Energy:", ElecE0)
    else:
        # Ask for ion information if not simulating electron bombardment
        massp = ut.type_check_input("Input incident atom mass in amu", float, ut.MASSP_DEFAULT)
        ZP = ut.type_check_input("Input incident charge Z in e", int, ut.ZP_DEFAULT)
        INELAB = ut.type_check_input("Input incident atom energy in keV", float, ut.INELAB_DEFAULT)
        if ut.DEBUG:
            print("Mass:", massp)
            print("Charge:", ZP)
            print("Energy:", INELAB)

    # Continue on to inputting substrate details
    MASSSUB = ut.type_check_input("Input substrate mass in amu", float, ut.MASSSUB_DEFAULT)
    ZSUB = ut.type_check_input("Input substrate atom charge Z in e", int, ut.ZSUB_DEFAULT)
    SUBDENSITY = ut.type_check_input("Input substrate density in atoms/cc", float, ut.SUBDENSITY_DEFAULT)
    SUBDENSITY = SUBDENSITY/1E24
    SUBWINDOW = ut.type_check_input("Input substrate window depth in Angstroms", int, ut.SUBWINDOW_DEFAULT)
    SIMULS = ut.type_check_input("How many simulations do you want?", int, ut.SIMULS_DEFAULT)

    if ut.DEBUG:
        print("Substrate Mass:", MASSSUB)
        print("Substrate Charge:", ZSUB)
        print("Substrate Density:", SUBDENSITY)
        print("Depth:", SUBWINDOW)
        print("Simulations:", SIMULS)

    # More declarations are made in the code at this point, not necessary in python