"""
File: main.py
Author:         Nathaniel Thomas
Date created:   February 9th, 2024
Last modified:  February 9th, 2024
Group:          Shao Group, Texas A&M Nuclear Engineering

Description:
-   This is the primary file for Dr. Shao's particle collision simulation 
    software. TODO: Provide a more detailed description
"""

import numpy as np
import shao_utils as ut
import os
import random as rd

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

Resolution = 10
scandelta = 20
EMAX = 10000

# Switch value. '0' is for an explicit method, '1' is for the middle-point
# method, and '2' is for back (implicit) method
switch = 1

# The code has several type declaration statements at this point. Python
# does not require this.

Ltotal = 0.0

# Matrix for storing Mott potentials (TODO: Include all elements up to 118)
mott = np.zeros((96,5,6))

# TODO: Write code for reading mott potential values in from a text file.

# Matrix for storing screening potentials up to element 92. TODO: Change to
# read up to element 118
scr = np.zeros((92,6))

# TODO: Write code for reading screening potential values

if __name__ == "__main__":
    current_directory = os.getcwd()

    if os.path.basename(current_directory) != "python":
        try:
            new_directory = os.path.join(current_directory, "python")
            os.chdir(new_directory)
        except Exception as e:
            print("Error: ", e)
            exit()



    # Various output files
    try:
        current_directory = os.getcwd()
        output_directory = os.path.join(current_directory, "output")
        out_1 = open(os.path.join(output_directory, "PRANGEE40.dat"), 'w')
        out_2 = open(os.path.join(output_directory, "simulEd40.dat"), 'w')
        out_3 = open(os.path.join(output_directory,
                                  "simulvimageEd40.dat"), 'w')
        out_4 = open(os.path.join(output_directory, "DINTLimageEd40.dat"), 'w')
        out_5 = open(os.path.join(output_directory,
                                  "eRangeimageEd40.dat"), 'w')
        out_6 = open(os.path.join(output_directory,
                                  "ProjectedeRangeimageEd40.dat"), 'w')
        out_7 = open(os.path.join(output_directory, "e-outimageEd40.dat"), 'w')
        out_9 = open(os.path.join(output_directory,
                                  "Projected_e_rangeImageEd40.txt"), 'w')
        out_10 = open(os.path.join(output_directory,"20imageEd40.txt"), 'w')
        out_11 = open(os.path.join(output_directory,
                                   "crossimageEd40.txt"), 'w')
        out_12 = open(os.path.join(output_directory, "imageEd40.txt"), 'w')
        out_13 = open(os.path.join(output_directory, "VacancyEd40.txt"), 'w')
        out_14 = open(os.path.join(output_directory, "dedx.txt"), 'w')
    except Exception as e:
        print("Error:", e)
        ignore_error = ut.type_check_input("Continue despite error?", bool,
                                           False)
        if not ignore_error:
            exit()

    # iontype determines whether it is an electron bombardment or ion
    # bombardment. The default is 0
    is_electron = ut.type_check_input("Is this an electron bombardment?", bool,
                                      ut.IS_ELECTRON_DEFAULT)

    if ut.DEBUG:
        print("Is electron:", is_electron)
    # In the case of electron bombardment, ask for energy:
    if is_electron:
        # Ask for the electron input if electron bombardment is being simulated
        ElecE0 = ut.type_check_input("Input electron energy in keV", float,
                                     ut.ELECE0_DEFAULT)
        if ut.DEBUG:
            print("Energy:", ElecE0)
    else:
        # Ask for ion information if not simulating electron bombardment
        massp = ut.type_check_input("Input incident atom mass in amu", float,
                                    ut.MASSP_DEFAULT)
        ZP = ut.type_check_input("Input incident charge Z in e", int,
                                 ut.ZP_DEFAULT)
        INELAB = ut.type_check_input("Input incident atom energy in keV",
                                     float, ut.INELAB_DEFAULT)
        if ut.DEBUG:
            print("Mass:", massp)
            print("Charge:", ZP)
            print("Energy:", INELAB)

    # Continue on to inputting substrate details
    MASSSUB = ut.type_check_input("Input substrate mass in amu", float,
                                  ut.MASSSUB_DEFAULT)
    ZSUB = ut.type_check_input("Input substrate atom charge Z in e", int,
                               ut.ZSUB_DEFAULT)
    SUBDENSITY = ut.type_check_input("Input substrate density in atoms/cc",
                                     float, ut.SUBDENSITY_DEFAULT)
    SUBDENSITY = SUBDENSITY/1E24
    SUBWINDOW = ut.type_check_input("Input substrate window depth in Angstroms"
                                    , int, ut.SUBWINDOW_DEFAULT)
    SIMULS = ut.type_check_input("How many simulations do you want?", int,
                                 ut.SIMULS_DEFAULT)

    if ut.DEBUG:
        print("Substrate Mass:", MASSSUB)
        print("Substrate Charge:", ZSUB)
        print("Substrate Density:", SUBDENSITY)
        print("Depth:", SUBWINDOW)
        print("Simulations:", SIMULS)

    # More declarations are made in the code at this point, not
    #    necessary in python
    
    # Initialize arrays TODO: Change this to allow for arrays of
    # arbitrary size
    array_size: int = 2000

    # Represents x position of ion/electron projectiles or target atoms,
    # up to 2000 in a single bombardment
    XP = np.zeros(array_size)
    # Represents y position of ion/electron projectiles or target atoms,
    # up to 2000 in a single bombardment
    YP = np.zeros(array_size)
    # Represents z position of ion/electron projectiles or target atoms,
    # up to 2000 in a single bombardment
    ZP = np.zeros(array_size)
    # Save charge number of `ion/electron projectiles or target atoms
    CP = np.zeros(array_size, dtype = int)
    # Save mass of ion/electron projectiles or target atoms
    mp = np.zeros(array_size, dtype = float)
    # Save angle with respect to the Z axis
    THETAP = np.zeros(array_size)
    # Save angle with respect to the X axis
    ALPHAP = np.zeros(array_size)
    # Save energy of ion/electron projectile or target atoms
    EP = np.zeros(array_size)
    # Save the particle type. 0 = ion, 1 = electron
    Ptype = np.zeros(array_size)
    # Save the profile of ion/electron projectiles along the z axis.
    ePR = np.zeros(array_size)
    
    # TODO: Figure out exactly what THETA2 and ALPHA2 do
    THETA2 = np.zeros(array_size)
    ALPHA2 = np.zeros(array_size)
    
    # For creating 2D plot of back-scattered electrons
    Scan = np.zeros((ut.sizeScan + 1, ut.sizeScan + 1, 8))
    
    # The QB Code has statements here that zero all the vectors. There
    # is no need for that thanks to np.zeros()
    SPUTTER: int = 0
    """Represents the number of projectiles that sputter"""

    # This seeds the random number generator using the system time.
    rd.seed(None)

    for SIMUL in range(1, SIMULS):
        """
        QB code resets the random seed every 200 iterations. The
        Mersenne Twister random number generator (used in the python
        'random' library has a period of 2^199937 - 1, making this
        step unnecessary)
        """
        # Account for the case of an electron
        if is_electron:
            # Set electron charge
            CP[1] = ut.ELECCP
            # Set electron masss
            mp[1] = ut.ELECMASSP
            # Set electron energy
            EP[1] = ElecE0
            # The particle type is an electron
            Ptype[1] = 1
            # TODO: Figure out what the correctFactor is
            correctFactor = 1
        # Account for the case of an ion
        else:   
            CP[1] = ZP
            mp[1] = massp
            EP[1] = INELAB
        
        """
        ================================================================
        SCANNING
        ================================================================
        """

        for i in range(1, ut.sizeScan):
            for j in range(1, ut.sizeScan):
                Scan[i,j,1] = scandelta*(i - ut.sizeScan/2)
                Scan[i,j,2] = scandelta*(i - ut.sizeScan/2)
                Scan[i,j,6] = ElecE0
                # All other entires remain zero
        
        """
        ================================================================
        END SCANNING
        ================================================================
        """

        # The first simulation of a new bombardment event
        rum = 1

        # Substrate atomic mass
        MASS2 = MASSSUB

        # Substrate atomic number
        Z2 = ZSUB
