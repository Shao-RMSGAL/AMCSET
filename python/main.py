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

# Standard library imports
import math
# from matplotlib import pyplot as plt
import numpy as np
import os
import random

# Local imports
import shao_utils as ut


# Depth of interest in Angstroms
# e_range = 80000000

# Interval number for the region of interest
e_interval = 30

# Range to averaging neighboring data
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
flyjudge = 1

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
mott = np.zeros((96, 5, 6))

# TODO: Write code for reading mott potential values in from a text file.

# Matrix for storing screening potentials up to element 92. TODO: Change to
# read up to element 118
scr = np.zeros((92, 6))

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
        out_10 = open(os.path.join(output_directory, "20imageEd40.txt"), 'w')
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
        ZP1 = ut.type_check_input("Input incident charge Z in e", int,
                                  ut.ZP1_DEFAULT)
        INELAB = ut.type_check_input("Input incident atom energy in keV",
                                     float, ut.INELAB_DEFAULT)
        if ut.DEBUG:
            print("Mass:", massp)
            print("Charge:", ZP1)
            print("Energy:", INELAB)

    # Continue on to inputting substrate details
    MASSSUB = ut.type_check_input("Input substrate mass in amu", float,
                                  ut.MASSSUB_DEFAULT)
    ZSUB = ut.type_check_input("Input substrate atom charge Z in e", int,
                               ut.ZSUB_DEFAULT)
    SUBDENSITY = ut.type_check_input("Input substrate density in atoms/cc",
                                     float, ut.SUBDENSITY_DEFAULT)
    SUBDENSITY = SUBDENSITY/1E24
    SUBWINDOW = ut.type_check_input("Input substrate window depth in "
                                    "Angstroms", int, ut.SUBWINDOW_DEFAULT)
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
    # arbitrary size. Currently depth of collision creation is limited to 2000
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
    CP = np.zeros(array_size, dtype=int)
    # Save mass of ion/electron projectiles or target atoms
    mp = np.zeros(array_size, dtype=float)
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
    # Electron data, third index 1 is for location, 2 is ionization energy loss
    # 3 is scattering, and 4 is vacancy production
    e_output = np.zeros((e_interval, e_interval, 4))
    RL = np.zeros(flyC + 1)
    RN = np.zeros(flyC)

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
    random.seed(None)

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
            CP[1] = ZP1
            mp[1] = massp
            EP[1] = INELAB

        """
        ================================================================
        SCANNING
        ================================================================
        """
        if is_electron:
            for i in range(1, ut.sizeScan):
                for j in range(1, ut.sizeScan):
                    Scan[i, j, 1] = scandelta*(i-ut.sizeScan/2)
                    Scan[i, j, 2] = scandelta*(i-ut.sizeScan/2)
                    Scan[i, j, 6] = ElecE0
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
        # Substrate density
        DENSITY = SUBDENSITY

        # LABEL: 100
        # To simulate a GOTO 600 statement, this will skip everything
        # that is before label 600.
        # TODO: This is recursive behavior, replace with recursive function
        while rum != 0:
            # If the particle flies out of the window, stop the
            # simulation and go to the next collision
            if ZP[rum] >= SUBWINDOW:
                rum = rum - 1
                # We want to go back to LABEL 100, so we continue to the
                # next loop iteration
                continue
            # Back-scattered electron
            # TODO: Replace with enumeration
            elif ZP[rum] < 0:
                if Ptype[rum] == 0:
                    SPUTTER = SPUTTER+1
                    rum -= 1
                    continue
            if Ptype[rum] == 1:
                # Electron-specific energy checks
                if EP[rum] < e_stopping:
                    # Compute radial distance from the z axis
                    # TODO: Pick more descriptive intermediate variable names,
                    E0 = math.sqrt(XP[rum]**2 + YP[rum]**2)
                    E1 = ZP[rum]
                    # Check that electrons are in a valid region
                    if E0 < SUBWINDOW and E1 > 0 and E1 < SUBWINDOW:
                        lat = round(1+E0/e_deltaLat)
                        depth = round(1+E1/e_deltaDep)
                        e_output[lat, depth, 1] += 1
                    # Electron has stopped, move on.
                    rum -= 1
                    continue
                elif EP[rum] < flyjudge * e_stopping:
                    # Prevents negative electron energy
                    # TODO: The QB code claims that flyjudge is set in the
                    # input. Right now, it is not.
                    flyC = fly0
                elif EP[rum] <= EMAX:
                    # Allows for normal flying distances
                    flyC = fly1
                elif EP[rum] > EMAX:
                    # For electrons with excess energy, avoid too-high
                    # flying distances
                    flyC = fly0
                # At this point, the QB code has a separate if-statement for
                # more electron bombardment, it will be continued here instead
                for trial2 in range(1, flyC):
                    # Rank random number from low to high
                    RL[trial2] = random.random()
                RL = np.sort(RL)

                TotalCross = 0
                for ANGLE1 in range(1, Divisor-1):
                    # Divide 180 degrees by Divisor, integrating
                    # Obtain differential cross section at specific angle
                    DiffCross = ut.DSCMOTT(mott, scr, EP[rum]*correctFactor,
                                           ZSUB,
                                           math.pi/(1 - 10)
                                           * (1 - 10**(ANGLE1/Divisor)))
                    if ANGLE1 == 1:
                        # Integration of the first angle, starting at 0
                        TotalCross += DiffCross*2*math.pi\
                            * math.sin(math.pi/(1-10)*(1-10**(1/Divisor)))\
                            * ((math.pi/2/(1-10)*(1 - 10**(2/Divisor)))
                                + (math.pi/2/(1-10)*(1-10**(1/Divisor))))
                    elif ANGLE1 < Divisor-1:
                        # integration for angle intervals excluding two
                        # boundary points
                        TotalCross += DiffCross*2*math.pi * \
                            math.sin(math.pi/(1-10)*(1-10**(ANGLE1/Divisor)))\
                            * ((math.pi/2/(1-10)*(1-10**((ANGLE1+1)/Divisor)))
                                - (math.pi/2/(1-10)*(1-10**((ANGLE1-1) /
                                                            Divisor))))
                    elif ANGLE1 == Divisor-1:
                        # integration concerning the last angle point
                        # concerning the boundary at 180 degree
                        TotalCross += DiffCross*2*math.pi\
                                * math.sin(math.pi
                                           / (1-10)*(1-10**((Divisor-1)
                                                            / Divisor)))\
                                * (math.pi-(math.pi/2/(1-10)*(1-10
                                                              ** ((Divisor-1)
                                                                  / Divisor)))
                                    - (math.pi/2/(1-10)*(1-10**((Divisor-2)
                                                                / Divisor))))
                # Prepare for flying distance assessment
                Lselected = 0
                # Start flying distance assessment. flyC is the group size
                for trial1 in range(1, flyC):
                    # TODO: Alter the range of random() to exclude 0
                    ForLselect = random.random()
                    if ForLselect == 0:
                        ForLselect = 1E-10
                    # Calculate flying distance using TotalCross
                    RN[trial1] = -math.log(ForLselect)\
                        / (SUBDENSITY*1E24*TotalCross)+(SUBDENSITY * 1E24)\
                        ** (-1/3)
                    # Add all free flying distances together
                    Lselected += RN[trial1]
                    # Prepare to identify scattering angle
                    ForThetaSelect0 = 0
            else:
                # Account for the case of an ion
                # TODO: This code is reused. This can be refactored.
                if EP[rum] < ion_stopping:
                    # Calculate radial distance from the z axis
                    E0 = math.sqrt(XP[rum]**2 + YP[rum]**2)
                    # Longitudinal depth at stopping position
                    E1 = ZP[rum]
                    if E0 < SUBWINDOW and E1 > 0 and E1 < SUBWINDOW:
                        lat = round(1 + E0/e_deltaLat)
                        long = round(1 + E1/e_deltaDep)
                        e_output[lat, long, 4] += 1
                    rum -= 1
                    continue

                # Electron stopping parameter
                KL = 1.212*CP[rum]**(7/6)*Z2/(CP[rum]**(2/3)+Z2**(2/3))**(3/2)\
                    / math.sqrt(mp[rum])
                # Electron stopping energy at EP, units of EV
                SE = KL*math.sqrt(EP[rum]*1000)
                L = DENSITY**(-1/3)
                # Using random number generator to select collision parameter
                # TODO: Check that the QuickBasic RND function is between 1-0
                P = math.sqrt(random.random()/(math.pi*DENSITY**(2/3)))
                """
                ========================================================
                COLLISION PARAMETER
                ========================================================
                """
                # Get updated energy after electron energy loss
                EP[rum] = EP[rum] - 1.59*L*DENSITY*SE/1000

                if EP[rum] < ion_stopping:
                    # Check again after electronic losses to see if ion stops
                    rum -= 1
                    continue

                # Call the TMAGIC function to get recoil energy
                # TODO: Understand exactly what TMAGIC and AMAGIC d0
                (THETA1RELATIVE, THETA2RELATIVE, RE) = ut.TMAGIC(mp[rum],
                                                                 CP[rum],
                                                                 MASS2, Z2,
                                                                 EP[rum], P)
                # Obtain new direction and angles
                (THETA1, ALPHA1, THETA2, ALPHA2) = ut.AMAGIC(THETAP[rum],
                                                             ALPHAP[rum],
                                                             THETA1RELATIVE,
                                                             THETA2RELATIVE)
                """
                ========================================================
                TEMPORARY SAVE INFORMATION
                TODO: Understand what exactly this section is doing
                ========================================================
                """
                # Assign new depth to target
                ZP[rum+1] = ZP[rum]+L*math.cos(THETAP[rum])
                XP[rum+1] = XP[rum]+L*math.sin(THETAP[rum])*math.cos(
                                                            ALPHAP[rum])
                YP[rum+1] = YP[rum]+L*math.sin(THETAP[rum])*math.sin(
                                                            ALPHAP[rum])
                # Transfer charge information to target.
                CP[rum+1] = Z2
                # Transfer target mass information
                mp[rum+1] = MASS2
                # Transfer angle information of recoiled target
                THETAP[rum+1] = THETA2
                # Transfer angle information of recoiled target
                ALPHAP[rum+1] = ALPHA2
                # Transfer recoil energy to target atom
                EP[rum+1] = RE

                # Plotting
                # TODO: Implement plotting
                if ut.DEBUG:
                    print("Simulation: ", SIMUL)

                if mp[rum] == massp:
                    # TODO: Plot the relevant line
                    if ut.DEBUG:
                        print("Line drawn")
                if mp[rum] == MASSSUB:
                    # TODO: Plot the relevant line
                    if ut.DEBUG:
                        print("Alt line drawn")
                if RE <= ion_Ed:
                    # If target atom energy is too low, there's no displacement
                    # Return coordinates back to projectile as an update
                    ZP[rum] = ZP[rum+1]
                    XP[rum] = XP[rum+1]
                    YP[rum] = YP[rum+1]
                    # Update new direction, from AMAGIC
                    THETAP[rum] = THETA1
                    ALPHAP[rum] = ALPHA1
                    # Update energy after recoil loss
                    EP[rum] -= RE
                    continue

                # If RE > 0.02, a target displacement occurred.
                # Projectile has the same coordinate data as the target
                XP[rum] = XP[rum+1]
                YP[rum] = YP[rum+1]
                ZP[rum] = ZP[rum+1]
                # Update angle obtained from AMAGIC. '1' is for the projectile
                THETAP[rum] = THETA1
                ALPHAP[rum] = ALPHA1
                # Subtract energy loss
                EP[rum] -= RE

                # TODO: Replace with redundancy code
                if ZP[rum+1] > 0 and ZP[rum+1] < SUBWINDOW:
                    # Due to displacement, save vacancy information
                    # Distance away from z axis
                    E0 = math.sqrt(XP[rum]**2 + YP[rum]**2)
                    # Depth along z axis
                    E1 = ZP[rum]
                    long = round(1 + E0/e_deltaLat)
                    lat = round(1 + E1/e_deltaDep)
                    e_output[long, lat, 4] += 1
                # Move on to next particle
                rum += 1
                break
