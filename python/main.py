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

# TODO: FIX INDEXING TO START FROM 0!!!

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
mott = np.zeros((96+1, 5+1, 6+1))

# TODO: Write code for reading mott potential values in from a text file.

# Matrix for storing screening potentials up to element 92. TODO: Change to
# read up to element 118
scr = np.zeros((92+1, 6+1))

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
    XP = np.zeros(array_size+1)
    # Represents y position of ion/electron projectiles or target atoms,
    # up to 2000 in a single bombardment
    YP = np.zeros(array_size+1)
    # Represents z position of ion/electron projectiles or target atoms,
    # up to 2000 in a single bombardment
    ZP = np.zeros(array_size+1)
    # Save charge number of `ion/electron projectiles or target atoms
    CP = np.zeros(array_size+1, dtype=int)
    # Save mass of ion/electron projectiles or target atoms
    mp = np.zeros(array_size+1, dtype=float)
    # Save angle with respect to the Z axis
    THETAP = np.zeros(array_size+1)
    # Save angle with respect to the X axis
    ALPHAP = np.zeros(array_size+1)
    # Save energy of ion/electron projectile or target atoms
    EP = np.zeros(array_size+1)
    # Save the particle type. 0 = ion, 1 = electron
    Ptype = np.zeros(array_size+1)
    # Save the profile of ion/electron projectiles along the z axis.
    ePR = np.zeros(array_size+1)
    # Electron data, third index 1 is for location, 2 is ionization energy loss
    # 3 is scattering, and 4 is vacancy production
    e_output = np.zeros((e_interval+1, e_interval+1, 4+1))
    Backup_e_output = np.zeros((e_interval+1, e_interval+1, 4+1))
    RL = np.zeros(flyC + 1+1)
    RN = np.zeros(flyC+1)
    RA = np.zeros(flyC+1)
    RE = np.zeros(flyC+1)

    # TODO: Figure out exactly what THETA2 and ALPHA2 do
    THETA2 = np.zeros(array_size+1)
    ALPHA2 = np.zeros(array_size+1)

    # For creating 2D plot of back-scattered electrons
    Scan = np.zeros((ut.sizeScan + 1+1, ut.sizeScan + 1+1, 8+1))

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
                    else:
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
                # Pointer needs to be updated
                trial5 = 1
                # Assign an unrealistic angle to the last scattering angle
                RL[flyC+1] = 1234
                # Go through each angle from 0 to 180
                for ANGLE2 in range(1, Divisor-1):
                    DiffCross = ut.DSCMOTT(mott, scr, EP[rum]
                                            * correctFactor,  ZSUB, math.pi
                                            / (1-10)*(1-10**(ANGLE2 /
                                                            Divisor)))
                    # Integration concerning first point at angle = 0 needs
                    # to be specially treated. 
                    if ANGLE2 == 1:
                        ForThetaSelect1 = ForThetaSelect0+DiffCross\
                            / TotalCross*2*math.pi\
                            * math.sin(math.pi/(1-10)*(1-10**(1/Divisor)))\
                            * ((math.pi/2/(1-10)*(1-10**(2/Divisor)))
                                + (math.pi/2/(1-10)*(1-10**(1/Divisor))))
                    elif ANGLE2 < Divisor - 1:
                        # integration for the middle points without two
                        # angle boundaries
                        ForThetaSelect1 = ForThetaSelect0 + DiffCross\
                            / TotalCross*2*math.pi\
                            * math.sin(math.pi/(1-10)
                                            * (1-10**(ANGLE2 / Divisor)))\
                            * ((math.pi/2/(1-10)*(1-10
                                ** ((ANGLE2+1)/Divisor)))-(math.pi/2/(1-10)
                                                            * (1 - 10**((
                                                                ANGLE2 - 1)
                                                                / Divisor)))
                                )
                    else:
                        # integration concerning the last boundary needs to
                        # be specially treated
                        ForThetaSelect1 = ForThetaSelect0+DiffCross\
                            / TotalCross*2*math.pi*math.sin(math.pi/(1-10)
                                                            * (1-10
                                                            ** ((Divisor-1)
                                                                / Divisor)
                                                                ))\
                            * (math.pi-(math.pi/2/(1-10)*(1-10**((Divisor
                                                                    - 1)
                                                                    / Divisor
                                                                    )))
                                        - (math.pi/2/(1-10)*(1-10
                                                            ** ((Divisor-2)
                                                                / Divisor))
                                            ))
                    # If the integration cross section is larger than the
                    # first random number, do the following
                    while ForThetaSelect1 > RL[trial5]:
                        if ANGLE2 == 1:
                            # ‘if it happens to the first angle interval,
                            # special treatment needs since interval width
                            # differs from middle points. RA is the
                            # scattering angle selected
                            RA[trial5] = math.pi/2/(1-10)*(1-10**(2/Divisor))+math.pi/2/(1-10)*(1-10**(1/Divisor))-(ForThetaSelect1-RL[trial5])*(math.pi/2/(1-10)*(1-10**(2/Divisor))+math.pi/2/(1-10)*(1-10**(1/Divisor)))/(ForThetaSelect1-ForThetaSelect0)
                        elif ANGLE2 < Divisor - 1:
                            # if it happens for the middle point,
                            # corresponding angle is read from the
                            # proportionality, judged by the distance from
                            # the right side boundary point.  RA is the
                            # scattering angle selected.
                            RA[trial5] = math.pi/2/(1-10)*(1-10**((ANGLE2+ 1)/Divisor))+math.pi/2/(1-10)*(1-10**(ANGLE2/Divisor))-(ForThetaSelect1-RL[trial5])*(math.pi/2/(1-10)*(1-10**((ANGLE2+1)/Divisor))-math.pi/2/(1-10)*(1-10^((ANGLE2-1)/Divisor)))/(ForThetaSelect1-ForThetaSelect0)
                        else:
                            # if it happens to the last angle interval,
                            # special treatment needs since the interval
                            # width differs from middle points. RA is the
                            # scattering angle selected.
                            RA[trial5] = math.pi-(ForThetaSelect1-RL[trial5])*(math.pi-math.pi/2/(1-10)*(1-10**((Divisor-1)/Divisor))-math.pi/2/(1-10)*(1-10**((Divisor-2)/Divisor)))/(ForThetaSelect1-ForThetaSelect0)
                            # moves to the next random point. The pointer
                            # is increased by 1
                            trial5 += 1
                    # For integration. The integrated value is saved as the
                    # base line, to be added with the increased value from
                    # the next angle interval, to be calculated in the next
                    # angle point.
                    ForThetaSelect0 = ForThetaSelect1
                # Go through all distances from the last one to the first one
                for trial99 in range(flyC, 2):
                    # Randomly pick a number smaller than trial99
                    jj = round(random.random*trial99)+1
                    # Save scattering angle of trial99 temporarily 
                    qqq = RA[trial99]
                    # transfer randomly picked scattering angle to trial99
                    RA[trial99] = RA[jj]
                    # transfer temporarily save value to the randomly
                    # picked free flying distance. Hence swapping finishes
                    RA[jj] = qqq
                # Prepare for energy loss calculation
                totalRE = 0

                # Go through each free flying distance in the group
                for trial6 in range(1, flyC):
                    # For calculation of energy transfer
                    ElecREup = ((EP[rum]+511)*(math.sin(RA[trial6]))**2+MASSSUB*931*1000*(1-math.cos(RA[trial6])))*EP[rum]*(EP[rum]+2*511)
                    # For calculation of energy transfer
                    ElecREdown = (EP[rum]+MASSSUB*931.5*1000)**2-EP[rum]*(EP[rum]+2*511)*(math.cos(RA[trial6]))**2
                    # For calculation of energy transfer
                    RE[trial6] = ElecREup/ElecREdown
                    # Adding all energy loss with the group
                    totalRE = totalRE+RE[trial6]
                """"
                ====================================================
                ADD ALL FREE-FLYING DISTANCES TOGETHER
                """
                # adding all free flying distances to get total flying
                # distance for the whole group 
                Ltotal += Lselected
                # Calculate energy loss due to ionization 
                (IoniEloss1, IoniEloss2) = ut.IoniElecLoss(ZSUB, SUBDENSITY*1E24, EP[rum]*correctFactor)
                # Make sure the value is positive 
                IoniEloss2 = (IoniEloss2 + abs(IoniEloss2))/2.0
                # Calculate energy loss due to braking irradiation 
                BEloss = ut.BremsEloss(ZSUB, SUBDENSITY*1E24, EP[rum]*correctFactor)
                # Make sure the value is positive
                BEloss = (BEloss+abs(BEloss))/2

                # Energy of electron after both non-Mott and Mott
                # scattering energy loss
                Energy1 = EP[rum]-Lselected*(IoniEloss2+BEloss)-totalRE

                # No energy correction is needed when switch=0
                if switch == 0:
                    NotImplemented
                # turn on energy correction is switch=1, correction below 
                # follows midpoint approximation
                elif switch == 1:
                    # ”1” means the correction was not performed yet since
                    # “1” is the preassigned value, do the following
                    if correctFactor == 1:
                        # the ratio of middle energy to the starting energy
                        # for the flying distance group
                        correctFactor = (Energy1+EP[rum])/2/EP[rum]
                        break

                    # Calculate the correction factor of the current group, 
                    # and use it for the next group
                    correctFactor = (Energy1+EP[rum])/2/EP[rum]
                # Switch == 2, turn on the correction but the energy correction follows the implicit method
                else:
                    # ”1” means the correction was not performed yet since “1” is the preassigned value, do the following
                    if correctFactor == 1: 
                        # Ratio is the final energy to the initial energy of the group. 
                        correctFactor = Energy1/EP[rum]
                        # Repeating the calculation for the first group. 
                        break
                    # Calculate the correction factor of the current group, and use it for the next group
                    correctFactor = Energy1/EP[rum]
                
                # Assign the final energy of the group as an updated energy as the starting energy of the next group
                EP[rum] = Energy1

                # Lateral distance from z axis
                E0 = math.sqrt(XP[rum]**2+YP[rum]**2)
                # Depth
                E1 = ZP[rum]                           
                # ionization energy loss rate
                E2 = Lselected*IoniEloss2
                # Mott scattering energy loss
                E3 = totalRE

                # for point with the valid region between two boundaries
                if E0 < SUBWINDOW and E1 > 0 and E1 < SUBWINDOW:
                    long = round(1+E0/e_deltaLat)
                    lat = round(1+E1/e_deltaDep)
                    # Saving ionization energy to 2-D position matrix
                    e_output[long, lat, 2] += E2
                    # saving Mott scattering energy loss to 2-D position matrix
                    e_output[long, lat, 3] += E3
                
                # for each collision after each flying distance, calculate the direction 
                for trial9 in range(1, flyC):
                    # TODO: Remember to correct indexing
                    index = rum + trial9 -1 -1
                    ZP[rum+trial9] = ZP[index] + RN[trial9]*10**8*math.cos(THETAP[index])
                    XP[rum + trial9] = XP[index] + RN[trial9]*10**8*math.sin(THETAP[index])*math.cos(ALPHAP[index])
                    YP[rum + trial9] = YP[index] + RN[trial9]*10**8*math.sin(THETAP[index])*math.sin(ALPHAP[index])
                    # calculate scattering angles after each collision within the free flying distance group
                    # scattering angle from each Mott scattering with respect to electron flying direction prior to collision 
                    THETA1RELATIVE = RA[trial9]
                    # approximation for target atom
                    THETA2RELATIVE = (math.pi-RA[trial9])/2
                    # Convert to angles with respect to xyz coordinate
                    (THETA1, ALPHA1, THETA2, ALPHA2) = ut.AMAGIC(THETAP[index], ALPHAP[index], THETA1RELATIVE, THETA2RELATIVE)

                    # assign scattering angles to projectile and target atoms for each collision
                    # for electron
                    THETAP[rum + trial9] = THETA1
                    # for electron
                    THETAP[rum + trial9] = ALPHA1
                    # for target atom
                    THETA2[rum + trial9] = THETA2
                    # for target atom
                    ALPHA2[rum + trial9] = ALPHA2

                print(" SIMULATION:", SIMUL, "(", SIMULS, ")")

                # TODO: Implement plotting

                NewDeltaZ = ZP[rum+flyC]-ZP[rum]
                NewDeltaX = XP[rum+flyC]-XP[rum]
                NewDeltaY = YP[rum+flyC]-YP[rum]
                NewTheta = THETAP[rum+flyC]
                NewAlpha = ALPHAP[rum+flyC]
                NewEP = EP[rum]

                # Update position after each flying distance group, transfer other information needed.  
                ZP[rum] = ZP[rum+flyC]
                XP[rum] = XP[rum+flyC]
                YP[rum] = YP[rum+flyC]
                THETAP[rum] = THETAP[rum+flyC]
                ALPHAP[rum] = ALPHAP[rum+flyC]

                mmmm = 1
                for trial8 in range(1, flyC):
                    if RE[trial8] >= ion_Ed:
                        ZP[rum+mmmm] = ZP[rum+trial8]
                        XP[rum+mmmm] = XP[rum+trial8]
                        YP[rum+mmmm] = YP[rum+trial8]
                        CP[rum+mmmm] = Z2
                        mp[rum+mmmm] = MASS2
                        THETAP[rum+mmmm] = THETA2[rum+trial8]
                        ALPHAP[rum+mmmm] = ALPHA2[rum+trial8]
                        Ptype[rum+mmmm] = 0
                        EP[rum+mmmm] = RE(trial8)
                        E0 = math.sqrt(XP[rum+mmmm]**2 + YP[rum + mmmm]**2)
                        E1 = ZP[rum + mmmm]
                        if E0 < SUBWINDOW and E1 > 0 and E1 < SUBWINDOW:
                            lat = round(1 + E0 / e_deltaLat)
                            dep = round(1 + E1 / e_deltaDep)
                            e_output[lat, dep, 4] += 1
                            # The above is to record one vacancy created by electron
                            E4V = EP[rum]
                            # TODO: Figure out what print statement in QB code does
                        
                        mmmm += 1
                
                rum += mmmm-1
                break

            else:
                """
                ========================================================
                ION BOMBARDMENT
                ========================================================
                """
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
                        e_output[lat-1, long-1, 4] += 1
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
    """
    ====================================================================
    SIMULATION COMPLETION
    ====================================================================
    """

    for N0 in range(1,e_interval):
        for N1 in range(1,e_interval):
            for N2 in range(1,4):
                Backup_e_output[N0, N1, N2] = e_output(N0, N1, N2)
    
    for R0 in range(1 + e_ave_range, e_interval - e_ave_range):
        for R1 in range(1 + e_ave_range, e_interval - e_ave_range):
            combine1 = 0
            combine2 = 0
            combine3 = 0
            combine4 = 0

            for R2 in range(-e_ave_range, e_ave_range):
                for R22 in range(-e_ave_range, e_ave_range):

                    combine1 = Backup_e_output[R0 + R2, R1 + R22, 1] + combine1
                    combine2 = Backup_e_output[R0 + R2, R1 + R22, 2] + combine2
                    combine3 = Backup_e_output[R0 + R2, R1 + R22, 3] + combine3
                    combine4 = Backup_e_output[R0 + R2, R1 + R22, 4] + combine4
                
            e_output[R0, R1, 1] = combine1/(2*e_ave_range+1)**2
            e_output[R0, R1, 2] = combine2/(2*e_ave_range+1)**2
            e_output[R0, R1, 3] = combine3/(2*e_ave_range+1)**2
            e_output[R0, R1, 4] = combine4/(2*e_ave_range+1)**2

    for W0 in range(1, e_interval):
        accum = 0
        for w1 in range(1, e_interval):
            accum = accum + e_output(w1, W0, 1)
        # TODO: Implement print statement

    test_e = 0
    test_IoniE = 0
    test_RecE = 0
    test_Disp = 0


    for R3 in range(1, e_interval):
        for R4 in range(1, e_interval):
            # unit is micron^3
            unit_volume = 2*math.pi*(R3+0.5)*e_deltaLat*e_deltaLat*e_deltaDep/10**12
            # TODO: Implement print statement
            test_e += unit_volume*e_output[R3, R4, 1]/unit_volume/SIMULS
            test_IoniE += e_output[R3,R4,2]/unit_volume/SIMULS*unit_volume
            test_RecE += e_output(R3, R4, 3) / unit_volume / SIMULS * unit_volume
            test_Disp += e_output(R3, R4, 4) / unit_volume / SIMULS * unit_volume
        
    print("integrated total electron inside =", test_e)
    print("integrtated total ioniztion (inelastic) energy loss=", test_IoniE)
    print("integrated total elastic energy loss=", test_RecE)
    print("integrated total displacement creation=", test_Disp)