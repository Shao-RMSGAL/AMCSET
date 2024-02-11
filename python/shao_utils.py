"""
File: shao_utils.py
Author:         Nathaniel Thomas
Date created:   February 9th, 2024
Last modified:  February 9th, 2024
Group:          Shao Group, Texas A&M Nuclear Engineering

Description:
-   This file contains useful constants and functions for use in Dr.
    Shao#s particle collision simulation software
"""

from typing import TypeVar, Union
import scipy as sc
import math
import random

# Constants

DEBUG: bool = True
"""Debug boolean. If True, then enable debug mode. Otherwise, disable
debug mode.
"""

DOSE: float = 0
"""Represents dose"""

ELECCP: int = -1
"""This is the charge of an electron in elementary charge units"""

ELECMASSP: float = sc.constants.electron_mass/sc.constants.atomic_mass
"""This is the mass of an electron in amu."""

# TODO: See what this value represents, and if there is a more robust way of
# implementing it.
EPRBIN: float = 1.5E6
"""This is the depth profile interval for an electron in Angstroms."""

sizeScan = 100
"""Represents the size of a scan used in electron back-scattering plot"""

WINDOWX: int = 630
"""Screen pixel width"""

WINDOWY: int = 330
"""Screen pixel height"""


# Input default values

ELECE0_DEFAULT: float = float(10000)
"""This is the default value for electron energy in keV."""

INELAB_DEFAULT: float = float(50)
"""This is the default ion incident energy in keV."""

IS_ELECTRON_DEFAULT: bool = False
"""This represents the default value for ion type. False means that the
ion isnot an electron, while True means that the ion is an electron.
"""

MASSP_DEFAULT: float = 27.9
"""This is the default value for ion mass in amu. Default mass is that
of a silicon ion.
"""

MASSSUB_DEFAULT: float = 55.85
"""This is the default value for substrate mass in amu. Default mass is
that of an iron atom.
"""

SIMULS_DEFAULT: int = 1000
"""This is the default number of simulations to carry out. Each
simulation is that of a single particle.
"""

SUBDENSITY_DEFAULT: float = 8.382E22
"""This is the default value for the density of the substrate in
atoms/cubic centimeters. Default density is that of iron.
"""

SUBWINDOW_DEFAULT: int = 80000000
"""This is the default value for the depth of the window in Angstroms."""

ZP1_DEFAULT: int = 14
"""This is the default charge number of the incident atom. Default
charge is that of silicon.
"""

ZSUB_DEFAULT: int = 26
"""This is the default charge number of the substrate atom. Default
charge is that of iron.
"""

# Functions
def sgn(x):
    if x > 0:
        return 1
    elif x == 0:
        return 0
    else:
        return -1

def AMAGIC(THETA0, ALPHA0, THETA1RELATIVE, THETA2RELATIVE):
    # RESULT TO MAIN PROGRAM, DEFLECTION ANGLE TO ORIGINAL C.
    THETA1 = -1
    # RESULT TO MAIN PROGRAM
    ALPHA1 = -1
    THETA2 = -1
    ALPHA2 = -1

    ALPHA1RELATIVE = random.random()*2*math.pi
    ALPHA2RELATIVE = ALPHA1RELATIVE+math.pi
    
    """
    ====================================================================
    for PRIMARY ION
    ====================================================================
    """
    X1 = math.sin(THETA1RELATIVE)*math.cos(ALPHA1RELATIVE)
    Y1 = math.sin(THETA1RELATIVE)*math.sin(ALPHA1RELATIVE)
    Z1 = math.cos(THETA1RELATIVE)

    Y0 = Y1*math.cos(THETA0)+Z1*math.sin(THETA0)
    Z0 = -Y1*math.sin(THETA0)+Z1*math.cos(THETA0)
    X0 = X1

    Z = Z0
    X = X0*math.sin(ALPHA0)+Y0*math.cos(ALPHA0)
    Y = -X0*math.cos(ALPHA0)+Y0*math.sin(ALPHA0)

    if Z > 0:
        THETA1 = math.atan(math.sqrt(X**2+Y**2)/Z)
    if Z == 0:
        THETA1 = math.pi / 2
    if Z < 0:
        THETA1 = math.pi + math.atan(math.sqrt(X ** 2 + Y ** 2) / Z)

    if math.sin(THETA1) != 0 :
        if X > 0:
            ALPHA1 = math.atan(Y / X)
        if X == 0:
            ALPHA1 = math.pi - sgn(Y) * math.pi / 2
        if X < 0:
            ALPHA1 = math.pi + math.atan(Y / X)
    else:
        ALPHA1 = 0

    """
    ====================================================================
    for RECOIL TARGET ION
    ====================================================================
    """
    
    X1 = math.sin(THETA2RELATIVE) * math.cos(ALPHA2RELATIVE)
    Y1 = math.sin(THETA2RELATIVE) * math.sin(ALPHA2RELATIVE)
    Z1 = math.cos(THETA2RELATIVE)

    Y0 = Y1 * math.cos(THETA0) + Z1 * math.sin(THETA0)
    Z0 = -Y1 * math.sin(THETA0) + Z1 * math.cos(THETA0)
    X0 = X1

    Z = Z0
    X = X0 * math.sin(ALPHA0) + Y0 * math.cos(ALPHA0)
    Y = -X0 * math.cos(ALPHA0) + Y0 * math.sin(ALPHA0)

    if Z > 0:
        THETA2 = math.atan(math.sqrt(X ** 2 + Y ** 2) / Z)
    if Z == 0:
        THETA2 = math.pi / 2
    if Z < 0:
        THETA2 = math.pi + math.atan(math.sqrt(X ** 2 + Y ** 2) / Z)

    if math.sin(THETA2) != 0:
        if X > 0:
            ALPHA2 = math.atan(Y / X)
        if X == 0:
            ALPHA2 = math.pi - sgn(Y) * math.pi / 2
        if X < 0:
            ALPHA2 = math.pi + math.atan(Y / X)
    else:
        ALPHA2 = 0
    return THETA1, ALPHA1, THETA2, ALPHA2


def BremsEloss(ZSUB, DENSITY, ElecE):
    #   SHARED BEloss AS DOUBLE
    BEloss = 1.4*10**(-4)*DENSITY/(6.022*10**23)*ZSUB*(ZSUB+1)*(ElecE+511)*(4*math.log(2*(ElecE+511)/511)-4/3)
    #PRINT "what the help density="; DENSITY
    # PRINT "what the help ElecE="; ElecE
    #    PRINT "BEloss inside="; BEloss
    return BEloss


def DF(X, COLUMBIAVK, Z1, Z2, AU):
    result = F(X,COLUMBIAVK,Z1,Z2,AU)/X+COLUMBIAVK/X*(.35*math.exp(-.3/X/AU)*.3/AU+.55*math.exp(-1.2/X/AU)*1.2/AU+.1*math.exp(-6/X/AU)*6/AU)
    return result


def DSCMOTT(mott, scr, ElectronEnergy, SubstrateZ, Theta):
    DiffCross = -1
    PI = 3.1415926#
    beta = (1-(ElectronEnergy/511+1)**(-2))**0.5
    betaAVE = 0.7181287
    gamma = ElectronEnergy/511+1
    elecmass = 1 #specific for a.u. unit,  only for Feq calcu.
    Vc = 137 #speed of light for a.u. unit, only for Feq calcu.
    re1 = 2.817938E-13 #in the unit of cm)

    Rmott = 0
    for xx in range(1 , 5):
        alpha1 = 0
        for yy in range(1 , 6):
            alpha1 = alpha1+mott[SubstrateZ,xx,yy]*(beta-betaAVE)**(yy-1)
        Rmott = Rmott+alpha1*(1-math.cos(Theta))**((xx-1)/2)
    DSC = (SubstrateZ*re1)**2*(1-beta**2)/beta**4*(1-math.cos(Theta))**(-2) #Differential Ruth cross section
    moment = 2*beta*gamma*elecmass*Vc*math.sin(Theta/2)
    Feq = 0
    for ww in range(1 , 3):
        Feq = Feq+scr[SubstrateZ,ww]*scr[SubstrateZ,ww+3]**2/(scr[SubstrateZ,ww+3]**2+moment**2)
    DiffCross = Rmott*DSC*(1-Feq)**2

    #PRINT "Rmott="; Rmott
    #PRINT "DSC="; DSC
    #PRINT "Feq="; Feq
    #PRINT "moment="; moment
    #PRINT "SubstrateZ="; SubstrateZ
    #PRINT "scr(SubstrateZ,1)="; scr(SubstrateZ, 1)
    #  PRINT "DiffCross="; DiffCross

    #INPUT XXXXXX
    return DiffCross


def EMAGIC(X):
    # TODO: Figure out what ATOMD is, as well as what this function does
    ATOMD = 0
    # X IS THE CURVE FITTING FROM TRIM and UNIT IS EV, X UNIT IS KEV
    if X*1000>1 and X*1000<=10:
        ELOSS = ATOMD * (.7122 + .1026 * X * 1000) / 1000
    if X * 1000 > 10 and X * 1000 <= 100 :
        ELOSS = ATOMD * (1.671 + .0252 * X * 1000) / 1000
    if X * 1000 > 100 and X * 1000 <= 400 :
        ELOSS = ATOMD * (5.0767 + .0044 * X * 1000) / 1000
    if X * 1000 > 1000 and X * 1000 <= 10000 :
        ELOSS = ATOMD * (9.28 + .00144 * X * 1000) / 1000
    if X * 1000 > 10000 :
        ELOSS = ATOMD * (9.28 + .00144 * X * 1000) / 1000
    return ATOMD, ELOSS


def F(X, COLUMBIAVK, Z1, Z2, AU):
    #===========================UNIVERAIAL SCREENING POTENTIAL=================
    if X == 0 :
        result = 0
    else:
        result = COLUMBIAVK * X * (.35 * math.exp(-.3 / X / AU) + .55 * math.exp(-1.2 / X / AU) + .1 * math.exp(-6 / X / AU))
    #=========================COLUMBIA POTENTIAL==============================
    #F = COLUMBIAVK * X
    return result


def IMAGE(Scan, DeltaZima, DeltaXima, DeltaYima, Thetaima, Alphaima, Eima,
          Resolution, size):
    m = 1

    for i in range(0 , size):
        for j in range(0 , size):
            for k in range(1 , Resolution):
                Scan[i, j, 1] = Scan[i, j, 1] + DeltaXima / Resolution
                Scan[i, j, 2] = Scan[i, j, 2] + DeltaYima / Resolution
                Scan[i, j, 3] = Scan[i, j, 3] + DeltaZima / Resolution
                999

                while (Scan[i, j, 1] ** 2 + Scan[i, j, 2] ** 2 + (Scan[i, j, 3] - 200) ** 2) ** 0.5 < 50 :
                    Scan[i, j, 1] = Scan[i, j, 1] + DeltaXima / Resolution
                    Scan[i, j, 2] = Scan[i, j, 2] + DeltaYima / Resolution
                    Scan[i, j, 3] = Scan[i, j, 3] + DeltaZima / Resolution
            Scan[i, j, 4] = Thetaima
            Scan[i, j, 5] = Alphaima
            Scan[i, j, 6] = Eima


            if Scan[i, j, 3] < 0 and Scan[i, j, 7] == 0:
                Scan[i, j, 8] = Scan[i, j, 8] + 1
                Scan[i, j, 7] = 1
    return -1


def IoniElecLoss(ZSUB, DENSITY, ElecE):
    IoniEloss1 = -1
    IoniEloss2 = -1
    # SHARED IoniEloss2 AS DOUBLE
    #  SHARED IoniEloss1 AS DOUBLE

    ##  PRINT "ElecE="; ElecE

    beta = (1 - (ElecE / 511 + 1) ** (-2)) ** 0.5
    KforLoss = 0.734 * ZSUB ** 0.037
    if ZSUB < 13 :
        JforLoss = 0.0115 * ZSUB
    if ZSUB >= 13 :
        JforLoss=0.00976*ZSUB+0.0585*ZSUB**(-0.19)

    JforLoss = JforLoss / (1 + KforLoss * JforLoss / ElecE)

    IoniEloss2 = 153.55*DENSITY*ZSUB/(6.022*10**23)/beta**2*(math.log(511*beta**2*ElecE/(2*JforLoss**2*(1-beta**2)))-math.log(2)*(2*(1-beta**2)**0.5-1+beta**2)+1-beta**2+1/8*(1-(1-beta**2)**0.5)**2)

    #  JforLoss = JforLoss / (1 + KforLoss * JforLoss / ElecE)
    #  IoniEloss1 = 7.85 * (10 ** 4) * DENSITY * ZSUB / (6.022 * 10 ** 23) / ElecE * math.log(1.166 * ElecE / JforLoss)
    #  if IoniEloss2 < IoniEloss1 : IoniEloss2 = IoniEloss1

    return IoniEloss1, IoniEloss2


def TMAGIC(MASS1, Z1, MASS2, Z2, INELAB, P):
    #=============================INITIAL PARAMETER============================

    COLUMBIAVK = .0143992 # * Z1 * Z2 #POTENTIAL V=COLUMBIAVK(ANSTRON*KEV)*Z1*Z2/R #
    MC = MASS1 * MASS2 / (MASS1 + MASS2) #MC--REDUCED MASS IN CENTER-MASS COORDINATOR
    INVLAB = math.sqrt(INELAB * 2 / MASS1) #INVLAB--INCIDENT VOLOCITY IN LAB C.
    EC = 1 / 2 * MC * INVLAB ** 2 #EC--INITIAL ENERGY IN CM
    AU = .8854 * .529 / (Z1 ** .5 + Z2 ** .5) ** (2 / 3)
    ELINHARD = EC * AU / COLUMBIAVK

    #==============FIND RMIN for DifFERENT ENERGY======================
    AA = P ** 2
    if AA == 0: 
        AA = .00001
    BB = COLUMBIAVK / EC
    CC = -1
    COLUMRMIN = 1 / 2 / AA * (-BB + math.sqrt(BB ** 2 - 4 * AA * CC)) #COLUMRMIN IS 1/RMIN, RMIN IS THE MINIMUM R UNDER COLUMBIA POTENTIAL

    CALTIME = 1
    while True:
        RMINTRY1 = COLUMRMIN

        DV = abs(-DF(RMINTRY1, COLUMBIAVK, Z1, Z2, AU) / EC - 2 * P ** 2 * RMINTRY1)
        if abs(DV) < .000001:
            DV = .1
        RMINTRY2 = RMINTRY1 + (1 - F(RMINTRY1, COLUMBIAVK, Z1, Z2, AU) / EC - P ** 2 * RMINTRY1 ** 2) / DV
        COLUMRMIN = RMINTRY2
        CALTIME = CALTIME + 1
        if CALTIME > 10000:
            break
        if abs(RMINTRY2 - RMINTRY1) > .00001 :
            continue
    
    RMIN = (RMINTRY2 + RMINTRY1) / 2
    #=====================CALCULATE DEFELCTION ANGLE=========================
    RBIERSACK = 2 * (EC - F(RMIN, COLUMBIAVK, Z1, Z2, AU)) / RMIN ** 2 / DF(RMIN, COLUMBIAVK, Z1, Z2, AU)

    BBIERSACK = P / AU
    ROBIERSACK = 1 / (RMIN * AU)
    RCBIERSACK = RBIERSACK / AU

    C1BIERSACK = .6743
    C2BIERSACK = .009611
    C3BIERSACK = .005175
    C4BIERSACK = 10
    C5BIERSACK = 6.314
    ALTHABIERSACK = 1 + C1BIERSACK * ELINHARD ** (-1 / 2)
    BELTABIERSACK = (C2BIERSACK + ELINHARD ** (1 / 2)) / (C3BIERSACK + ELINHARD ** (1 / 2))
    GAMABIERSACK = (C4BIERSACK + ELINHARD) / (C5BIERSACK + ELINHARD)
    ABIERSACK = 2 * ALTHABIERSACK * ELINHARD * BBIERSACK ** BELTABIERSACK
    GBIERSACK = GAMABIERSACK * 1 / ((1 + ABIERSACK ** 2) ** (1 / 2) - ABIERSACK)
    DELTABIERSACK = ABIERSACK * (ROBIERSACK - BBIERSACK) / (1 + GBIERSACK)
    # PRINT "COS(H/2)="; (BBIERSACK + RCBIERSACK + DELTABIERSACK) / (ROBIERSACK + RCBIERSACK)
    # END
    if P == 0:
        CALPHA1 = math.pi
    else:
        CALPHA1 = 2 * math.atan(math.sqrt((ROBIERSACK + RCBIERSACK) ** 2 / (BBIERSACK + RCBIERSACK + DELTABIERSACK) ** 2 - 1))
    # ===================DEFLECTION ANGLE IN LAB COORDINATOR===========
    COSPLUSMASS = math.cos(CALPHA1) + MASS1 / MASS2
    if COSPLUSMASS == 0:
         THETA1RELATIVE = math.pi/2
    if COSPLUSMASS < 0:
         THETA1RELATIVE = math.pi + math.atan(math.sin(CALPHA1) / COSPLUSMASS)
    if COSPLUSMASS > 0:
         THETA1RELATIVE = math.atan(math.sin(CALPHA1) / COSPLUSMASS)

    #===================CALCUTE T & DIRECTION OF TARGET ATOM=============
    RE = 4 * EC * MC / MASS2 * math.sin(CALPHA1 / 2) * math.sin(CALPHA1 / 2)
    #RE IS ENERGY TRANSIMITED TO TARGET ATOM'
    #===================RECOILED DIRECTION==================================
    # RECOILED DIRECTION
    THETA2RELATIVE = (math.pi - CALPHA1) / 2
    #PRINT "RE="; RE, "CALPHA1="; CALPHA1
    return THETA1RELATIVE, THETA2RELATIVE, RE


T = TypeVar('T')


def type_check_input(prompt: str, type: Union[type, T], default: T) -> T:
    """
    Performs type-checking for simple inputs. Inputs can either be a
    boolean, a integer, or float

    Args:
        prompt (str):   Prompt shown to the user
        type (str):     The datatype that is being input
        default (any):  The default value
    """
    if not isinstance(default, type):
        raise ValueError("Default type does not match specified type")
    if type is bool:
        try:
            mod_prompt = prompt + " (y/n - default = " + (" y" if default else
                                                          "n") + "): "
            output = input(mod_prompt)
            if output == 'y':
                output = True
            elif output == 'n':
                output = False
            elif output == "":
                output = default
            else:
                raise ValueError("Invalid option: " + output)
        except ValueError as e:
            print("Error:", e)
            print("Please enter either y or n")
            return type_check_input(prompt, type, default)
    elif type is int:
        try:
            mod_prompt = prompt + " (integer, default = " + str(default)\
                                + "): "
            output = input(mod_prompt)
            if output == "":
                output = default
            else:
                output = int(output)
        except ValueError as e:
            print("Error:", e)
            print("Please enter an integer")
            return type_check_input(prompt, type, default)
    elif type is float:
        try:
            mod_prompt = prompt + " (floating point number, default = " + \
                        str(default) + "): "
            output = input(mod_prompt)
            if output == "":
                output = default
            else:
                output = float(output)
        except ValueError as e:
            print("Error:", e)
            print("Please enter a floating-point number")
            return type_check_input(prompt, type, default)
    else:
        raise ValueError("Invalid type")
    return output


if __name__ == "__main__":
    """Only run this module as main for debugging purposes"""

    print("NOT THE MAIN FUNCTION, PLEASE RUN main.py!")
    print("Electron mass: ", ELECMASSP)
