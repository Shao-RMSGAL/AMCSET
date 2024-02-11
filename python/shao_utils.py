"""
File: shao_utils.py
Author:         Nathaniel Thomas
Date created:   February 9th, 2024
Last modified:  February 9th, 2024
Group:          Shao Group, Texas A&M Nuclear Engineering

Description:
-   This file contains useful constants and functions for use in Dr.
    Shao's particle collision simulation software
"""

from typing import TypeVar, Union
import scipy as sc

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


def AMAGIC(THETAO, ALPHAO, THETA1RELATIVE, THETA2RELATIVE):
    THETA1 = -1
    ALPHA1 = -1
    THETA2 = -1
    ALPHA2 = -1
    return THETA1, ALPHA1, THETA2, ALPHA2


def BremsEloss(ZSUB, DENSITY, ElecE):
    return -1


def DF(X, COLUMBIAVK, Z1, Z2, AU):
    result = -1
    return result


def DSCMOTT(mott, scr, ElectronEnergy, SubstrateZ, Theta):
    DiffCross = -1
    return DiffCross


def EMAGIC(X):
    return -1


def F(X, COLUMBIAVK, Z1, Z2, AU):
    return -1


def IMAGE(Scan, NewDeltaZ, NewDeltaX, NewDeltaY, NewTheta, NewAlpha, NewEP,
          Resolution, size):
    return -1


def IoniElecLoss(ZSUB, DENSITY, EleE):
    return -1


def TMAGIC(MASS1, MASS2, Z2, INELAB, EP, P):
    return -1, -1, -1   


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
