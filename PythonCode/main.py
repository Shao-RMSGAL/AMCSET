## ===============================
## AUTHOR           : Lin Shao
## CREATE DATE      :
## PURPOSE          : Collision Simulation
## SPECIAL NOTES    : Currently under the process of translation
## ===============================
## Change History:
##  Python translation currently in progress.
##==================================

from random import random
from inputs import *
from SubRoutineFunctions import *

##################INPUT_STARTING_DATA##################

massp,ZProj,MASSSUB,ZSUB,SUBDENSITY,INELAB,SUBWINDOW,SIMULS = retrieveInputs()

#########################DEFINING_LISTS#########################

# Primary Ion Range
PR = []

# Recoiled Atom Range
RR = []

# Vacancy Distribution
V = []

# Interstitial Distribution
I = []

# X Position
XP = []

# Y Position
YP = []

# Z Position
ZP = []

# Charge Number of Primary Ion
CP = []

# Mass of Primary Ion
mp = []

# Angle with Z Axis
THETAP = []

# Angle with X Axis
ALPHAP = []

# Energy of Primary Ion
EP = []

# Minimum Energy for Atom Displacement (20 eV)
EBARRIER = 0.02

###################SETTING_UP_STARTING_PROJECTILE###################

XP.append(0)
YP.append(0)
ZP.append(random() * SUBDENSITY**(1/3))
CP.append(ZProj)
mp.append(massp)
THETAP.append(0)
ALPHAP.append(0)

# INELAB is the implantation Energy of the incident ion
EP.append(INELAB)

# rum is the index. In the original code it began at 1 due to 
# how arrays used to work in Quick Basic
rum = 0
MASS2 = MASSSUB
Z2 = ZSUB
DENSITY = SUBDENSITY
