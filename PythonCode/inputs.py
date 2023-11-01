## ===============================
## AUTHOR           : Lin Shao
## CREATE DATE      :
## PURPOSE          : Collision Simulation
## SPECIAL NOTES    : Currently under the process of translation
## ===============================
## Change History:
##  Python translation currently in progress.
##  This file is for the inputs that happen at the beginning of the program.
##==================================

def retrieveInputs():
##################INPUT_STARTING_DATA##################
# Mass of the incident atom
    massp = int(input("Input incident atom mass:") or  "14")

# Charge of incident atom
# In the original code this variable was called "ZP" but 
# in this code it would clash with the list that holds the Z-position
# of the atom being tracked.
    ZProj = int(input("Input incident atom charge Z:") or "14")

# Substrate atom mass
    MASSSUB = int(input("Input substrate atom mass:") or "28")

# Charge of target atom
    ZSUB = int(input("Input substrate atom charge Z:") or "14")

# Substrate density immediately converted to cubic angstrom
    SUBDENSITY = float(input("Input substrate density (atoms/cc):") or "5E22")/ float("1E24")

# Incident atom energy
    INELAB = float(input("Input incident atom energy (keV):") or "50")

# Substrate Window
    SUBWINDOW = int(input("Substrate Window (A):") or "2000")

# Amount of Simulations
    SIMULS = int(input("How many simulations:") or "100")

    return massp, ZProj, MASSSUB, ZSUB, SUBDENSITY, INELAB, SUBWINDOW, SIMULS