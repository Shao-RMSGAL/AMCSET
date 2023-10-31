## ===============================
## AUTHOR           : Lin Shao
## CREATE DATE      :
## PURPOSE          : Collision Simulation
## SPECIAL NOTES    : Currently under the process of translation
## ===============================
## Change History:
##  Python translation currently in progress.
##==================================


##################INPUT_STARTING_DATA###################################################

# Mass of the incident atom
massp = int(input("Input incident atom mass:") or  "14")

# Charge of incident atom
ZP = int(input("Input incident atom charge Z:") or "14")

# Substrate atom mass
MASSSUB = int(input("Input substrate atom mass:") or "28")

# Charge of target atom
ZSUB = int(input("Input substrate atom charge Z:") or "14")

# Substrate density immediately converted to cubic angstrom
SUBDENSITY = float(input("Input substrate density (atoms/cc):") or "5E22")/ float("1E24")

# Incident atom energy
INELAB = float(input("Input incident atom energy (keV):") or "50")

# Substrate Window
SUBWINDOW = int(input("Substrate Window (A)") or "2000")

# Amount of Simulations
SIMULS = int(input("How many simulations:") or "100")


#########################DEFINING_LISTS###########################################

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
