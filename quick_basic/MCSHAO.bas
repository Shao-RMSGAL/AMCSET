DECLARE SUB EMAGIC (X!)
DECLARE SUB TMAGIC (MASS1!, Z1!, MASS2!, Z2!, INELAB!, P!)
DECLARE SUB AMAGIC (THETAO!, ALPHAO!, THETA1RELATIVE!, THETA2RELATIVE!)
DECLARE FUNCTION DF! (X!, COLUMBIAVK!, Z1!, Z2!, AU!)
DECLARE FUNCTION F! (X!, COLUMBIAVK, Z1, Z2, AU)
DECLARE SUB DSCMOTT (mott(), scr(), ElectronEnergy, SubstrateZ, Theta)
DECLARE SUB IoniElecLoss (ZSUB, DENSITY, ElecE)
DECLARE SUB BremsEloss (ZSUB, DENSITY, ElecE)
DECLARE SUB IMAGE(Scan(),NewDeltaZ, NewDeltaX, NewDeltaY, NewTheta, NewAlpha, NewEP, Resolution,size)




e_range = 80000000 '7000000 'interested region depth in the unit of Angsrom
e_interval = 30 'interval number for the interest region
e_ave_range = 1 'range to averageing neighboring data
e_deltaLat = e_range / e_interval
e_deltaDep = e_range / e_interval
e_stopping = 1 '0.01 'electrons are assumed to stop at energies of 100 eV
ion_stopping = 0.04 'ions are assumed to stop at energies of 20 eV
ion_Ed = 0.04 'threshold displacement energy
Divisor = 1000 'the number of angle intervals from 0 to 180 degrees for electron bombardment
fly0 = 1000
fly1 = 1000
flyC = fly1
flyjudge = 1 '20 '100
sizeScan = 100
Resolution = 10
scandelta = 20
EMAX = 10000
switch = 1 'switch 0 is for explicit method, switch 1 is for middlepoint method, and switch 2 si for back (implicit) method




RANDOMIZE TIMER / 3

DIM RN(flyC) AS DOUBLE
DIM RL(flyC + 1)
DIM RA(flyC) AS DOUBLE
DIM RE(flyC) AS DOUBLE
DIM totalRE AS DOUBLE
DIM ElecRE AS DOUBLE
DIM EP AS DOUBLE
DIM ElecREup AS DOUBLE
DIM ElecREdown AS DOUBLE
DIM L AS DOUBLE
DIM Lselected AS DOUBLE
DIM TotalCross AS DOUBLE
DIM SHARED DiffCross AS DOUBLE
DIM SHARED IoniEloss1 AS DOUBLE
DIM SHARED IoniEloss2 AS DOUBLE
DIM SHARED BEloss AS DOUBLE
DIM E0 AS SINGLE
DIM E1 AS SINGLE
DIM E2 AS SINGLE
DIM E3 AS SINGLE
DIM E4V AS SINGLE
DIM e_output(e_interval, e_interval, 4) AS SINGLE '1 for electron location, 2 for ionization energy loss, 3 for elastic scattering, 4 for vacancy production
DIM Backup_e_output(e_interval, e_interval, 4) AS SINGLE
DIM Ltotal AS DOUBLE

Ltotal = 0.0


DIM SHARED mott(96, 5, 6) ' Creating a matrix of size 96x5x6 to save parameters for Mott scattering. 96 is the number of elements (from Z=1 to Z=96) 
ROW = 1   ‘starting row number for data reading
OPEN "C:\Users\lshao.AUTH\Downloads\Mott.txt" FOR INPUT AS #1 ‘Reading the input file for Mott scattering
DO WHILE NOT EOF(1)   ‘Continue to read until the end of the input file
    INPUT #1, mott(ROW, 1, 1), mott(ROW, 1, 2), mott(ROW, 1, 3), mott(ROW, 1, 4), mott(ROW, 1, 5), mott(ROW, 1, 6)  ‘reading first row for Z=1, numbers are for mott(Z=1, 1, 1 to 6)
    INPUT #1, mott(ROW, 2, 1), mott(ROW, 2, 2), mott(ROW, 2, 3), mott(ROW, 2, 4), mott(ROW, 2, 5), mott(ROW, 2, 6) ‘reading first row for Z=1, numbers are for mott(Z=1, 2, 1 to 6)
    INPUT #1, mott(ROW, 3, 1), mott(ROW, 3, 2), mott(ROW, 3, 3), mott(ROW, 3, 4), mott(ROW, 3, 5), mott(ROW, 3, 6) ‘reading first row for Z=1, numbers are for mott(Z=1, 3, 1 to 6)
    INPUT #1, mott(ROW, 4, 1), mott(ROW, 4, 2), mott(ROW, 4, 3), mott(ROW, 4, 4), mott(ROW, 4, 5), mott(ROW, 4, 6) ‘reading first row for Z=1, numbers are for mott(Z=1, 4, 1 to 6)
    INPUT #1, mott(ROW, 5, 1), mott(ROW, 5, 2), mott(ROW, 5, 3), mott(ROW, 5, 4), mott(ROW, 5, 5), mott(ROW, 5, 6) ‘reading first row for Z=1, numbers are for mott(Z=1, 5, 1 to 6)
    ROW = ROW + 1   ‘after finishing reading ROW=1 for Z=1,  ROW changes to 2 for Z=2

LOOP  ‘reading continues and restarts from line 75, with ROW increased by 1
CLOSE #1 ' Close the file after reading the last input line

PRINT "mott(26,1,1)="; mott(26, 1, 1) ‘print on the screen to check if the reading is accurate

DIM SHARED scr(92, 6) ‘Creating a matrix of size 96x6 to save parameters for charge screening effect. 96 is the number of elements (from Z=1 to Z=96)
ROW = 1 ‘starting row number for data reading
OPEN "C:\Users\lshao.AUTH\Downloads\screen factor.txt" FOR INPUT AS #2 'Reading the input file for charge screening effect
DO WHILE NOT EOF(2) ‘Continue to read until the end of the input file
    INPUT #2, scr(ROW, 1), scr(ROW, 2), scr(ROW, 3), scr(ROW, 4), scr(ROW, 5), scr(ROW, 6) ‘reading first row for Z=1, from scr(Z=1, 1) to scr(Z=1, to 6)

    ROW = ROW + 1  ‘after finishing reading ROW=1 for Z=1,  ROW changes to 2 for Z=2
LOOP ‘reading continues and restarts from line 97, with ROW increased by 1
CLOSE #2 ' Close the file after reading the last input line

PRINT "scr(26,1)="; scr(26, 5) ) ‘print on the screen to check if the reading is accurate

‘Lines 108 to 120 are for creating various output files to save data calculated. 
OPEN "C:\Users\lshao.AUTH\Downloads\PRANGEEd40.DAT" FOR OUTPUT AS #1
OPEN "C:\Users\lshao.AUTH\Downloads\simulEd40.DAT" FOR OUTPUT AS #2
OPEN "C:\Users\lshao.AUTH\Downloads\simulvimageEd40" FOR OUTPUT AS #3
OPEN "C:\Users\lshao.AUTH\Downloads\3DINTLimageEd40.DAT" FOR OUTPUT AS #4
OPEN "C:\Users\lshao.AUTH\Downloads\eRangeimageEd40.dat" FOR OUTPUT AS #5
OPEN "C:\Users\lshao.AUTH\Downloads\ProjectedeRangeimageEd40" FOR OUTPUT AS #6
OPEN "C:\Users\lshao.AUTH\Downloads\e-outimageEd40.dat" FOR OUTPUT AS #7
OPEN "C:\Users\lshao.AUTH\Downloads\Projected_e_rangeImageEd40.txt" FOR OUTPUT AS #9 '
OPEN "C:\Users\lshao.AUTH\Downloads\20imageEd40.txt" FOR OUTPUT AS #10
OPEN "C:\Users\lshao.AUTH\Downloads\crossimageEd40.txt" FOR OUTPUT AS #11
OPEN "C:\Users\lshao.AUTH\Downloads\imageEd40.txt" FOR OUTPUT AS #12
OPEN "C:\Users\lshao.AUTH\Downloads\VancancyEd40.txt" FOR OUTPUT AS #13
OPEN "C:\Users\lshao.AUTH\Downloads\dedx.txt" FOR OUTPUT AS #14 '''''''''''''''''


INPUT "Is this an electron bombardment (1=yes, 0=no)", iontype  ‘iontype decides whether it is an electron bombardment or ion bombardment. Default without typing = 0
IF iontype = 1 THEN ‘Lines 126 to 132 are for the case of electron bombardment
    INPUT "Input electron energy(default=10 MeV)", ElecE0 ‘if it is electron bombardment, what is the electron energy? The energy is saved as ElecE0
    IF ElecE0 = 0 THEN ElecE0 = 10000  ‘if no typing of ElecE0, it becomes default (=1MeV=10000keV)
    Elecmassp = 5.45E-4 'mass of electron in the unit of amu
    ElecCP = -1 'charge of electron
    ePRbin = 1.5E6 'depth profile interval for electron in the unit of Angstroms

ELSE      ‘if it is not electron bombardment, do nothing here, moves to the next line
END IF  ‘the end of IF command 

IF iontype = 0 THEN ‘lines 137 to 143 are for ion bombardment
    INPUT " Input incident atom mass?  (28)  ", massp 'MASSP--MASS OF INCIDENT ATOM
    IF massp = 0 THEN massp = 27.9  ‘if no typing, default is for silicon atomic mass
    INPUT " Input incident atom charge Z ?  (14)", CP 'ZP--CHARGE NUMBER OF INCIDENT ATOM
    IF CP = 0 THEN CP = 14 ‘if no typing, default is for silicon atomic number
    INPUT " Input incident atom energy (kev, default=50 keV)? ", INELAB 'IME-INCIDENT ATOM ENERGY (KEV) IN LAB COORDINATOR
    IF INELAB = 0 THEN INELAB = 50 ‘if not typing, default is 50 keV

ELSE
END IF

INPUT " Input substrate atom mass? (55.85)", MASSSUB 'MASSSUB--MASS OF SUBSTRATE ATOM
IF MASSSUB = 0 THEN MASSSUB = 55.85 ‘if no typing, default is for pure Fe substrate of atomic mass 55.85 amu
INPUT " Input substrate atom charge Z? (26)", ZSUB 'ZSUB--ATOMIC NUMBER OF TARGET ATOM
IF ZSUB = 0 THEN ZSUB = 26 ‘if not typing, default is for pure Fe of atomic number 26
INPUT "Input substrate density (atoms/cc) (for C d=8.482E22/cc)", SUBDENSITY ‘input for substrate density in the unit of atoms per cc
IF SUBDENSITY = 0 THEN SUBDENSITY = 8.482E+22 ‘if no input, atomic density is 8.482E22 per cc, for Fe
SUBDENSITY = SUBDENSITY / 1E+24 ‘the atomic density changes to the unit of atoms per angstrom^3

INPUT "SUBSTRATE WINDOW (A)", SUBWINDOW   ‘input for the depth region (in Angstrom) for plotting on the screen
IF SUBWINDOW = 0 THEN SUBWINDOW = e_range ‘if no input, the plotting window width is the full range of interest. The value of e_range is defined already. 
WINDOWX = 630   ‘Relevant to screen pixel for width, corresponding to depth in the longitudinal direction
WINDOWY = 330   ‘relevant to screen pixel for height, corresponding to range in the lateral direction


INPUT "HOW MANY SIMULATION YOU WANT(50)?", SIMULS  ‘input for how many particles to be simulated
IF SIMULS = 0 THEN SIMULS = 1000 ‘if there is no input, 1000 particles are to be simulated


DIM XP(2000) AS DOUBLE 'saving X POSTION of ion/electron projectiles and target atoms, up to 2000 in a single bombardment
DIM YP(2000) AS DOUBLE 'saving Y POSITION of ion/electron projectile or target atoms
DIM ZP(2000) AS DOUBLE 'saving Z(DEPTH) POSITION  of ion/electron projectile or target atoms
DIM CP(2000) 'saving CHARGE NUMBER of ion/electron projectile or target atoms
DIM mp(2000) 'saving MASS of ion/electron projectile or target atoms
DIM THETAP(2000) 'saving ANGLE with respect to Z AXIS
DIM ALPHAP(2000) 'saving ANGLE WITH respect to X AXIS on the YZ PLANE)
DIM EP(2000) AS DOUBLE 'saving ENERGY of ion/electron projectile or target atoms
DIM Ptype(2000) 'saving the particle type, type 0 =ion, type 1=electron
DIM ePR(2000) ‘saving profile of ion/electron projectiles along z axis (one dimension) 
DIM THETA2(2000)
DIM ALPHA2(2000)
DIM SHARED Scan(sizeScan + 1, sizeScan + 1, 8) 'AS DOUBLE 'For 2D plot of backscattered electrons

‘lines 187 to 199 for setting initial values of matrix as zero
FOR rum = 1 TO 2000
    XP(rum) = 0
    ZP(rum) = 0
    YP(rum) = 0
    CP(rum) = 0
    mp(rum) = 0
    ALPHAP(rum) = 0
    THETAP(rum) = 0
    EP(rum) = 0
    Ptype(rum) = 0
    ePR(rum) = 0
NEXT rum



SCREEN 9 ‘setting screen color
VIEW (1, 1)-(WINDOWX, WINDOWY), 15, 15 'setting screen size and color for plotting

PI = 3.1415926#   ‘setting value for PI

DOSE = 0
SPUTTER = 0    ‘setting initial value of atoms being sputtered as zero


FOR i = 0 TO sizeScan
    FOR j = 0 TO sizeScan
        Scan(i, j, 8) = 0 'zero the initial value for backscattered electron count
    NEXT j
NEXT i


FOR SIMUL = 1 TO SIMULS    ‘starting the bombardment 

    IF SIMUL - 200 * INT(SIMUL / 200) = 0 THEN RANDOMIZE TIMER / 3  ‘change random number seed every 200 particles


    IF iontype = 1 THEN 'iontype=1 means electron irradiation, lines 225 to 236 are for initial projectile conditions
        XP(1) = 0 ‘projectile starting location is at the original point with X=0
        YP(1) = 0 ‘projectile starting location is at the original point with Y=0
        ZP(1) = 0 'projectile starting location is at the original point with Z=0
        CP(1) = ElecCP   ‘electron charge
        mp(1) = Elecmassp 'electron mass
        THETAP(1) = 0 ‘initial incident direction is along the z axis
        ALPHAP(1) = 0 ‘initial projected direction vector is zero degree from x-axis
        EP(1) = ElecE0 ‘initial electron energy in keV
        Ptype(1) = 1  ‘particle type is electron
        correctFactor = 1  

    ELSE
    END IF

    IF iontype = 0 THEN 'iontype=0 means ion irradiation
        XP(1) = 0 ‘projectile starting location is at the original point with X=0
        YP(1) = 0 ‘projectile starting location is at the original point with Y=0
        ZP(1) = 0  ‘projectile starting location is at the original point with Z=0
        CP(1) = CP ‘projectile charge, i.e. CP=14 for silicon 
        mp(1) = massp ‘projectile mass, i.e. massp=28 for silicon
        THETAP(1) = 0  ‘initial incident direction is along the z axis
        ALPHAP(1) = 0  ‘initial projected direction vector is zero degree from x-axis
        EP(1) = INELAB 'initial ion bombardment energy in keV
        Ptype(1) = 0 ‘particle type is ion
    ELSE
    END IF
    'scanning ===============================

    FOR i = 0 TO sizeScan 'scanning
        FOR j = 0 TO sizeScan 'scanning
            '   Scan(i, j, 1) = e_deltaLat * (i - sizeScan / 2) 'scanning
            '  Scan(i, j, 2) = e_deltaLat * (i - sizeScan / 2) 'scanning
            Scan(i, j, 1) = scandelta * (i - sizeScan / 2) 'scanning
            Scan(i, j, 2) = scandelta * (j - sizeScan / 2) 'scanning

            Scan(i, j, 3) = 0 'depth scanning
            Scan(i, j, 4) = 0 'theta scanning
            Scan(i, j, 5) = 0 'alpha scanning
            Scan(i, j, 6) = ElecE0 'E scanning
            Scan(i, j, 7) = 0

        NEXT j 'scanning
    NEXT i 'scanning
    'scanning ======================


    rum = 1        ‘the first simulation of a new bombardment event 
    MASS2 = MASSSUB  ‘substrate atomic mass
    Z2 = ZSUB  ‘substrate atomic number (charge)
    DENSITY = SUBDENSITY ‘substrate density

    100
    IF rum = 0 THEN GOTO 600 ' If all collisions caused by one bombardment are finished,  start a new bombardment by jumping to 600

    IF ZP(rum) >= SUBWINDOW THEN 
        rum = rum – 1 ‘if electron/ion fly out of the window, stop the simulation and point to the next collision saved but not finished. 
        GOTO 100 ‘jumping to the saved collisions not yet finished. If one ion bombardment creates 900 displacements and if rum= #900 finished the collision, the pointer moves to #899 to finish the rest of collisions. #0 after rum=rum-1 means all collisions have finished. 
    ELSE
    END IF
    IF ZP(rum) < 0 THEN   ‘for the case that the particle is backscattered
        IF Ptype(rum) = 0 THEN SPUTTER = SPUTTER + 1  ‘if backscattered particle is ion, sputtering number increases by 1
        rum = rum – 1 ‘with backscattering or sputtering, no need to continue the collision, pointer moves to the next saved collision not yet finished. 
        GOTO 100


    ELSE
    END IF


    IF Ptype(rum) = 1 AND EP(rum) < flyjudge * e_stopping THEN flyC = fly0   ‘if the projectile is an electron and if the electron energy is smaller than a value, the number of flying distances to be combined for evaluation is set to be fly0.  The critical value is flyjudge X e_stopping.  Flyjudge is a number defined in the input. e_stopping is the value simulations stop when electron energy is reduced to this value.  This command is used to avoid the situation that final electron energy become negative when the flying distance combination overestimates the energy loss. 
    IF Ptype(rum) = 1 AND EP(rum) >= flyjudge * e_stopping THEN flyC = fly1 ‘if the projectile is an electron and if electron energy is above flyjudge X e_stopping, the number of flying distances to be combined is fly1. This is to allow flying distance combination if electron energy is not too low.
    IF Ptype(rum) = 1 AND EP(rum) > EMAX THEN flyC = fly0  “if the energy of the electron is above EMAX, the number of combined free flying distance is fly0. This is used to avoid the situation that flying distance combination become too long as very high energy since a single free flying distance is large at high energy. 


    ''    IF mp(rum) = MASS2 THEN PRINT #4, USING "##.###^^^^^ ##.###^^^^^ ##.###^^^^^"; XP(rum); YP(rum); ZP(rum) 'modified 01/12/01
    ''    I(INT(ZP(rum))) = I(INT(ZP(rum))) + 1
    ''   IF mp(rum) = massp THEN PR(INT(ZP(rum))) = PR(INT(ZP(rum))) + 1
    '       IF Ptype(rum) = 1 THEN PRINT #5, USING "##.###^^^^^ ##.###^^^^^ ##.###^^^^^"; XP(rum); YP(rum); ZP(rum)
    '       IF Ptype(rum) = 1 THEN ePR(1 + INT(ZP(rum) / ePRbin)) = ePR(1 + INT(ZP(rum) / ePRbin)) + 1
    IF Ptype(rum) = 1 AND EP(rum) < e_stopping THEN    ‘saving information if electrons stop. e_stopping is the threshold energy

        E0 = (XP(rum) ^ 2 + YP(rum) ^ 2) ^ 0.5 ‘radial distance from z axis
        E1 = ZP(rum)   ‘longitudinal depth of electrons when stopping
        IF E0 < SUBWINDOW AND E1 > 0 AND E1 < SUBWINDOW THEN ‘if electrons are in the valid region, not sputtered, and not beyond the window
            e_output(INT(1 + E0 / e_deltaLat), INT(1 + E1 / e_deltaDep), 1) = e_output(INT(1 + E0 / e_deltaLat), INT(1 + E1 / e_deltaDep), 1) + 1  ‘adding one to the saved counts at the 2-D position matrix of lateral distance and longitudinal depth. Both distances are divided by lateral direction interval (e_deltaLat) and longitudinal direction interval (e_deltaDep). Both intervals are predefined. INT() is to obtain the integral portion of a real number. 
        ELSE
        END IF
        rum = rum – 1   ‘after the stopped electron finishes the position counting, move to the previously saved not-yet-finished ion/atom, by pointing the index number to rum-1. 
        GOTO 100 ‘go back to start the next particle (other particles not-yet-finished, but produced from the same bombardment)


    ELSE
    END IF

    IF Ptype(rum) = 0 AND EP(rum) < ion_stopping THEN    ‘check if the new particle is an atom and the energy is below the threshold value

        E0 = (XP(rum) ^ 2 + YP(rum) ^ 2) ^ 0.5   ‘if yes, calculate the radial distance from z axis
        E1 = ZP(rum)  ‘longitudinal depth of electrons at the stopping position
        IF E0 < SUBWINDOW AND E1 > 0 AND E1 < SUBWINDOW THEN ‘judge if the ion stoops between two boundaries, surface and the maxim depth
            e_output(INT(1 + E0 / e_deltaLat), INT(1 + E1 / e_deltaDep), 4) = e_output(INT(1 + E0 / e_deltaLat), INT(1 + E1 / e_deltaDep), 4) + 1  ‘if yes, saving one count to 2-D position matrix (lateral position, longitudinal position). Both positions are divided by an distance interval. e_deltaLat is for the lateral and e_deltaDep is for the longitudinal 
        ELSE
        END IF
        rum = rum – 1  ‘after the stopped atom finishes the position counting, move to the previously saved not-yet-finished collision, by pointing the index number to rum-1. 
        GOTO 100 ‘go back to start the next particle (other particles not-yet-finished, but produced from the same bombardment)


    ELSE
    END IF


    '================== BEGIN A NEW COLLISION

    IF Ptype(rum) = 0 THEN ‘check if the new collision is for an atom, instead of an electron
        KL = 1.212 * CP(rum) ^ (7 / 6) * Z2 / (CP(rum) ^ (2 / 3) + Z2 ^ (2 / 3)) ^ (3 / 2) / (mp(rum) ^ .5) ‘parameter for electronic stopping
        PL = .5  ‘another parameter for electronic stopping
        SE = KL * (EP(rum) * 1000) ^ PL 'electronic stopping at energy EP(), in the unit of eV
        L = DENSITY ^ (-1 / 3) ‘average atomic distance of the target, in the unit of Angstrom
        P = (RND / (PI * DENSITY ^ (2 / 3))) ^ .5  ‘using random number RND to select a collision parameter


        '================== COLLISON PARAMETER

        EP(rum) = EP(rum) - 1.59 * L * DENSITY * SE / 1000 ‘get an updated energy after considering electron energy loss. L is the average atomic distance selected as the step length. DENSITY is the atomic density in the unit of atoms per angstrom  ^3. SE is the electronic stopping power. 1/1000 is to convert the energy loss from eV to keV

        IF EP(rum) < ion_stopping THEN ‘check again if adding electronic stopping leads to  ion stopping. Ion_stopping is the energy criteria to stop. The stopping criteria is not necessary to be threshold displacement energy
            rum = rum – 1 ‘after the stopped atom finishes the position counting, move to the previously saved not-yet-finished collision, by pointing the index number to rum-1. 

            GOTO 100 ‘go back to start the next particle (other particles not-yet-finished, but produced from the same bombardment)
        ELSE
        END IF
        CALL TMAGIC(mp(rum), CP(rum), MASS2, Z2, EP(rum), P)  
        'TO get recoil energy RE, deviation angle from the incident direction for the projectile THETA1RELATIVE, , recoil direction of target atom with respect to the incident direction THETA2RELATIVE
        CALL AMAGIC(THETAP(rum), ALPHAP(rum), THETA1RELATIVE, THETA2RELATIVE)  ‘obtain the new directions with respect to the original xyz coordinate for the projectile and target atom (THETAP is the angle away from the z axis, ALPHAP is the angle away from the x axis on YZ plane); THETAP(RUM), ALPHA(RUM) are angles before collision; THETA1 and THETA2 are new angles relative to the direction before the collision

        '=====TEMPERARY SAVE INFORMATION.  Using matrix at index rum+1 to save everything about the target. Note the information about project is save with index rum

        ZP(rum + 1) = ZP(rum) + L * COS(THETAP(rum))   ‘assign new depth to target
        XP(rum + 1) = XP(rum) + L * SIN(THETAP(rum)) * COS(ALPHAP(rum)) ‘assign new X to target 
        YP(rum + 1) = YP(rum) + L * SIN(THETAP(rum)) * SIN(ALPHAP(rum)) ‘assign new Y to target
        CP(rum + 1) = Z2 'CP(rum + 1) = CP(rum)   ‘transfer target charge information
        mp(rum + 1) = MASS2 'mp(RUM + 1) = mp(RUM)  ‘transfer target mass information
        THETAP(rum + 1) = THETA2   ‘transfer angle information of Recoiled target. Note “2” for target. THETA2 is a shared parameter obtained from subroutine AMAGIC 
        ALPHAP(rum + 1) = ALPHA2 'transfer angle information of  Recoiled target. Note “2” for target. ALPHA2 is a shared parameter obtained from subroutine AAMAGIC
        EP(rum + 1) = RE  ‘transfer recoil energy to target atom energy. RE is a shared parameter obtained from subroutine TMAGIC


        WINDOW SCREEN(0, 0)-(WINDOWX, WINDOWY)  ‘define/draw a window

        LOCATE 4, 47 ‘locate position for words below
        PRINT " SIMULATION:"; SIMUL; "("; SIMULS; ")" 'provide updates on how many finished out of the total


        IF mp(rum) = massp THEN LINE (INT(ZP(rum) / SUBWINDOW * WINDOWX), INT((0.5 + XP(rum) / SUBWINDOW) * WINDOWY))-(INT(ZP(rum + 1) / SUBWINDOW * WINDOWX), INT((0.5 + XP(rum + 1) / SUBWINDOW) * WINDOWY)), 4 'draw a short line from previous position to the current position. Note all positions are normalized by the window size and then change to correct pixel number. The drawing is for projectile ion
        IF mp(rum) = MASSSUB THEN LINE (INT(ZP(rum) / SUBWINDOW * WINDOWX), INT((0.5 + XP(rum) / SUBWINDOW) * WINDOWY))-(INT(ZP(rum + 1) / SUBWINDOW * WINDOWX), INT((0.5 + XP(rum + 1) / SUBWINDOW) * WINDOWY)), 1 'similar to the above but the drawing is for target atoms with a different color



        IF RE <= ion_Ed THEN ‘if target atom energy is too low. No need to save information for target since there is no displacement. Saved information in rum+1 for the target is “selectively” transferred back to rum for the projectile.  
            ZP(rum) = ZP(rum + 1)  ‘return new Z back to projectile as an update
            XP(rum) = XP(rum + 1) ‘return new X back to projectile as an update
            YP(rum) = YP(rum + 1)  ‘return new Y back to projectile as an update
            mp(rum) = mp(rum)  ‘keep the original projectile mass
            CP(rum) = CP(rum) ‘keep the original projectile charge
            THETAP(rum) = THETA1  ‘update on the new direction for projectile, obtained from AMAGIC
            ALPHAP(rum) = ALPHA1 ‘update on the new direction for projectile, obtained from AMAGIC
            EP(rum) = EP(rum) - RE ‘update energy for projectile, considering energy transfer to the target. RE obtained from TMAIG
            GOTO 100    ‘return for the next step with new position and new energy. Note the index is kept at rum. It means the pointer stays on the same projectile.  All previous saved information on rum+1 are not used if target atom does not become a displacement
        ELSE
        END IF

        '====IF RE>0.02  A target displacement is created. The new atom needs to be assigned with rum+1.  From lines 411 to 420, the information transfer for rum+1 already happened. So, only the projectile needs to be updated for rum

        XP(rum) = XP(rum + 1)    ‘projectile has the same X as the target
        YP(rum) = YP(rum + 1) 'YP(RUM) = YP(RUM + 1) ‘projectile has the same Y as the target
        ZP(rum) = ZP(rum + 1) 'ZP(RUM) = ZP(RUM + 1) ‘projectile has the same Z as the target
        mp(rum) = mp(rum) 'projectile keeps its original mass
        CP(rum) = CP(rum) 'projectile keeps its original charge
        THETAP(rum) = THETA1 'projectile has its updated angle, obtained from AMAGIC. “1” is for projectile. TEHTA1 is shared from AMAGIC
        ALPHAP(rum) = ALPHA1 'projectile has its updated angle, obtained from AMAGIC. “1” is for projectile. ALPHA1 is shared from AMAGIC
        EP(rum) = EP(rum) - RE 'new energy considering the recoil energy loss

        IF ZP(rum + 1) > 0 AND ZP(rum + 1) < SUBWINDOW THEN  ‘Since one displacement occurs, vacancy information needs to be saved, if the displacement position is within two boundaries: surface and backside.
            '            V(INT(ZP(rum + 1))) = V(INT(ZP(rum + 1))) + 1  ‘optional for vacancy counting 
            E0 = (XP(rum) ^ 2 + YP(rum) ^ 2) ^ 0.5   ‘distance away from z axis
            E1 = ZP(rum)  ‘depth along z axis
            e_output(INT(1 + E0 / e_deltaLat), INT(1 + E1 / e_deltaDep), 4) = e_output(INT(1 + E0 / e_deltaLat), INT(1 + E1 / e_deltaDep), 4) + 1  ‘adding one vacancy into output matrix for e_output (Lateral, Longitudinal, 4). “4” is for vacancy. e_deltaLat and e_deltaDep are distance interval. INT is for taking integral number. 

        ELSE
        END IF
        rum = rum + 1   ‘since a new energetic atom is produced and taking the index rum+1, the calculation pointer is re-appointed to rum+1
        GOTO 100   ‘restart a new calculation with pointer at rum+1. If a displacement is created. The newly displaced atom will finish the simulation. Then the pointer moves back to the projectile to continue
    ELSE
    END IF

    'Below is for electron irradiation

    IF Ptype(rum) = 1 THEN   ‘”1” means the particle is electron

        FOR trial2 = 1 TO flyC   ‘create a chain of random number for free flying distance within the group.  flyC is the number of free flying distances in the group
            '    RANDOMIZE TIMER / 3   
            RL(trial2) = RND    ‘Assigned random number to random number matrix with ID of trial2. Trial2 starts from 1, increase till flyC. 
        NEXT trial2     ‘repeat until all numbers are assigned

        FOR trial3 = 1 TO flyC – 1   ‘first step to rank random number from low to high. Here is to pick from the original order.
            FOR trial4 = trial3 + 1 TO flyC  ‘compare with the rest after the one being picked 
                IF RL(trial3) > RL(trial4) THEN  ‘judge whether there is random number behind is smaller than the picked one
                    temp = RL(trial3)  ‘assign a tempera saving for the originally picked number
                    RL(trial3) = RL(trial4)  ‘assigned the smaller random number to the originally picked one
                    RL(trial4) = temp  ‘transfer the save value of the original random to the one being swapped
                END IF
            NEXT trial4   ‘finish everyone behind trial3
        NEXT trial3  ‘finish all trial3 from 1 to fly-1

        
        TotalCross = 0      ‘preparing to calculate integrated cross section 
        FOR ANGLE1 = 1 TO (Divisor - 1) ‘dividing 180 degrees by Divisor, performing integration
            CALL DSCMOTT(mott(), scr(), EP(rum) * correctFactor, ZSUB, PI / (1 - 10) * (1 - 10 ^ (ANGLE1 / Divisor)))    ‘obtaining differential cross section at specific angle. The angle reading is not uniform. It has higher density close to 0.
            IF ANGLE1 = 1 THEN TotalCross = TotalCross + DiffCross * 2 * PI * SIN(PI / (1 - 10) * (1 - 10 ^ (1 / Divisor))) * ((PI / 2 / (1 - 10) * (1 - 10 ^ (2 / Divisor))) + (PI / 2 / (1 - 10) * (1 - 10 ^ (1 / Divisor))))  ‘integration concerning the first interval which starts from angle=0
            IF ANGLE1 > 1 AND ANGLE1 < (Divisor - 1) THEN TotalCross = TotalCross + DiffCross * 2 * PI * SIN(PI / (1 - 10) * (1 - 10 ^ (ANGLE1 / Divisor))) * ((PI / 2 / (1 - 10) * (1 - 10 ^ ((ANGLE1 + 1) / Divisor))) - (PI / 2 / (1 - 10) * (1 - 10 ^ ((ANGLE1 - 1) / Divisor))))  ‘integration for angle intervals excluding two boundary points
            IF ANGLE1 = (Divisor - 1) THEN TotalCross = TotalCross + DiffCross * 2 * PI * SIN(PI / (1 - 10) * (1 - 10 ^ ((Divisor - 1) / Divisor))) * (PI - (PI / 2 / (1 - 10) * (1 - 10 ^ ((Divisor - 1) / Divisor))) - (PI / 2 / (1 - 10) * (1 - 10 ^ ((Divisor - 2) / Divisor))))  ‘integration concerning the last angle point concerning the boundary at 180 degree
            '    P1 = DiffCross
            '   P2 = TotalCross
            '  PRINT #11, USING "##.###^^^^^ ##.###^^^^^ ##.###^^^^^"; 180 / PI * PI / (1 - 10) * (1 - 10 ^ (ANGLE1 / Divisor)); P1; P2 'optional for saving differential cross section and integrated cross section as a function of angle, for checking 
        NEXT ANGLE1  ‘finish from 0 to 180



        Lselected = 0    ‘prepare for flying distance assignment
        FOR trial1 = 1 TO flyC ‘starting flying distance assignment with the group. flyC is the group size

            ForLselect = RND    ‘assign random number to ForLselect
            IF ForLselect = 0 THEN ForLselect = 1E-10  ‘if random number is zero, change to a very small number
            RN(trial1) = -LOG(ForLselect) / (SUBDENSITY * 1E24 * TotalCross) + (SUBDENSITY * 1E24) ^ (-1 / 3)
            '   Random number is converted to free flying distance using total cross section obtained from integration (lines 521 to 541)
            Lselected = Lselected + RN(trial1)  ‘adding each free flying distance to get the total length of the whole group
        NEXT trial1  ‘repeat and go through all free flying distance in the group


        ForThetaSelect0 = 0  ‘prepare to identify scattering angle. ThetaSelect0 is the integrated cross section. It is zero before the integration starts 

        trial5 = 1     ‘pointer needs to be updated

        RL(flyC + 1) = 12345.   ‘assignment of a “ridiculous angle” to the last scattering angle in the group
        FOR ANGLE2 = 1 TO (Divisor - 1)   ‘go through each angle point from 0 to 180

            CALL DSCMOTT(mott(), scr(), EP(rum) * correctFactor, ZSUB, PI / (1 - 10) * (1 - 10 ^ (ANGLE2 / Divisor)))  ‘obtain Mott differential cross section as a specific angle

            IF ANGLE2 = 1 THEN ForThetaSelect1 = ForThetaSelect0 + DiffCross / TotalCross * 2 * PI * SIN(PI / (1 - 10) * (1 - 10 ^ (1 / Divisor))) * ((PI / 2 / (1 - 10) * (1 - 10 ^ (2 / Divisor))) + (PI / 2 / (1 - 10) * (1 - 10 ^ (1 / Divisor))))    ‘Integration concerning first point at angle=0 needs to be specially treated. 
            IF ANGLE2 > 1 AND ANGLE2 < (Divisor - 1) THEN ForThetaSelect1 = ForThetaSelect0 + DiffCross / TotalCross * 2 * PI * SIN(PI / (1 - 10) * (1 - 10 ^ (ANGLE2 / Divisor))) * ((PI / 2 / (1 - 10) * (1 - 10 ^ ((ANGLE2 + 1) / Divisor))) - (PI / 2 / (1 - 10) * (1 - 10 ^ ((ANGLE2 - 1) / Divisor)))) ‘integration for the middle points without two angle boundaries
            IF ANGLE2 = (Divisor - 1) THEN ForThetaSelect1 = ForThetaSelect0 + DiffCross / TotalCross * 2 * PI * SIN(PI / (1 - 10) * (1 - 10 ^ ((Divisor - 1) / Divisor))) * (PI - (PI / 2 / (1 - 10) * (1 - 10 ^ ((Divisor - 1) / Divisor))) - (PI / 2 / (1 - 10) * (1 - 10 ^ ((Divisor - 2) / Divisor))))  ‘integration concerning the last boundary needs to be specially treated

            150
            IF ForThetaSelect1 > RL(trial5) THEN  ‘if the integration cross section is larger than the first random number, do the following
                IF ANGLE2 > 1 AND ANGLE2 < Divisor - 1 THEN RA(trial5) = PI / 2 / (1 - 10) * (1 - 10 ^ ((ANGLE2 + 1) / Divisor)) + PI / 2 / (1 - 10) * (1 - 10 ^ (ANGLE2 / Divisor)) - (ForThetaSelect1 - RL(trial5)) * (PI / 2 / (1 - 10) * (1 - 10 ^ ((ANGLE2 + 1) / Divisor)) - PI / 2 / (1 - 10) * (1 - 10 ^ ((ANGLE2 - 1) / Divisor))) / (ForThetaSelect1 - ForThetaSelect0)   ‘if it happens for the middle point,  corresponding angle is read from the proportionality, judged by the distance from the right side boundary point.  RA is the scattering angle selected.
                IF ANGLE2 = 1 THEN RA(trial5) = PI / 2 / (1 - 10) * (1 - 10 ^ (2 / Divisor)) + PI / 2 / (1 - 10) * (1 - 10 ^ (1 / Divisor)) - (ForThetaSelect1 - RL(trial5)) * (PI / 2 / (1 - 10) * (1 - 10 ^ (2 / Divisor)) + PI / 2 / (1 - 10) * (1 - 10 ^ (1 / Divisor))) / (ForThetaSelect1 - ForThetaSelect0)  ‘if it happens to the first angle interval, special treatment needs since interval width differs from middle points. RA is the scattering angle selected
                IF ANGLE2 = Divisor - 1 THEN RA(trial5) = PI - (ForThetaSelect1 - RL(trial5)) * (PI - PI / 2 / (1 - 10) * (1 - 10 ^ ((Divisor - 1) / Divisor)) - PI / 2 / (1 - 10) * (1 - 10 ^ ((Divisor - 2) / Divisor))) / (ForThetaSelect1 - ForThetaSelect0) ‘if it happens to the last angle interval, special treatment needs since the interval width differs from middle points. RA is the scattering angle selected.
                trial5 = trial5 + 1    ‘moves to the next random point. The pointer is increased by 1 

                GOTO 150     ‘go back and repeat. The last pointer value (trial5) becomes flyC+1. Random number assigned for flyC+1 in the matrix RL(flyC+1) is ‘ridiculously large’ number, to make sure “GOTO 150 will not happen after flyC+1. 
            ELSE
            END IF

            ForThetaSelect0 = ForThetaSelect1   ‘for integration. The integrated value is saved as the base line, to be added with the increased value from the next angle interval, to be calculated in the next angle point. 

        NEXT ANGLE2

‘lines 615 to 622 is to go through the flying distance group, from the last one to the first one, randomly pick one before the current one and switch the value. This is used to disorder the ordered scattering angle and randomly assign it them to different free flying distances. 
        FOR trial99 = flyC TO 2 STEP -1   ‘go through all distances from the last one to the first one
            jj = INT(RND * trial99) + 1    ‘randomly pick a number smaller than trial99
            qqq = RA(trial99)  ‘save scattering angle of trial99 temporarily 
            RA(trial99) = RA(jj)   ‘transfer randomly picked scattering angle to trial99
            RA(jj) = qqq  ‘transfer temporarily save value to the randomly picked free flying distance. Hence swapping finishes
        NEXT trial99  ‘go through all random number in the group


        totalRE = 0   ‘prepare for the energy loss calculation 

        FOR trial6 = 1 TO flyC  ‘go through each free flying distance in the group

            ElecREup = ((EP(rum) + 511) * (SIN(RA(trial6))) ^ 2 + MASSSUB * 931 * 1000 * (1 - COS(RA(trial6)))) * EP(rum) * (EP(rum) + 2 * 511)   ‘for calculation of energy transfer 
            ElecREdown = (EP(rum) + MASSSUB * 931.5 * 1000) ^ 2 - EP(rum) * (EP(rum) + 2 * 511) * (COS(RA(trial6))) ^ 2 ‘for calculation of energy transfer
            RE(trial6) = ElecREup / ElecREdown ‘for calculation of energy transfer
            totalRE = totalRE + RE(trial6)  ‘adding all energy loss with the group
        NEXT trial6  ‘go through the whole free flying distance group 

        ''''''''''''''''''''''''''''''''''''''''''''''
        Ltotal = Ltotal + Lselected    ‘adding all free flying distances to get total flying distance for the whole group 
        CALL IoniElecLoss(ZSUB, SUBDENSITY * 1E24, EP(rum) * correctFactor)  ‘calculate energy loss due to ionization 

        IoniEloss2 = (IoniEloss2 + ABS(IoniEloss2)) / 2.0  ‘make sure the value is positive 


        CALL BremsELoss(ZSUB, SUBDENSITY * 1E24, EP(rum) * correctFactor) ‘calculate energy loss due to braking irradiation 
        BEloss = (BEloss + ABS(BEloss)) / 2   ‘make sure the value is positive

        ''   w10 = IoniEloss2 + BEloss 'total non-Mott scattering energy loss
        ''  w11 = EP(rum) 'electron energy 
        ''  PRINT #14, USING "##.###^^^^^ ##.###^^^^^"; w11; w10 'optional to get non-Mott energy loss as a function of energy

        DIM Energy1 AS DOUBLE
        DIM Energy2 AS DOUBLE

        Energy1 = EP(rum) - Lselected * (IoniEloss2 + BEloss) – totalRE   ‘energy of electron after both non-Mott and Mott scattering energy loss




        IF switch = 0 THEN    ‘No energy correction is needed when switch=0
        ELSE
        END IF



        IF switch = 1 THEN    ‘turn on energy correction is switch=1, correction below follows midpoint approximation
            IF correctFactor = 1 THEN    ‘”1” means the correction was not performed yet since “1” is the preassigned value, do the following
                correctFactor = (Energy1 + EP(rum)) / 2 / EP(rum)  ‘the ratio of middle energy to the starting energy for the flying distance group
                GOTO 100  ‘repeat the calculation used modified energy, utilizing correctFactor.  This means the first flying distance group will be re-peated with energy correction.  Once it is repeated,  correctFactor is not one anymore, and repeating will not happen. 
            ELSE
            END IF


            correctFactor = (Energy1 + EP(rum)) / 2 / EP(rum)   ‘calculate the correction factor of the current group, and use it for the next group
        ELSE
        END IF


        


        IF switch = 2 THEN ‘turn on the correction but the energy correction follows the implicit method
            IF correctFactor = 1 THEN ‘”1” means the correction was not performed yet since “1” is the preassigned value, do the following
                correctFactor = Energy1 / EP(rum) ‘ratio is the final energy to the initial energy of the group. 
                GOTO 100 ‘repeating the calculation for the first group. 
            ELSE
            END IF
            correctFactor = Energy1 / EP(rum) ‘calculate the correction factor of the current group, and use it for the next group

        ELSE
        END IF

        EP(rum) = Energy1   ‘assign the final energy of the group as an updated energy as the starting energy of the next group 



        '  PRINT #10, USING "##.###^^^^^"; correctFactor ‘optional for saving correctFactor


        E0 = (XP(rum) ^ 2 + YP(rum) ^ 2) ^ 0.5    ‘lateral distance from z axis
        E1 = ZP(rum)          ‘depth                           
        E2 = Lselected * IoniEloss2  ‘ionization energy loss rate
        E3 = totalRE ‘Mott scattering energy loss


        IF E0 < SUBWINDOW AND E1 > 0 AND E1 < SUBWINDOW THEN  ‘for point with the valid region between two boundaries
            e_output(INT(1 + E0 / e_deltaLat), INT(1 + E1 / e_deltaDep), 2) = e_output(INT(1 + E0 / e_deltaLat), INT(1 + E1 / e_deltaDep), 2) + E2     ‘saving ionization energy to 2-D position matrix
            e_output(INT(1 + E0 / e_deltaLat), INT(1 + E1 / e_deltaDep), 3) = e_output(INT(1 + E0 / e_deltaLat), INT(1 + E1 / e_deltaDep), 3) + E3 ‘saving Mott scattering energy loss to 2-D position matrix
        ELSE
        END IF

        '     PRINT #7, USING "##.###^^^^^ ##.###^^^^^ ##.###^^^^^ ##.###^^^^^"; E0; E1; E2; E3 ‘optional


        FOR trial9 = 1 TO flyC  ‘for each collision after each flying distance, calculate the direction 
            ZP(rum + trial9) = ZP(rum + trial9 - 1) + RN(trial9) * 10 ^ 8 * COS(THETAP(rum + trial9 - 1))
            XP(rum + trial9) = XP(rum + trial9 - 1) + RN(trial9) * 10 ^ 8 * SIN(THETAP(rum + trial9 - 1)) * COS(ALPHAP(rum + trial9 - 1))
            YP(rum + trial9) = YP(rum + trial9 - 1) + RN(trial9) * 10 ^ 8 * SIN(THETAP(rum + trial9 - 1)) * SIN(ALPHAP(rum + trial9 - 1))



‘calculate scattering angles after each collision within the free flying distance group
            THETA1RELATIVE = RA(trial9)  ‘scattering angle from each Mott scattering with respect to electron flying direction prior to collision 
            THETA2RELATIVE = (PI - RA(trial9)) / 2 ‘approximation for target atom
            CALL AMAGIC(THETAP(rum + trial9 - 1), ALPHAP(rum + trial9 - 1), THETA1RELATIVE, THETA2RELATIVE)  ‘convert to angles with respect to xyz coordinate

‘assign scattering angles to projectile and target atoms for each collision
            THETAP(rum + trial9) = THETA1   ‘for electron
            ALPHAP(rum + trial9) = ALPHA1  ‘for electron
            THETA2(rum + trial9) = THETA2 ‘for target atom
            ALPHA2(rum + trial9) = ALPHA2 ‘for target atom

        NEXT trial9    ‘finish all distance within the group


‘print on the screen, number of electrons simulated out of the total to be calculated 
        LOCATE 4, 47 'modified 11/21/23
        PRINT " SIMULATION:"; SIMUL; "("; SIMULS; ")" 'modified 01/12/01
        


        IF Ptype(rum) = 1 THEN LINE (INT(ZP(rum) / SUBWINDOW * WINDOWX), INT((0.5 + XP(rum) / SUBWINDOW) * WINDOWY))-(INT(ZP(rum + flyC) / SUBWINDOW * WINDOWX), INT((0.5 + XP(rum + flyC) / SUBWINDOW) * WINDOWY)), 0 'plot a line connecting the current and next electron position
        '    IF Ptype(rum) = 0 THEN LINE (INT(ZP(rum) / SUBWINDOW * WINDOWX), INT((0.5 + XP(rum) / SUBWINDOW) * WINDOWY))-(INT(ZP(rum + flyC) / SUBWINDOW * WINDOWX), INT((0.5 + XP(rum + flyC) / SUBWINDOW) * WINDOWY)), 1 'Option to plot for target atom


        NewDeltaZ = ZP(rum + flyC) - ZP(rum)
        NewDeltaX = XP(rum + flyC) - XP(rum)
        NewDeltaY = YP(rum + flyC) - YP(rum)
        NewTheta = THETAP(rum + flyC)
        NewAlpha = ALPHAP(rum + flyC)
        NewEP = EP(rum)

        '   IF NewTheta > PI / 2 AND ZP(rum + flyC) < 0 THEN PRINT " Theta="; NewTheta; "ZP(rum)="; ZP(rum + flyC)

        '    IF NewTheta > PI / 2 AND ZP(rum + flyC) < 0 THEN
        '   PRINT "NewTheta="; NewTheta; "ZP(rum+flyC)="; ZP(rum + flyC)
        '    INPUT xxx
        '  ELSE
        ' END IF


        '''''''''    CALL IMAGE(Scan(), NewDeltaZ, NewDeltaX, NewDeltaY, NewTheta, NewAlpha, NewEP, Resolution, sizeScan)

        '    INPUT bnb
        '      SUB IMAGE (Scan(),DeltaZima, DeltaXima, DeltaYima, Thetaima, Alphaima, Eima, Resolution, size)


‘lines 797 to 805, update position after each flying distance group, transfer other information needed.  
        ZP(rum) = ZP(rum + flyC)
        XP(rum) = XP(rum + flyC)
        YP(rum) = YP(rum + flyC)
        THETAP(rum) = THETAP(rum + flyC)
        ALPHAP(rum) = ALPHAP(rum + flyC)
        CP(rum) = CP(rum) 'CP(rum + 1) = CP(rum)
        mp(rum) = mp(rum) 'mp(RUM + 1) = mp(RUM)
        EP(rum) = EP(rum) '- ElecRE - Lselected * (SUBDENSITY * 1E24) * TotalelasticE
        Ptype(rum) = Ptype(rum)



        mmmm = 1
        FOR trial8 = 1 TO flyC

            IF RE(trial8) >= ion_Ed THEN
                ZP(rum + mmmm) = ZP(rum + trial8)
                XP(rum + mmmm) = XP(rum + trial8)
                YP(rum + mmmm) = YP(rum + trial8)
                CP(rum + mmmm) = Z2
                mp(rum + mmmm) = MASS2
                THETAP(rum + mmmm) = THETA2(rum + trial8)
                ALPHAP(rum + mmmm) = ALPHA2(rum + trial8)
                Ptype(rum + mmmm) = 0
                EP(rum + mmmm) = RE(trial8)



                E0 = (XP(rum + mmmm) ^ 2 + YP(rum + mmmm) ^ 2) ^ 0.5
                E1 = ZP(rum + mmmm)
                IF E0 < SUBWINDOW AND E1 > 0 AND E1 < SUBWINDOW THEN
                    e_output(INT(1 + E0 / e_deltaLat), INT(1 + E1 / e_deltaDep), 4) = e_output(INT(1 + E0 / e_deltaLat), INT(1 + E1 / e_deltaDep), 4) + 1
                    'the above is to record one vacancy created by electron
                    E4V = EP(rum)
                    PRINT #13, USING "##.###^^^^^ ##.###^^^^^"; E4V; 1
                ELSE
                END IF
                '       INPUT XXX
                mmmm = mmmm + 1
            ELSE
            END IF

        NEXT trial8







        rum = rum + (mmmm - 1)



        '  B1 = Ltotal
        '   B2 = EP(rum)
        '     B3 = THETAP(rum)
        '     B4 = Lselected

        '     PRINT #10, USING "##.###########^^^^^ ##.###^^^^^ ##.###^^^^^"; B2; B4; B4 / B2


        '      IF THETAP(rum) = 0 THEN PRINT "something is wrong here and check thetap, why =0"
        GOTO 100

    ELSE
    END IF
    '  ELSE
    ' END IF



    'end of electron irradiation
    600
NEXT SIMUL


'''''''FOR m = 0 TO sizeScan
''''''FOR mm = 0 TO sizeScan
'''''''PRINT #12, USING "##.###^^^^^ ##.###^^^^^ ##.###^^^^^"; scandelta * (m - sizeScan / 2); scandelta * (mm - sizeScan / 2); Scan(m, mm, 8)
'''''NEXT mm
''''''''NEXT m








'''''FOR RRR = 0 TO SUBWINDOW
'''''PRINT #1, USING "##.###^^^^^ ##.###^^^^^"; RRR; PR(RRR)
'''''NEXT RRR
''''FOR RRRR = 0 TO SUBWINDOW
'''''PRINT #2, USING "##.###^^^^^ ##.###^^^^^"; RRRR; I(RRRR)
'''''PRINT #3, USING "##.###^^^^^ ##.###^^^^^"; RRRR; V(RRRR)

'''''NEXT RRRR

FOR N0 = 1 TO e_interval
    FOR N1 = 1 TO e_interval
        FOR N2 = 1 TO 4
            Backup_e_output(N0, N1, N2) = e_output(N0, N1, N2)
        NEXT N2
    NEXT N1
NEXT N0


FOR R0 = 1 + e_ave_range TO e_interval - e_ave_range
    FOR R1 = 1 + e_ave_range TO e_interval - e_ave_range
        combine1 = 0
        combine2 = 0
        combine3 = 0
        combine4 = 0

        FOR R2 = -e_ave_range TO e_ave_range
            FOR R22 = -e_ave_range TO e_ave_range

                combine1 = Backup_e_output(R0 + R2, R1 + R22, 1) + combine1
                combine2 = Backup_e_output(R0 + R2, R1 + R22, 2) + combine2
                combine3 = Backup_e_output(R0 + R2, R1 + R22, 3) + combine3
                combine4 = Backup_e_output(R0 + R2, R1 + R22, 4) + combine4
            NEXT R22
        NEXT R2
        e_output(R0, R1, 1) = combine1 / (2 * e_ave_range + 1) ^ 2
        e_output(R0, R1, 2) = combine2 / (2 * e_ave_range + 1) ^ 2
        e_output(R0, R1, 3) = combine3 / (2 * e_ave_range + 1) ^ 2
        e_output(R0, R1, 4) = combine4 / (2 * e_ave_range + 1) ^ 2
    NEXT R1
NEXT R0


FOR W0 = 1 TO e_interval
    accum = 0
    FOR w1 = 1 TO e_interval
        accum = accum + e_output(w1, W0, 1)
    NEXT w1
    PRINT #9, USING "##.###^^^^^ ##.###^^^^^"; (W0 - 1 + 0.5) * e_deltaDep / 10000; accum / SIMULS / ((e_deltaDep / 1E4) * 1E4 * 1E4)
NEXT W0


test_e = 0
test_IoniE = 0
test_RecE = 0
test_Disp = 0

FOR R3 = 1 TO e_interval
    FOR R4 = 1 TO e_interval
        unit_volume = 2 * PI * (R3 + 0.5) * e_deltaLat * e_deltaLat * e_deltaDep / 10 ^ 12 'unit is micron^3
        PRINT #7, USING "##.###^^^^^ ##.###^^^^^ ##.###^^^^^ ##.###^^^^^ ##.###^^^^^ ##.###^^^^^"; (R4 - 1 + 0.5) * e_deltaDep / 10000; (R3 - 1 + 0.5) * e_deltaLat / 10000; e_output(R3, R4, 1) / unit_volume / SIMULS; e_output(R3, R4, 2) / unit_volume / SIMULS; e_output(R3, R4, 3) / unit_volume / SIMULS; e_output(R3, R4, 4) / unit_volume / SIMULS
        test_e = test_e + unit_volume * e_output(R3, R4, 1) / unit_volume / SIMULS
        test_IoniE = test_IoniE + e_output(R3, R4, 2) / unit_volume / SIMULS * unit_volume
        test_RecE = test_RecE + e_output(R3, R4, 3) / unit_volume / SIMULS * unit_volume
        test_Disp = test_Disp + e_output(R3, R4, 4) / unit_volume / SIMULS * unit_volume
    NEXT R4
NEXT R3
PRINT "integrated total electron inside ="; test_e
PRINT "integrtated total ioniztion (inelastic) energy loss="; test_IoniE
PRINT "integrated total elastic energy loss="; test_RecE
PRINT "integrated total displacement creation="; test_Disp
END

SUB AMAGIC (THETAO, ALPHAO, THETA1RELATIVE, THETA2RELATIVE)
    SHARED PI
    SHARED THETA1 'RESULT TO MAIN PROGRAM, DEFLECTION ANGLE TO ORIGINAL C.
    SHARED ALPHA1 'RESULT TO MAIN PROGRAM
    SHARED THETA2 'RESULT TO MAIN PROGRAM
    SHARED ALPHA2 'RESULT TO MAIN PROGRAM
    ALPHA1RELATIVE = RND * 2 * PI
    ALPHA2RELATIVE = ALPHA1RELATIVE + PI
    '=====FOR PRIMARY ION
    X1 = SIN(THETA1RELATIVE) * COS(ALPHA1RELATIVE)
    Y1 = SIN(THETA1RELATIVE) * SIN(ALPHA1RELATIVE)
    Z1 = COS(THETA1RELATIVE)

    Y0 = Y1 * COS(THETAO) + Z1 * SIN(THETAO)
    Z0 = -Y1 * SIN(THETAO) + Z1 * COS(THETAO)
    X0 = X1

    Z = Z0
    X = X0 * SIN(ALPHAO) + Y0 * COS(ALPHAO)
    Y = -X0 * COS(ALPHAO) + Y0 * SIN(ALPHAO)

    IF Z > O THEN THETA1 = ATN(SQR(X ^ 2 + Y ^ 2) / Z)
    IF Z = 0 THEN THETA1 = PI / 2
    IF Z < 0 THEN THETA1 = PI + ATN(SQR(X ^ 2 + Y ^ 2) / Z)

    IF SIN(THETA1) <> 0 THEN
        IF X > 0 THEN ALPHA1 = ATN(Y / X)
        IF X = 0 THEN ALPHA1 = PI - SGN(Y) * PI / 2
        IF X < 0 THEN ALPHA1 = PI + ATN(Y / X)
    ELSE
        ALPHA1 = 0
    END IF

    '========FOR RECOIL TARGET ION

    X1 = SIN(THETA2RELATIVE) * COS(ALPHA2RELATIVE)
    Y1 = SIN(THETA2RELATIVE) * SIN(ALPHA2RELATIVE)
    Z1 = COS(THETA2RELATIVE)

    Y0 = Y1 * COS(THETAO) + Z1 * SIN(THETAO)
    Z0 = -Y1 * SIN(THETAO) + Z1 * COS(THETAO)
    X0 = X1

    Z = Z0
    X = X0 * SIN(ALPHAO) + Y0 * COS(ALPHAO)
    Y = -X0 * COS(ALPHAO) + Y0 * SIN(ALPHAO)

    IF Z > O THEN THETA2 = ATN(SQR(X ^ 2 + Y ^ 2) / Z)
    IF Z = 0 THEN THETA2 = PI / 2
    IF Z < 0 THEN THETA2 = PI + ATN(SQR(X ^ 2 + Y ^ 2) / Z)

    IF SIN(THETA2) <> 0 THEN
        IF X > 0 THEN ALPHA2 = ATN(Y / X)
        IF X = 0 THEN ALPHA2 = PI - SGN(Y) * PI / 2
        IF X < 0 THEN ALPHA2 = PI + ATN(Y / X)
    ELSE
        ALPHA2 = 0
    END IF


END SUB

FUNCTION DF (X, COLUMBIAVK, Z1, Z2, AU)
    DF = F(X, COLUMBIAVK, Z1, Z2, AU) / X + COLUMBIAVK / X * (.35 * EXP(-.3 / X / AU) * .3 / AU + .55 * EXP(-1.2 / X / AU) * 1.2 / AU + .1 * EXP(-6 / X / AU) * 6 / AU)
END FUNCTION

SUB EMAGIC (X)
    SHARED ATOMD
    SHARED ELOSS

    ' X IS THE CURVE FITTING FROM TRIM AND UNIT IS EV, X UNIT IS KEV
    IF X * 1000 > 1 AND X * 1000 <= 10 THEN ELOSS = ATOMD * (.7122 + .1026 * X * 1000) / 1000
    IF X * 1000 > 10 AND X * 1000 <= 100 THEN ELOSS = ATOMD * (1.671 + .0252 * X * 1000) / 1000
    IF X * 1000 > 100 AND X * 1000 <= 400 THEN ELOSS = ATOMD * (5.0767 + .0044 * X * 1000) / 1000
    IF X * 1000 > 1000 AND X * 1000 <= 10000 THEN ELOSS = ATOMD * (9.28 + .00144 * X * 1000) / 1000
    IF X * 1000 > 10000 THEN ELOSS = ATOMD * (9.28 + .00144 * X * 1000) / 1000
END SUB

FUNCTION F (X, COLUMBIAVK, Z1, Z2, AU)
    '===========================UNIVERAIAL SCREENING POTENTIAL=================
    IF X = O THEN
        F = 0
    ELSE
        F = COLUMBIAVK * X * (.35 * EXP(-.3 / X / AU) + .55 * EXP(-1.2 / X / AU) + .1 * EXP(-6 / X / AU))
    END IF
    '=========================COLUMBIA POTENTIAL==============================
    'F = COLUMBIAVK * X


END FUNCTION

SUB TMAGIC (MASS1, Z1, MASS2, Z2, INELAB, P)

    SHARED ATOMD 'WHAT SHARED IS ALL THE NEEDED DATA BUT NOT DEFINED HERE
    'AND ALL THE RESULT NEEDING TRANSFERED TO MAIN PROGRAM
    SHARED THETA1RELATIVE 'RESULT TO MAIN PROGRAM
    SHARED THETA2RELATIVE 'RESULT TO MAIN PROGRAM
    SHARED RE 'RESULT TO MAIN PROGRAM

    PI = 3.1415926#
    '=============================INITIAL PARAMETER============================

    COLUMBIAVK = .0143992# * Z1 * Z2 'POTENTIAL V=COLUMBIAVK(ANSTRON*KEV)*Z1*Z2/R '
    MC = MASS1 * MASS2 / (MASS1 + MASS2) 'MC--REDUCED MASS IN CENTER-MASS COORDINATOR
    INVLAB = SQR(INELAB * 2 / MASS1) 'INVLAB--INCIDENT VOLOCITY IN LAB C.
    EC = 1 / 2 * MC * INVLAB ^ 2 'EC--INITIAL ENERGY IN CM
    AU = .8854 * .529 / (Z1 ^ .5 + Z2 ^ .5) ^ (2 / 3)
    ELINHARD = EC * AU / COLUMBIAVK

    '==============FIND RMIN FOR DIFFERENT ENERGY======================
    AA = P ^ 2
    IF AA = 0 THEN AA = .00001
    BB = COLUMBIAVK / EC
    CC = -1
    COLUMRMIN = 1 / 2 / AA * (-BB + SQR(BB ^ 2 - 4 * AA * CC)) 'COLUMRMIN IS 1/RMIN, RMIN IS THE MINIMUM R UNDER COLUMBIA POTENTIAL

    CALTIME = 1
    300 RMINTRY1 = COLUMRMIN
    DV = ABS(-DF(RMINTRY1, COLUMBIAVK, Z1, Z2, AU) / EC - 2 * P ^ 2 * RMINTRY1)
    IF ABS(DV) < .000001 THEN DV = .1
    RMINTRY2 = RMINTRY1 + (1 - F(RMINTRY1, COLUMBIAVK, Z1, Z2, AU) / EC - P ^ 2 * RMINTRY1 ^ 2) / DV
    COLUMRMIN = RMINTRY2
    CALTIME = CALTIME + 1
    IF CALTIME > 10000 THEN GOTO 350
    IF ABS(RMINTRY2 - RMINTRY1) > .00001 THEN GOTO 300
    350 RMIN = (RMINTRY2 + RMINTRY1) / 2




    '=====================CALCULATE DEFELCTION ANGLE=========================



    RBIERSACK = 2 * (EC - F(RMIN, COLUMBIAVK, Z1, Z2, AU)) / RMIN ^ 2 / DF(RMIN, COLUMBIAVK, Z1, Z2, AU)

    BBIERSACK = P / AU
    ROBIERSACK = 1 / (RMIN * AU)
    RCBIERSACK = RBIERSACK / AU

    C1BIERSACK = .6743
    C2BIERSACK = .009611
    C3BIERSACK = .005175
    C4BIERSACK = 10!
    C5BIERSACK = 6.314
    ALTHABIERSACK = 1 + C1BIERSACK * ELINHARD ^ (-1 / 2)
    BELTABIERSACK = (C2BIERSACK + ELINHARD ^ (1 / 2)) / (C3BIERSACK + ELINHARD ^ (1 / 2))
    GAMABIERSACK = (C4BIERSACK + ELINHARD) / (C5BIERSACK + ELINHARD)
    ABIERSACK = 2 * ALTHABIERSACK * ELINHARD * BBIERSACK ^ BELTABIERSACK
    GBIERSACK = GAMABIERSACK * 1 / ((1 + ABIERSACK ^ 2) ^ (1 / 2) - ABIERSACK)
    DELTABIERSACK = ABIERSACK * (ROBIERSACK - BBIERSACK) / (1 + GBIERSACK)
    'PRINT "COS(H/2)="; (BBIERSACK + RCBIERSACK + DELTABIERSACK) / (ROBIERSACK + RCBIERSACK)
    'END
    IF P = 0 THEN
        CALPHA1 = PI
    ELSE
        CALPHA1 = 2 * ATN(SQR((ROBIERSACK + RCBIERSACK) ^ 2 / (BBIERSACK + RCBIERSACK + DELTABIERSACK) ^ 2 - 1))
    END IF
    ' ===================DEFLECTION ANGLE IN LAB COORDINATOR===========
    COSPLUSMASS = COS(CALPHA1) + MASS1 / MASS2
    IF COSPLUSMASS = 0 THEN THETA1RELATIVE = PI / 2
    IF COSPLUSMASS < 0 THEN THETA1RELATIVE = PI + ATN(SIN(CALPHA1) / COSPLUSMASS)
    IF COSPLUSMASS > 0 THEN THETA1RELATIVE = ATN(SIN(CALPHA1) / COSPLUSMASS)

    '===================CALCUTE T & DIRECTION OF TARGET ATOM=============
    RE = 4 * EC * MC / MASS2 * SIN(CALPHA1 / 2) * SIN(CALPHA1 / 2)
    'RE IS ENERGY TRANSIMITED TO TARGET ATOM'
    '===================RECOILED DIRECTION==================================
    THETA2RELATIVE = (PI - CALPHA1) / 2 'RECOILED DIRECTION
    'PRINT "RE="; RE, "CALPHA1="; CALPHA1
END SUB

SUB DSCMOTT (mott(), scr(), ElectronEnergy, SubstrateZ, Theta)


    SHARED DiffCross AS DOUBLE

    PI = 3.1415926#
    beta = (1 - (ElectronEnergy / 511 + 1) ^ (-2)) ^ 0.5
    betaAVE = 0.7181287
    gamma = ElectronEnergy / 511 + 1
    elecmass = 1 'specific for a.u. unit,  only for Feq calcu.
    Vc = 137 'speed of light for a.u. unit, only for Feq calcu.
    re1 = 2.817938E-13 'in the unit of cm)

    Rmott = 0
    FOR xx = 1 TO 5
        alpha1 = 0
        FOR yy = 1 TO 6
            alpha1 = alpha1 + mott(SubstrateZ, xx, yy) * (beta - betaAVE) ^ (yy - 1)
        NEXT yy
        Rmott = Rmott + alpha1 * (1 - COS(Theta)) ^ ((xx - 1) / 2)
    NEXT xx
    DSC = (SubstrateZ * re1) ^ 2 * (1 - beta ^ 2) / beta ^ 4 * (1 - COS(Theta)) ^ (-2) 'Differential Ruth cross section
    moment = 2 * beta * gamma * elecmass * Vc * SIN(Theta / 2)
    Feq = 0
    FOR ww = 1 TO 3
        Feq = Feq + scr(SubstrateZ, ww) * scr(SubstrateZ, ww + 3) ^ 2 / (scr(SubstrateZ, ww + 3) ^ 2 + moment ^ 2)
    NEXT ww
    DiffCross = Rmott * DSC * (1 - Feq) ^ 2

    'PRINT "Rmott="; Rmott
    'PRINT "DSC="; DSC
    'PRINT "Feq="; Feq
    'PRINT "moment="; moment
    'PRINT "SubstrateZ="; SubstrateZ
    'PRINT "scr(SubstrateZ,1)="; scr(SubstrateZ, 1)
    '  PRINT "DiffCross="; DiffCross

    'INPUT XXXXXX

END SUB

SUB IoniElecLoss (ZSUB, DENSITY, ElecE)
    ' SHARED IoniEloss2 AS DOUBLE
    '  SHARED IoniEloss1 AS DOUBLE

    ''  PRINT "ElecE="; ElecE

    beta = (1 - (ElecE / 511 + 1) ^ (-2)) ^ 0.5
    KforLoss = 0.734 * ZSUB ^ 0.037
    IF ZSUB < 13 THEN JforLoss = 0.0115 * ZSUB
    IF ZSUB >= 13 THEN JforLoss = 0.00976 * ZSUB + 0.0585 * ZSUB ^ (-0.19)

    JforLoss = JforLoss / (1 + KforLoss * JforLoss / ElecE)

    IoniEloss2 = 153.55 * DENSITY * ZSUB / (6.022 * 10 ^ 23) / beta ^ 2 * (LOG(511 * beta ^ 2 * ElecE / (2 * JforLoss ^ 2 * (1 - beta ^ 2))) - LOG(2) * (2 * (1 - beta ^ 2) ^ 0.5 - 1 + beta ^ 2) + 1 - beta ^ 2 + 1 / 8 * (1 - (1 - beta ^ 2) ^ 0.5) ^ 2)

    '  JforLoss = JforLoss / (1 + KforLoss * JforLoss / ElecE)
    '  IoniEloss1 = 7.85 * (10 ^ 4) * DENSITY * ZSUB / (6.022 * 10 ^ 23) / ElecE * LOG(1.166 * ElecE / JforLoss)
    '  IF IoniEloss2 < IoniEloss1 THEN IoniEloss2 = IoniEloss1

END SUB

SUB BremsELoss (ZSUB, DENSITY, ElecE)
    '   SHARED BEloss AS DOUBLE
    BEloss = 1.4 * 10 ^ (-4) * DENSITY / (6.022 * 10 ^ 23) * ZSUB * (ZSUB + 1) * (ElecE + 511) * (4 * LOG(2 * (ElecE + 511) / 511) - 4 / 3)
    'PRINT "what the help density="; DENSITY
    ' PRINT "what the help ElecE="; ElecE
    '    PRINT "BEloss inside="; BEloss
END SUB


SUB IMAGE (Scan(), DeltaZima, DeltaXima, DeltaYima, Thetaima, Alphaima, Eima, Resolution, size)
    '  Scan(sizeScan + 1, sizeScan + 1, 6)
    SHARED Scan()
    '  PRINT "ScanZ="; Scan(20, 20, 3)
    ' INPUT xx

    PI = 3.14159
    m = 1

    FOR i = 0 TO size
        FOR j = 0 TO size
            FOR k = 1 TO Resolution
                Scan(i, j, 1) = Scan(i, j, 1) + DeltaXima / Resolution
                Scan(i, j, 2) = Scan(i, j, 2) + DeltaYima / Resolution
                Scan(i, j, 3) = Scan(i, j, 3) + DeltaZima / Resolution
                999

                IF (Scan(i, j, 1) ^ 2 + Scan(i, j, 2) ^ 2 + (Scan(i, j, 3) - 200) ^ 2) ^ 0.5 < 50 THEN

                    Scan(i, j, 1) = Scan(i, j, 1) + DeltaXima / Resolution
                    Scan(i, j, 2) = Scan(i, j, 2) + DeltaYima / Resolution
                    Scan(i, j, 3) = Scan(i, j, 3) + DeltaZima / Resolution


                    GOTO 999

                ELSE
                END IF
            NEXT k
            Scan(i, j, 4) = Thetaima
            Scan(i, j, 5) = Alphaima
            Scan(i, j, 6) = Eima

            '  PRINT "Scan(Z)="; Scan(i, j, 3); "DeltaZima="; DeltaZima; "Thetaima="; Thetaima; "Resolution="; Resolution; "size="; size
            ' INPUT xx


            '  IF Scan(i, j, 3) < 0 AND Scan(i, j, 7) = 0 AND Scan(i, j, 4) > 150 * PI / 180 AND Scan(i, j, 4) < 170 * PI / 180 THEN
            ' PRINT #12, USING "##.###^^^^^ ##.###^^^^^ ##.###^^^^^ ##.###^^^^^ ##.###^^^^^"; Scan(i, j, 1); Scan(i, j, 2); Scan(i, j, 4); Scan(i, j, 5); Scan(i, j, 6)
            ' Scan(i, j, 7) = 1
            '        PRINT "something occurs here"
            'ELSE
            'END IF
            '  INPUT bb


            IF Scan(i, j, 3) < 0 AND Scan(i, j, 7) = 0 THEN
                Scan(i, j, 8) = Scan(i, j, 8) + 1
                Scan(i, j, 7) = 1
            ELSE
            END IF

            ' PRINT #12, USING "##.###^^^^^ ##.###^^^^^ ##.###^^^^^ ##.###^^^^^ ##.###^^^^^"; Scan(i, j, 1); Scan(i, j, 2); Scan(i, j, 4); Scan(i, j, 5); Scan(i, j, 6)
            ' Scan(i, j, 7) = 1
            '        PRINT "something occurs here"
            'ELSE
            'END IF



        NEXT j
    NEXT i
    '    PRINT "Scan(20,20,3)="; Scan(20, 20, 3); "Scan(20,20,4)="; Scan(20, 20, 4)
END SUB


'       Scan(i, j, 1) = e_deltaLat * (i - sizeScan / 2) 'scanning
'      Scan(i, j, 2) = e_deltaLat * (i - sizeScan / 2) 'scanning
'     Scan(i, j, 3) = 0 'depth scanning
'    Scan(i, j, 4) = 0 'theta scanning
'   Scan(i, j, 5) = 0 'alpha scanning
'  Scan(i, j, 6) = ElecE0 'E scanning

' FOR i = 0 TO sizeScan 'scanning
'    FOR j = 0 TO sizeScan 'scanning


