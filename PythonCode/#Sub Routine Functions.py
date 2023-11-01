#Sub Routine Functions
#Kenneth

import numpy as np
import random as random


############## EMagic Function ###############
def EMagic(ATOMD,X):
    ELOSS = 0

    if 1 <= X * 1000 <= 10:
        ELOSS = ATOMD * (0.7122 + 0.1026 * X * 1000) / 1000
    elif 10 < X * 1000 <= 100:
        ELOSS = ATOMD * (1.671 + 0.0252 * X * 1000) / 1000
    elif 100 < X * 1000 <= 400:
        ELOSS = ATOMD * (5.0767 + 0.0044 * X * 1000) / 1000
    elif 1000 < X * 1000 <= 10000:
        ELOSS = ATOMD * (9.28 + 0.00144 * X * 1000) / 1000
    elif X * 1000 > 10000:
        ELOSS = ATOMD * (9.28 + 0.00144 * X * 1000) / 1000
    return ELOSS


############## TMagic Function ###############
def TMagic(ATOMD,MASS1, Z1, MASS2, Z2, INELAB, P):
    PI = math.pi

    COLUMBIAVK = 0.0143992 * Z1 * Z2
    MC = MASS1 * MASS2 / (MASS1 + MASS2)
    INVLAB = math.sqrt(2 * INELAB * 1000 / MASS1)
    EC = 0.5 * MC * INVLAB ** 2
    AU = 0.8854 * 0.529 / (Z1 ** 0.5 + Z2 ** 0.5) ** (2 / 3)
    ELINHARD = EC * AU / COLUMBIAVK

    AA = P ** 2
    if AA == 0:
        AA = 0.00001
    BB = COLUMBIAVK / EC
    CC = -1
    COLUMRMIN = 1 / (2 * AA) * (-BB + math.sqrt(BB ** 2 - 4 * AA * CC))

    CALTIME = 1
    RMINTRY1 = COLUMRMIN
    DV = abs(-DF(RMINTRY1, COLUMBIAVK, Z1, Z2, AU) / EC - 2 * P ** 2 * RMINTRY1)
    if abs(DV) < 0

############## AMagic Function ###############
def AMagic(THETAO, ALPHAO, THETA1RELATIVE, THETA2RELATIVE):
    PI = math.pi
    ALPHA1RELATIVE = random.random() * 2 * PI
    ALPHA2RELATIVE = PI + ALPHA1RELATIVE

    X1 = math.sin(THETA1RELATIVE) * math.cos(ALPHA1RELATIVE)
    Y1 = math.sin(THETA1RELATIVE) * math.sin(ALPHA1RELATIVE)
    Z1 = math.cos(THETA1RELATIVE)

    Y0 = Y1 * math.cos(THETAO) + Z1 * math.sin(THETAO)
    Z0 = -Y1 * math.sin(THETAO) + Z1 * math.cos(THETAO)
    X0 = X1

    Z = Z0
    X = X0 * math.sin(ALPHAO) + Y0 * math.cos(ALPHAO)
    Y = -X0 * math.cos(ALPHAO) + Y0 * math.sin(ALPHAO)

    if Z > 0:
        THETA1 = math.atan(math.sqrt(X**2 + Y**2) / Z)
    elif Z == 0:
        THETA1 = PI / 2
    else:
        THETA1 = PI + math.atan(math.sqrt(X**2 + Y**2) / Z)

    if math.sin(THETA1) != 0:
        if X > 0:
            ALPHA1 = math.atan(Y / X)
        elif X == 0:
            ALPHA1 = PI - math.copysign(PI / 2, Y)
        else:
            ALPHA1 = PI + math.atan(Y / X)
    else:
        ALPHA1 = 0

    X1 = math.sin(THETA2RELATIVE) * math.cos(ALPHA2RELATIVE)
    X2 = math.sin(THETA2RELATIVE) * math.sin(ALPHA2RELATIVE)
    Z1 = math.cos(THETA2RELATIVE)

    Y0 = Y1 * math.cos(THETAO) + Z1 * math.sin(THETAO)
    Z0 = -Y1 * math.sin(THETAO) + Z1 * math.cos(THETAO)
    X0 = X1

    Z = Z0
    X = X0 * math.sin(ALPHAO) + Y0 * math.cos(ALPHAO)
    Y = -X0 * math.cos(ALPHAO) + Y0 * math.sin(ALPHAO)

    if Z > 0:
        THETA2 = math.atan(math.sqrt(X**2 + Y**2) / Z)
    elif Z == 0:
        THETA2 = PI / 2
    else:
        THETA2 = PI + math.atan(math.sqrt(X**2 + Y**2) / Z)

    if math.sin(THETA2) != 0:
        if X > 0:
            ALPHA2 = math.atan(Y / X)
        elif X == 0:
            ALPHA2 = PI - math.copysign(PI / 2, Y)
        else:
            ALPHA2 = PI + math.atan(Y / X)
    else:
        ALPHA2 = 0

    THETA1RELATIVE = THETA1
    ALPHAO = ALPHA1
    THETA2RELATIVE = THETA2

    return THETA1RELATIVE, ALPHAO, THETA2RELATIVE, ALPHA2

############## F Function ###############
def F (X, COLUMBIAVK, Z1, Z2, AU,Screening):
    #UNIVERAIAL SCREENING POTENTIAL Screening=0
    if Screening == 0:
        if X == O:
            return 0
        else:
            return COLUMBIAVK * X * (.35 * np.exp(-.3 / X / AU) + .55 * np.exp(-1.2 / X / AU) + .1 * np.exp(-6 / X / AU))
    
    #COLUMBIA POTENTIAL Screening
    else:
        return COLUMBIAVK * X

############## DF Function ###############
def DF(X, COLUMBIAVK, Z1, Z2, AU):
    return F(X, COLUMBIAVK, Z1, Z2, AU) / X + COLUMBIAVK / X * (.35 * EXP(-.3 / X / AU) * .3 / AU + .55 * EXP(-1.2 / X / AU) * 1.2 / AU + .1 * EXP(-6 / X / AU) * 6 / AU)

#########################################################
