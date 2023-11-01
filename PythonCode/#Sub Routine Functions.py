#Sub Routine Functions
#Kenneth

import numpy as np
import random as random



############## EMagic Function ###############
def EMagic(X,ATOMD):
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
def TMagic(MASS1, Z1, MASS2, Z2, INELAB, P):


############## AMagic Function ###############
def AMagic(THETAO, ALPHAO, THETA1RELATIVE, THETA2RELATIVE):


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
