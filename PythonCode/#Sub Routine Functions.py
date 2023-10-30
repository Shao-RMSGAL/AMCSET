#Sub Routine Functions
#Kenneth

import numpy as np
import random as random



############## EMagic Function ###############
def EMagic(X):

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
    
    #COLUMBIA POTENTIAL Screening=0
    else:
        return COLUMBIAVK * X

############## DF Function ###############
def DF(X, COLUMBIAVK, Z1, Z2, AU):
    return F(X, COLUMBIAVK, Z1, Z2, AU) / X + COLUMBIAVK / X * (.35 * EXP(-.3 / X / AU) * .3 / AU + .55 * EXP(-1.2 / X / AU) * 1.2 / AU + .1 * EXP(-6 / X / AU) * 6 / AU)

#########################################################
