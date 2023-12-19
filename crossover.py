from Bio.Seq import Seq
import numpy as np
import random
import matplotlib.pyplot as plt

import variables
import criteria
import mutation

R = []
i = 0
for i in range(4):
    R.append(random.randint(0,1))

Ptx = [145,164,989,1011]
Pty = [150,173,1066,1091]
X = []
Y = []
X = variables.createPtPrime(Ptx)
Y = variables.createPtPrime(Pty)


def crossoverPt(X,Y):
    Xdf = X
    Ydf = Y
    if R[0] == 1:
        Xdf[0] = Y[0]
        Ydf[0] = X[0]
    elif R[1] == 1:
        Xdf[1] = Y[1]
        Ydf[1] = X[1]
    elif R[2] == 1:
        Xdf[2] = Y[2]
        Ydf[2] = X[2]
    elif R[3] == 1:
        Xdf[3] = Y[3]
        Ydf[3] = X[3]
    return Xdf, Ydf

crossoverPt(X,Y)