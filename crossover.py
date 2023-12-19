from Bio.Seq import Seq
import numpy as np
import random
import matplotlib.pyplot as plt

import variables

def crossoverPt(PtX,PtY):
    # 어느 자리를 바꿀지 랜덤으로 정함
    R = []
    i = 0
    for i in range(4):
        R.append(random.randint(0,1))

    # Pt' 생성
    Xdf = variables.createPtPrime(Ptx)
    Ydf = variables.createPtPrime(Pty)
    oldXdf = variables.createPtPrime(Ptx)
    oldYdf = variables.createPtPrime(Pty)

    # Crossover 수행 
    if R[0] == 1:
        Xdf[0] = oldYdf[0]
        Ydf[0] = oldXdf[0]
    elif R[1] == 1:
        Xdf[1] = oldYdf[1]
        Ydf[1] = oldXdf[1]
    elif R[2] == 1:
        Xdf[2] = oldYdf[2]
        Ydf[2] = oldXdf[2]
    elif R[3] == 1:
        Xdf[3] = oldYdf[3]
        Ydf[3] = oldXdf[3]

    # 결과 반환     
    return Xdf, Ydf
