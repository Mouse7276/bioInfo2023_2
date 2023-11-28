from Bio.Seq import Seq
import numpy as np
import random
import matplotlib.pyplot as plt

import variables
import criteria
import mutation
import crossover

# Create variable
seq = variables.createDNASeq(500)
Pt = variables.createPt(seq)
PtPrime = variables.createPtPrime(Pt)
P1 = variables.getP1(seq, Pt)

# Print result
print('sequence: ', seq, '\n')
print('Pt: ', Pt, '\n')
print('Pt\': ', PtPrime, '\n')
print('P1: ', P1, '\n')
print('length: ', variables.getDNALength(P1), '\n')
print('Tm: ', variables.getDNAMeltingTemp(P1), '\n')
print('GC content: ', variables.getDNAGCContent(P1), '\n')
print('annealing position appearance: ', variables.checkDNAAnnealingPositionAppearance(P1, seq), '\n')
print('3\' end primer: ', variables.checkDNA3EndPrimer(P1), '\n')
