from Bio.Seq import Seq
import numpy as np
import random
import matplotlib.pyplot as plt

import criteria
import mutation
import crossover

# Definition of the proposed algorithm
# Function to create a random DNA sequence of a given length
def createDNASeq (length):
    choices = ['A', 'T', 'C', 'G']
    seq = ''
    for i in range(length):
        seq += random.choice(choices)
    result = Seq(seq)
    return result

# Funcion to create a Pt of a given DNA sequence
def createPt (sequence):
    forwardPrimerSeqLength = random.randint(10, 30)
    reversePrimerSeqLength = random.randint(10, 30)
    
    forwardPrimerStart = random.randint(0, len(sequence))
    while forwardPrimerStart + forwardPrimerSeqLength > len(sequence):
        forwardPrimerStart = random.randint(0, len(sequence))
    reversePrimerStart = random.randint(forwardPrimerStart + forwardPrimerSeqLength, len(sequence))

    result = [
        forwardPrimerStart,
        forwardPrimerStart + forwardPrimerSeqLength,
        reversePrimerStart,
        reversePrimerStart + reversePrimerSeqLength
    ]
    return result 

# Function to create a Pt' of a given DNA sequence
def createPtPrime (Pt):
    result = [
        Pt[0],
        Pt[1]-Pt[0],
        Pt[2]-Pt[0],
        Pt[3]-Pt[2]
    ]
    return result

# Function to get a P1 of a given DNA sequence
def getP1 (sequence, Pt):
    return sequence[Pt[0]:Pt[1]]

def getP2 (sequence, Pt):
    return sequence[Pt[2]:Pt[3]].reverse_complement()

# Function to get the length of a DNA sequence
def getDNALength (P1):
    return len(P1)

# Function to get the melting temperature of a DNA sequence
def getDNAMeltingTemp (P1):
    GCNumber = P1.count('G') + P1.count('C')
    ATNumber = P1.count('A') + P1.count('T')
    return 4 * GCNumber + 2 * ATNumber

# Function to get the GC content of a DNA sequence
def getDNAGCContent (P1):
    return 100* ((P1.count('G') + P1.count('C')) / len(P1))

# Function to check annealing position appearance of a DNA sequence only one time
def checkDNAAnnealingPositionAppearance (P1, sequence):
    result = 0 
    if sequence.count(P1) > 1:
        result = 1
    return result

# Function to check 3' end of the primer of a DNA sequence is G or C
def checkDNA3EndPrimer (P1):
    result = 1
    if P1[-1] == 'G' or P1[-1] == 'C':
        result = 0
    if P1[-2] == 'G' and P1[-1] == 'C':
        result = 0
    if P1[-2] == 'C' and P1[-1] == 'G':
        result = 0    
    return result