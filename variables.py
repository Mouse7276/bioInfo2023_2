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

# Function to create a forward primer index of a given DNA sequence 
def createFwdPrimerIndex (sequence):
    return

# Function to create a reverse primer index of a given DNA sequence
def createRevPrimerIndex (sequence):
    return

# Function to get the length of a DNA sequence
def getDNALength (sequence):
    return

# Function to get the melting temperature of a DNA sequence
def getDNAMeltingTemp (sequence):
    return

# Function to get the GC content of a DNA sequence
def getDNAGCContent (sequence):
    return

# Function to check annealing position appearance of a DNA sequence only one time
def checkDNAAnnealingPositionAppearance (sequence):
    return

# Function to check 3' end of the primer of a DNA sequence is G or C
def checkDNA3EndPrimer (sequence):
    return
