from Bio.Seq import Seq
import numpy as np
import random
import matplotlib.pyplot as plt

import variables

# The proposed algorithm for PCR the primer design
# Function to length check of Pt
def lengthCheckPt (Pt):
    absBf = abs(len(Pt[0]) - len(Pt[1]))
    absBr = abs(len(Pt[2]) - len(Pt[3]))
    if absBf >= 18 and absBr <= 26:
        return 0
    else:
        return 1 

# Function to length difference check of Pt
def lengthDifferenceCheckPt (Pt):
    absBf = abs(len(Pt[0]) - len(Pt[1]))
    absBr = abs(len(Pt[2]) - len(Pt[3]))
    if abs(absBf - absBr) <= 3:
        return 0
    else:
        return 1

# Function to Tmd check of Pt
def TmdCheckPt (Pt):
    absBf = abs(len(Pt[0]) - len(Pt[1]))
    absBr = abs(len(Pt[2]) - len(Pt[3]))
    return 
    

# Function to GC content check of Pt
def GCContentCheckPt (Pt):
    return

# Function to check specificity and termination of primer
def specificityAndTerminationCheckPt (Pt):
    return

# Function to check self complementarity of Pt
def selfComplementarityCheckPt (Pt):
    return

# Function to check cross complementarity of Pt
def crossComplementarityCheckPt (Pt):
    return

# Function to check restirction enzyme site of Pt
def restrictionEnzymeSiteCheckPt (Pt):
    return

# Function to get fitness value of Pt
def getFitnessValuePt (Pt):
    return
