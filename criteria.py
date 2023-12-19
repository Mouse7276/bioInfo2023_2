from Bio.Seq import Seq
import numpy as np
import random
import matplotlib.pyplot as plt

import variables
import mutation
import crossover

# The proposed algorithm for PCR the primer design
# Function to length check of Pt
def lengthCheckPt (Pt):
    lenBf = Pt[1] - Pt[0]
    lenBr = Pt[3] - Pt[2]
    flag = (18<=lenBf) and (lenBf<=26) and (18<=lenBr) and (lenBf<=26)
    return int(not flag)

# Function to length difference check of Pt
def lengthDifferenceCheckPt (Pt):
    lenBf = Pt[1] - Pt[0]
    lenBr = Pt[3] - Pt[2]
    flag = (abs(lenBf - lenBr) <= 3)
    return int(not flag)

# Function to Tmd check of Pt
def TmdCheckPt (sequence, Pt):
    Bfseq = variables.getP1(sequence, Pt)
    Brseq = variables.getP2(sequence, Pt)
    TmBf = variables.getDNAMeltingTemp(Bfseq)
    TmBr = variables.getDNAMeltingTemp(Brseq)
    flag = (abs(TmBf - TmBr) <= 5)
    return int(not flag)

# Function to GC content check of Pt
def GCContentCheckPt (sequence, Pt):
    Bfseq = variables.getP1(sequence, Pt)
    Brseq = variables.getP2(sequence, Pt)
    GCBf = variables.getDNAGCContent(Bfseq)
    GCBr = variables.getDNAGCContent(Brseq)
    flag = (40 <= GCBf) and (GCBf <= 60) and (40 <= GCBf) and (GCBr <= 60)
    return int(not flag)

# Function to check specificity and termination of primer
def specificityCheckPt (sequence, Pt):
    Bfseq = variables.getP1(sequence, Pt)
    Brseq = variables.getP2(sequence, Pt)
    UniBf = variables.checkDNAAnnealingPositionAppearance(Bfseq)
    UniBr = variables.checkDNAAnnealingPositionAppearance(Brseq)
    return UniBf + UniBr

def TerminationCheckPt (sequence, Pt):
    Bfseq = variables.getP1(sequence, Pt)
    Brseq = variables.getP2(sequence, Pt)
    TermBf = variables.checkDNA3EndPrimer(Bfseq)
    TermBr = variables.checkDNA3EndPrimer(Brseq)
    return TermBf + TermBr

# Function to check self complementarity of Pt
def selfComplementarityCheckPt (sequence, Pt):
    # hairpin check
    Bfseq = variables.getP1(sequence, Pt)
    Brseq = variables.getP2(sequence, Pt)
    for i in range(3, len(Bfseq)//2):
        tempseq1 = Bfseq[:i]
        tempseq2 = Bfseq[i:2*i].reverse_complement()
        diff = [ai==bi for ai, bi in zip(tempseq1, tempseq2)]
        if sum(diff) >= 3:
            return 1
    for i in range(3, len(Brseq)//2):
        tempseq1 = Brseq[:i]
        tempseq2 = Brseq[i:2*i].reverse_complement()
        diff = [ai==bi for ai, bi in zip(tempseq1, tempseq2)]
        if sum(diff) >= 3:
            return 1
    # self-self check
    tempseq = Bfseq.reverse()
    Bfcompseq = Bfseq.reverse_complement()
    Brcompseq = variables.getP1(sequence, Pt[2:])
    for i in range(3, len(Bfseq)):
        tempseq1 = Bfseq[len(Bfseq)-i:]
        tempseq2 = Bfcompseq[:i]
        diff = [ai==bi for ai, bi in zip(tempseq1, tempseq2)]
        if sum(diff) >= 3:
            return 1
    for i in range(3, len(Brseq)):
        tempseq1 = Brseq[len(Brseq)-i:]
        tempseq2 = Brcompseq[:i]
        diff = [ai==bi for ai, bi in zip(tempseq1, tempseq2)]
        if sum(diff) >= 3:
            return 1
    return 0

# Function to check cross complementarity of Pt
def crossComplementarityCheckPt (sequence, Pt):
    Bfseq = variables.getP1(sequence, Pt)
    Brseq = variables.getP1(sequence, Pt[2:])
    for i in range(3, min(len(Bfseq), len(Brseq))):
        tempseq1 = Bfseq[len(Bfseq)-i:]
        tempseq2 = Brseq[:i]
        diff = [ai==bi for ai, bi in zip(tempseq1, tempseq2)]
        if sum(diff) >= 3:
            return 1
    return 0

# Function to check restirction enzyme site of Pt
def restrictionEnzymeSiteCheckPt (sequence, Pt):
    ApaI = 'GGGCCC'
    AvrII = 'CCTAGG'
    BamHI = 'GGATCC'
    BglII = 'AGATCT'
    DraI = 'TTTAAA'
    enz_list = [ApaI, AvrII, BamHI, BglII, DraI]
    Bfseq = variables.getP1(sequence, Pt)
    Brseq = variables.getP2(sequence, Pt)
    for enz in enz_list:
        # Bfseq search
        for i in range(len(Bfseq) - 6):
            tempseq = Bfseq[i:i+6]
            diff = [ai==bi for ai, bi in zip(tempseq, enz)]
            Pm = sum(diff)
            if 3 <= Pm and Pm <= 6:
                return 0
        for i in range(len(Brseq) - 6):
            tempseq = Brseq[i:i+6]
            diff = [ai==bi for ai, bi in zip(tempseq, enz)]
            Pm = sum(diff)
            if 3 <= Pm and Pm <= 6:
                return 0
    return 1

# Function to get fitness value of Pt
def getFitnessValuePt (sequence, Pt):
    Fitness = lengthCheckPt(Pt) + 3*lengthDifferenceCheckPt(Pt) + 3*TmdCheckPt(sequence, Pt) + 3*GCContentCheckPt(sequence, Pt) + 3*TerminationCheckPt(sequence, Pt) + 50*specificityCheckPt(sequence, Pt) + 10*selfComplementarityCheckPt(sequence, Pt) + 10*crossComplementarityCheckPt(sequence, Pt) + restrictionEnzymeSiteCheckPt(sequence, Pt)
    return Fitness
