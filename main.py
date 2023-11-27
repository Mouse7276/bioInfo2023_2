from Bio.Seq import Seq
import numpy as np
import random
import matplotlib.pyplot as plt

import variables
import criteria
import mutation
import crossover

# Create variables
Gd = variables.createDNASeq(60)
Gd_ = Gd.complement()

# Test functions
print("Gd: " + str(Gd) + "\n")
print("Gd_: " + str(Gd_) + "\n")
