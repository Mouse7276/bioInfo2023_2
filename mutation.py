from Bio.Seq import Seq
import numpy as np
import random
import matplotlib.pyplot as plt

# Functions to make a mutation
def makeMutation(array):
    oldArray = array.copy()

    # Select two random rows in array
    row1Index = random.randint(0, len(array)-1)
    row2Index = random.randint(0, len(array)-1)
    while row1Index == row2Index:
        row2Index = random.randint(0, len(array)-1)
    row1 = array[row1Index].copy()
    row2 = array[row2Index].copy()

    # Select random mutation components
    mutationIndex = random.randint(0, 3)
    newValue = round((row1[mutationIndex] + row2[mutationIndex])/2)

    row1[mutationIndex] = newValue
    row2[mutationIndex] = newValue

    # Return the result
    oldArray[row1Index] = row1
    oldArray[row2Index] = row2

    return oldArray

array = [
    [1, 2, 3, 4],
    [5, 6, 7, 8],
    [9, 10, 11, 12],
    [13, 14, 15, 16]
]

newArray = makeMutation(array)

print(newArray)
