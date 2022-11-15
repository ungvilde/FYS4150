import numpy as np
import matplotlib.pyplot as plt


cm = 1 / 2.54

def readfile(filename):

    energy = []
    magnetisation = []

    with open(filename, 'r') as f:
        for line in f:
            energy.append(float(line.split()[0]))
            magnetisation.append(float(line.split()[1]))

    return np.array(energy), np.array(magnetisation) 

def heat_capacity(energy, T, L, kb=1):
    return 1 / (L*L) * 1 / (kb * T*T) * (np.var(energy))