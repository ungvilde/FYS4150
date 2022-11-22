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

def compute_heat_capacity(energy, T, L, kb=1):
    return 1 / (L*L) * 1 / (kb * T*T) * (np.var(energy))

def compute_m_abs(magnetisation, L):
    return np.mean(np.abs(magnetisation)) / (L*L)

def compute_epsilon(energy, L):
    return np.mean(energy) / (L*L)

def compute_susceptibility(magnetisation, T, L, kb=1):
    return 1 / (L*L) * 1 / (kb * T) * np.var(np.abs(magnetisation))

def compute_m_abs_error(magnetisation, L):
    return np.std(np.abs(magnetisation)/ (L*L)) 

def compute_epsilon_error(energy, L):
    return np.std(energy/ (L*L)) 

def readfile_num_vals(filename):

    energy = []
    magnetisation = []
    heat_capacity = []
    susceptibility =[]

    with open(filename, 'r') as f:
        for line in f:
            energy.append(float(line.split()[0]))
            heat_capacity.append(float(line.split()[1]))
            magnetisation.append(float(line.split()[2]))
            susceptibility.append(float(line.split()[3]))

    return energy, magnetisation, heat_capacity, susceptibility
