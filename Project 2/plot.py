import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def readfile(filename):
    """
    Read data file and and return a table of x and u values.
    Input:
    - filename: name of the file to read
    """
    N, num_iter = [], []
    with open(filename, 'r') as f:
        for line in f:
            N.append(float(line.split()[0]))
            num_iter.append(float(line.split()[1]))

    return np.array(N), np.array(num_iter)

N, num_iter = readfile("data/problem5.txt")
cm = 1/2.54
fig = plt.figure(figsize=(12*cm,10*cm))
ax = fig.add_subplot(111)

plt.loglog(N, num_iter, '-o')
ax.set_xscale('log', base=2)
ax.set_yscale('log', base=2)
ax.set_xlabel("$N$")
ax.set_ylabel("Iterations")
ax.grid()
fig.tight_layout()
plt.savefig("figs/problem5.pdf")