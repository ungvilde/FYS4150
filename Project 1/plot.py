#from cProfile import label
import numpy as np
import matplotlib.pyplot as plt

def readfile(filename):
    """
    Read data file and and return a table of x and u values.
    Input:
    - filename: name of the file to read
    """
    x, u = [], []
    with open(filename, 'r') as f:
        for line in f:
            x.append(float(line.split()[0]))
            u.append(float(line.split()[1]))

    return np.array(x), np.array(u)

if __name__ == '__main__':
    # load data
    xexact, uexact = readfile("data.txt")
    x1, v1 = readfile("approx10.txt")
    x2, v2 = readfile("approx100.txt")
    x3, v3 = readfile("approx1000.txt")
    x4, v4 = readfile("approx10000.txt")

    # make plot for Problem 2
    plt.figure()
    plt.plot(xexact, uexact, label="Exact $u(x)$")
    plt.xlabel(r'$x$')
    plt.ylabel(r'$u(x)$')
    plt.savefig('figs/problem2.pdf')

    # make plot for Problem 7
    plt.figure()
    plt.plot(xexact, uexact, label="Exact $u(x)$")
    plt.plot(x1, v1, label="$N=10$")
    plt.plot(x2, v2, label="$N=10^2$")
    plt.plot(x3, v3, label="$N=10^3$")
    plt.plot(x4, v4, label="$N=10^4$")
    plt.xlabel(r'$x$')
    plt.ylabel(r'$u(x)$')
    plt.legend()
    plt.savefig('figs/problem7.pdf')
