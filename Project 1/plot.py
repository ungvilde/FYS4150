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
    xexact, uexact = readfile("data/exact.txt")
    
    # make plot for Problem 2
    plt.figure()
    plt.plot(xexact, uexact)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$u(x)$')
    plt.savefig('figs/problem2.pdf')

    # make plot for Problem 7
    plt.figure()
    plt.plot(xexact, uexact, label="Exact $u(x)$")
    for i in range(1,4):
        N = 10**i
        filename = f"data/general_approx_N{N}.txt"
        x, v = readfile(filename)    
        plt.plot(x, v, label=f"$N=10^{i}$")
    plt.xlabel(r'$x$')
    plt.ylabel(r'$u(x)$')
    plt.legend()
    plt.savefig('figs/problem7.pdf')

    # make plot for Problem 9
    plt.figure()
    plt.plot(xexact, uexact, label="Exact $u(x)$")
    for i in range(1,4):
        N = 10**i
        filename = f"data/special_approx_N{N}.txt"
        x, v = readfile(filename)    
        plt.plot(x, v, label=f"$N=10^{i}$")
    plt.xlabel(r'$x$')
    plt.ylabel(r'$u(x)$')
    plt.legend()
    plt.savefig('figs/problem9.pdf')
    