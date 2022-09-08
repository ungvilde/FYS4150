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


def plot(x, u, save=False):
    """Plot an array of x and u values."""
    plt.plot(x, u)

    plt.xlabel(r'$x$')
    plt.ylabel(r'$u(x)$')
    if save:
        plt.savefig('problem_1.pdf')
    plt.show()


if __name__ == '__main__':
    #x, u = readfile("data.txt")
    x, v = readfile("approx.txt")
    h = 1/100
    plot(x, v)