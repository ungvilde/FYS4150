#from cProfile import label
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

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

def u_func(x): 
    #excact solution
    return 1 - (1 - np.exp(-10.))*x - np.exp(-10*x);

if __name__ == '__main__':
    # # load data
    xexact, uexact = readfile("data/exact.txt")
    
    # make plot for Problem 2
    plt.figure()
    plt.plot(xexact, uexact)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$u(x)$')
    plt.savefig('figs/problem2.pdf')

    # make plot for Problem 7
    plt.figure()
    for i in range(1,4):
        N = 10**i
        filename = f"data/general_approx_N{N}.txt"
        x, v = readfile(filename)    
        plt.plot(x, v, label=f"$N=10^{i}$")
    plt.plot(xexact, uexact, 'k--', label="Exact $u(x)$")
    plt.xlabel(r'$x$')
    plt.ylabel(r'$u(x)$')
    plt.legend()
    plt.savefig('figs/problem7.pdf')

    # make plot for Problem 9
    # plt.figure()
    # plt.plot(xexact, uexact, label="Exact $u(x)$")
    # for i in range(1,4):
    #     N = 10**i
    #     filename = f"data/special_approx_N{N}.txt"
    #     x, v = readfile(filename)    
    #     plt.plot(x, v, label=f"$N=10^{i}$")
    # plt.xlabel(r'$x$')
    # plt.ylabel(r'$u(x)$')
    # plt.legend()
    # plt.savefig('figs/problem9.pdf')
    
    # now we make error plots
    plt.figure()
    for i in range(1,6):
        N = 10**i
        filename = f"data/general_approx_N{N}.txt"
        x, v = readfile(filename)    
        x = x[1:-1] # remove boundary values to avoid 0-0 terms
        v = v[1:-1]
        Delta = np.abs(v - u_func(x)) # absolute error
        plt.plot(x, Delta, label=f"$N=N^{i}$")
    plt.yscale('log')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\Delta$')
    plt.legend()
    plt.savefig("figs/problem8a.pdf")

    plt.figure()
    for i in range(1,6):
        N = 10**i
        filename = f"data/general_approx_N{N}.txt"
        x, v = readfile(filename)    
        x = x[1:-1] # remove boundary values to avoid 0/0 terms
        v = v[1:-1]
        epsilon = np.abs((v - u_func(x))/u_func(x)) # absolute error
        plt.plot(x, epsilon, label=f"$N=N^{i}$")
    plt.yscale('log')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\Delta$')
    plt.legend()
    plt.savefig("figs/problem8b.pdf")

    plt.figure()
    eps_max = []
    for i in range(1,8):
        N = 10**i
        filename = f"data/general_approx_N{N}.txt"
        x, v = readfile(filename)   
        x = x[1:-1]
        v = v[1:-1]
        eps = np.abs((v - u_func(x))/u_func(x))
        eps_max.append(np.max(eps))
    plt.loglog([10**i for i in range(1,8)], eps_max, '-o')
    plt.xlabel("Number of steps")
    plt.ylabel(r"$\max(\epsilon)$")
    plt.savefig("figs/problem8c.pdf")

    df = pd.DataFrame({'N' : [f"10{i}" for i in range(1,8)], r'max(epsilon)' : np.round(eps_max, 8)})
    print(df.to_latex(index=False)) 

    # we make the timing table
    logN, time_general, time_special = [], [], []
    with open("data/timing.txt", 'r') as f:
        for line in f:
            logN.append(float(line.split()[0]))
            time_general.append(float(line.split()[1]))
            time_special.append(float(line.split()[2]))
    
    df = pd.DataFrame({
      'N':   [f"10{i}" for i in logN],
      'Time general': time_general,
      'Time special': time_special
    })

    print(df.to_latex(index=False)) 