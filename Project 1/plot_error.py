import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from plot import readfile

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

plt.figure()
for Nsteps in [10**i for i in range(1,5)]:
    x, v = readfile(f"approx{Nsteps}.txt")
    x = x[1:-1] # remove boundary values
    v = v[1:-1]
    Delta = np.abs(v - u_func(x))
    plt.plot(x, Delta, label=f"N={Nsteps}")
plt.yscale('log')
plt.xlabel(r'$x$')
plt.ylabel(r'$\Delta$')
plt.legend()
plt.savefig("figs/problem8a.pdf")

plt.figure()
for Nsteps in [10**i for i in range(1,5)]:
    x, v = readfile(f"approx{Nsteps}.txt")
    x = x[1:-1]
    v = v[1:-1]
    eps = np.abs((v - u_func(x))/u_func(x))

    plt.plot(x, eps, label=f"N={Nsteps}")
plt.yscale('log')
plt.legend()
plt.xlabel(r'$x$')
plt.ylabel(r'$\epsilon$')
plt.savefig("figs/problem8b.pdf")

plt.figure()
eps_max = []
for Nsteps in [10**i for i in range(1,8)]:
    x, v = readfile(f"approx{Nsteps}.txt")
    x = x[1:-1]
    v = v[1:-1]
    eps = np.abs((v - u_func(x))/u_func(x))
    eps_max.append(np.max(eps))

plt.loglog([10**i for i in range(1,8)], eps_max, '-o')
plt.xlabel("Number of steps")
plt.ylabel(r"$\max(\epsilon)$")
plt.savefig("figs/problem8c.pdf")

df = pd.DataFrame({'N' : [10**i for i in range(1,8)], r'max(epsilon)' : np.round(eps_max, 8)})
print(df.to_latex(index=False)) 
