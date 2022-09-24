import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def readfile(filename):

    N, num_iter = [], []
    with open(filename, 'r') as f:
        for line in f:
            N.append(float(line.split()[0]))
            num_iter.append(float(line.split()[1]))

    return np.array(N), np.array(num_iter)

def readfile1(filename):

    v1, v2, v3 = [], [], []
    with open(filename, 'r') as f:
        for line in f:
            v1.append(float(line.split()[0]))
            v2.append(float(line.split()[1]))
            v3.append(float(line.split()[2]))

    return v1, v2, v3

def readfile2(filename):

    eigenvals = []
    with open(filename, 'r') as f:
        for line in f:
            eigenvals.append(float(line.split()[0]))

    return eigenvals

# for creating plots with meaningful figsize
cm = 1/2.54

# problem 5
N, num_iter = readfile("data/problem5.txt")
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

# problem 6a
N = 10
xhat = np.linspace(0, 1, N+2) # include boundary points
v1, v2, v3 = readfile1("data/problem6_numerical_N10.txt")
a1, a2, a3 = readfile1("data/problem6_analytic_N10.txt")
eigenvals = readfile2("data/problem6_numerical_N10_eigenvals.txt")
for v in [v1, v2, v3]:
    v.insert(0, 0) #include boundary conditions
    v.append(0)

for a in [a1, a2, a3]:
    a.insert(0, 0) #include boundary conditions
    a.append(0)

fig, ax = plt.subplots(figsize=(18*cm, 6*cm), nrows=1, ncols=3, sharex=True, sharey=True)
ax[0].plot(xhat, a1, label = "$u(\hat x_i)$")
ax[0].plot(xhat, v1, '--', label = "$v_i$")
ax[0].hlines(y=0,xmin=0, xmax=1, linestyles='dotted', colors="black", linewidths=1)

ax[1].plot(xhat, a2, label = "$u(\hat x_i)$")
ax[1].plot(xhat, v2, '--', label = "$v_i$")
ax[1].hlines(y=0,xmin=0,xmax=1, linestyles='dotted', colors="black", linewidths=1)
ax[2].plot(xhat, a3, label = "$u(\hat x_i)$")
ax[2].plot(xhat, v3, '--', label = "$v_i$")
ax[2].hlines(y=0,xmin=0,xmax=1, linestyles='dotted', colors="black", linewidths=1)

for i in range(3):
    ax[i].set_title(f"$\lambda_{i+1}$ = {np.round(eigenvals[i])}")

ax[1].set_xlabel("$\hat x$")
ax[0].set_ylabel("$u$, $v$")
ax[0].legend()
#ax[1].legend()
#ax[2].legend()
plt.tight_layout()
plt.savefig("figs/problem6_N10.pdf")

# problem 6b
N = 100
xhat = np.linspace(0, 1, N+2) 
v1, v2, v3 = readfile1("data/problem6_numerical_N100.txt")
a1, a2, a3 = readfile1("data/problem6_analytic_N100.txt")
eigenvals = readfile2("data/problem6_numerical_N100_eigenvals.txt")

for v in [v1, v2, v3]:
    v.insert(0, 0) #include boundary conditions
    v.append(0)

for a in [a1, a2, a3]:
    a.insert(0, 0) #include boundary conditions
    a.append(0)

fig, ax = plt.subplots(figsize=(18*cm, 6*cm), nrows=1, ncols=3, sharex=True, sharey=True)
ax[0].plot(xhat, a1, label = "$u(\hat x_i)$")
ax[0].plot(xhat, v1, '--', label = "$v_i$")
ax[0].hlines(y=0,xmin=0, xmax=1, linestyles='dotted', colors="black", linewidths=1)

ax[1].plot(xhat, a2, label = "$u(\hat x_i)$")
ax[1].plot(xhat, v2, '--', label = "$v_i$")
ax[1].hlines(y=0,xmin=0,xmax=1, linestyles='dotted', colors="black", linewidths=1)
ax[2].plot(xhat, a3, label = "$u(\hat x_i)$")
ax[2].plot(xhat, v3, '--', label = "$v_i$")
ax[2].hlines(y=0,xmin=0,xmax=1, linestyles='dotted', colors="black", linewidths=1)

for i in range(3):
    ax[i].set_title(f"$\lambda_{i+1}$ = {np.round(eigenvals[i])}")

ax[1].set_xlabel("$\hat x$")
ax[0].set_ylabel("$u$, $v$")
ax[0].legend()
#ax[1].legend()
#ax[2].legend()
plt.tight_layout()
plt.savefig("figs/problem6_N100.pdf")