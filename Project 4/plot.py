import numpy as np
import matplotlib.pyplot as plt

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


energy, magnetisation = readfile("datasets/data_L20_initRand_MC50000_T2.4.txt")

L=20
n_cycles = 50_000
mean_eps = []
mean_m_abs = []
Cv = []
Chi = []

for k in range(1,n_cycles):
    mean_eps.append(np.mean(energy[:k]) / (L*L))
    mean_m_abs.append(np.mean(np.abs(magnetisation[:k])) / (L*L))

plt.title("$T = 2.4$")
plt.plot(np.arange(1,n_cycles), mean_eps)
plt.xlabel("MC cycles")
plt.ylabel(r"$\langle \epsilon \rangle$")
plt.show()

plt.title("$T = 2.4$")
plt.plot(np.arange(1,n_cycles), mean_m_abs)
plt.xlabel("MC cycles")
plt.ylabel(r"$\langle | m | \rangle$")
plt.show()

energy, magnetisation = readfile("datasets/data_L20_initRand_MC50000_T1.txt")

mean_eps = []
mean_m_abs = []
Cv = []
Chi = []

for k in range(1,n_cycles):
    mean_eps.append(np.mean(energy[:k]) / (L*L))
    mean_m_abs.append(np.mean(np.abs(magnetisation[:k])) / (L*L))

plt.plot(np.arange(1,n_cycles), mean_eps)
plt.title("$T = 1.0, Random initiasion$")
plt.xlabel("MC cycles")
plt.ylabel(r"$\langle \epsilon \rangle$")
plt.show()

plt.title("$T = 1.0$")
plt.plot(np.arange(1,n_cycles), mean_m_abs)
plt.xlabel("MC cycles")
plt.ylabel(r"$\langle | m | \rangle$")
plt.show()

energy, magnetisation = readfile("datasets/data_L20_initOrdered_MC50000_T1.txt")

mean_eps = []
mean_m_abs = []
Cv = []
Chi = []

for k in range(1,n_cycles):
    mean_eps.append(np.mean(energy[:k]) / (L*L))
    mean_m_abs.append(np.mean(np.abs(magnetisation[:k])) / (L*L))

plt.plot(np.arange(1,n_cycles), mean_eps)
plt.title("$T = 1.0, Ordered initiasion$")
plt.xlabel("MC cycles")
plt.ylabel(r"$\langle \epsilon \rangle$")
plt.show()

plt.title("$T = 1.0$")
plt.plot(np.arange(1,n_cycles), mean_m_abs)
plt.xlabel("MC cycles")
plt.ylabel(r"$\langle | m | \rangle$")
plt.show()