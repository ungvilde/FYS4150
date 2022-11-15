import numpy as np
import matplotlib.pyplot as plt

from common import *


energy_rand, magnetisation_rand = readfile("datasets/data_L20_initRand_MC50000_T2.4.txt")
energy_ordered, magnetisation_ordered = readfile("datasets/data_L20_initOrdered_MC50000_T2.4.txt")

L=20
# n_cycles = 50_000
# mean_eps = [[],[]]
# mean_m_abs = [[],[]]

# for k in range(1, n_cycles):
#     mean_eps[0].append(np.mean(energy_rand[:k]) / (L*L))
#     mean_eps[1].append(np.mean(energy_ordered[:k]) / (L*L))

#     mean_m_abs[0].append(np.mean(np.abs(magnetisation_rand[:k])) / (L*L))
#     mean_m_abs[1].append(np.mean(np.abs(magnetisation_ordered[:k])) / (L*L))

# plt.figure(figsize=(12*cm,10*cm))
# plt.plot(np.arange(1, n_cycles), mean_eps[0], label="Random init.")
# plt.plot(np.arange(1, n_cycles), mean_eps[1], label="Ordered init.")
# plt.xlabel("MC cycles")
# plt.legend()
# plt.ylabel(r"$\langle \epsilon \rangle$")
# plt.text(x=0.7, y=0.7, s=r"$T=2.4 \, J / k_B$", transform=plt.gca().transAxes)
# plt.tight_layout()
# plt.savefig("figs/eps_MC_cycles_random_ordered_T2.4.pdf")

# plt.figure(figsize=(12*cm,10*cm))
# plt.plot(np.arange(1, n_cycles), mean_m_abs[0], label="Random init.")
# plt.plot(np.arange(1, n_cycles), mean_m_abs[1], label="Ordered init.")
# plt.xlabel("MC cycles")
# plt.legend()
# plt.ylabel(r"$\langle |m| \rangle$")
# plt.text(x=0.7, y=0.7, s=r"$T=2.4 \, J  / k_B$", transform=plt.gca().transAxes)
# plt.tight_layout()
# plt.savefig("figs/m_abs_MC_cycles_random_ordered_T2.4.pdf")

epsilon = energy_rand[10000:] / (L*L)

#plt.hist(epsilon, density=True, bins=20)
#plt.xlabel("$\epsilon$")
#plt.show()
# ----------------------------- #

energy_rand, magnetisation_rand = readfile("datasets/data_L20_initRand_MC50000_T1.txt")
energy_ordered, magnetisation_ordered = readfile("datasets/data_L20_initOrdered_MC50000_T1.txt")

epsilon = energy_rand[10000:] / (L*L)

plt.hist(epsilon, density=True, bins=np.arange(-2.0, -1.94, step=0.01) )
plt.xlabel("$\epsilon$")
plt.show()

# mean_eps = [[],[]]
# mean_m_abs = [[],[]]

# for k in range(1, n_cycles):
#     mean_eps[0].append(np.mean(energy_rand[:k]) / (L*L))
#     mean_eps[1].append(np.mean(energy_ordered[:k]) / (L*L))

#     mean_m_abs[0].append(np.mean(np.abs(magnetisation_rand[:k])) / (L*L))
#     mean_m_abs[1].append(np.mean(np.abs(magnetisation_ordered[:k])) / (L*L))


# plt.figure(figsize=(12*cm,10*cm))
# plt.plot(np.arange(1, n_cycles), mean_eps[0], label="Random init.")
# plt.plot(np.arange(1, n_cycles), mean_eps[1], label="Ordered init.")
# plt.xlabel("MC cycles")
# plt.legend()
# plt.ylabel(r"$\langle \epsilon \rangle$")
# plt.text(x=0.7, y=0.7, s=r"$T=1.0 \, J / k_B$", transform=plt.gca().transAxes)
# plt.tight_layout()
# plt.savefig("figs/eps_MC_cycles_random_ordered_T1.pdf")

# plt.figure(figsize=(12*cm,10*cm))
# plt.plot(np.arange(1, n_cycles), mean_m_abs[0], label="Random init.")
# plt.plot(np.arange(1, n_cycles), mean_m_abs[1], label="Ordered init.")
# plt.xlabel("MC cycles")
# plt.legend()
# plt.ylabel(r"$\langle |m| \rangle$")
# plt.text(x=0.7, y=0.7, s=r"$T=1.0 \, J  / k_B$", transform=plt.gca().transAxes)
# plt.tight_layout()
# plt.savefig("figs/m_abs_MC_cycles_random_ordered_T1.pdf") 