import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme("notebook", "whitegrid", palette="colorblind")

from common import *
#sns.set_palette("husl")
params = {
    'legend.fontsize': 9,
    'font.size': 9,
    'figure.figsize': (8.647*cm, 8.0*cm),
    'axes.labelsize': 9,
    'axes.titlesize': 9,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9
    }

plt.rcParams.update(params)

# # plot to esimtate burn time when T=2.4

# # data collected using different initiations
# energy_rand, magnetisation_rand = readfile("datasets/data_L20_initRand_MC50000_T2.4.txt")
# energy_ordered, magnetisation_ordered = readfile("datasets/data_L20_initOrdered_MC50000_T2.4.txt")

# L=20
# n_cycles = 50_000
# mean_eps = [[],[]]
# mean_m_abs = [[],[]]

# for k in range(1, n_cycles):
#     mean_eps[0].append(np.mean(energy_rand[:k]) / (L*L))
#     mean_eps[1].append(np.mean(energy_ordered[:k]) / (L*L))

#     mean_m_abs[0].append(np.mean(np.abs(magnetisation_rand[:k])) / (L*L))
#     mean_m_abs[1].append(np.mean(np.abs(magnetisation_ordered[:k])) / (L*L))

# plt.figure()
# plt.plot(np.arange(1, n_cycles), mean_eps[0], label="Random initial state")
# plt.plot(np.arange(1, n_cycles), mean_eps[1], label="Ordered initial state")
# plt.xlabel("MC cycles")
# plt.legend()
# plt.ylabel(r"$\langle \epsilon \rangle$")
# plt.text(x=0.7, y=0.1, s=r"$T=2.4 \, J / k_B$", transform=plt.gca().transAxes)
# plt.tight_layout()
# plt.savefig("figs/eps_MC_cycles_random_ordered_T2.4.pdf")

# plt.figure()
# plt.plot(np.arange(1, n_cycles), mean_m_abs[0], label="Random init.")
# plt.plot(np.arange(1, n_cycles), mean_m_abs[1], label="Ordered init.")
# plt.xlabel("MC cycles")
# plt.legend()
# plt.ylabel(r"$\langle |m| \rangle$")
# plt.text(x=0.7, y=0.1, s=r"$T=2.4 \, J  / k_B$", transform=plt.gca().transAxes)
# plt.tight_layout()
# plt.savefig("figs/m_abs_MC_cycles_random_ordered_T2.4.pdf")

# ----------------------------------------------- #

# plot to estimate burn time when T=1

# energy_rand, magnetisation_rand = readfile("datasets/data_L20_initRand_MC50000_T1.txt")
# energy_ordered, magnetisation_ordered = readfile("datasets/data_L20_initOrdered_MC50000_T1.txt")

# mean_eps = [[],[]]
# mean_m_abs = [[],[]]

# for k in range(1, n_cycles):
#     mean_eps[0].append(np.mean(energy_rand[:k]) / (L*L))
#     mean_eps[1].append(np.mean(energy_ordered[:k]) / (L*L))

#     mean_m_abs[0].append(np.mean(np.abs(magnetisation_rand[:k])) / (L*L))
#     mean_m_abs[1].append(np.mean(np.abs(magnetisation_ordered[:k])) / (L*L))


# plt.figure()
# plt.plot(np.arange(1, n_cycles), mean_eps[0], label="Random init.")
# plt.plot(np.arange(1, n_cycles), mean_eps[1], label="Ordered init.")
# plt.xlabel("MC cycles")
# plt.legend()
# plt.ylabel(r"$\langle \epsilon \rangle$")
# plt.text(x=0.7, y=0.7, s=r"$T=1.0 \, J / k_B$", transform=plt.gca().transAxes)
# plt.tight_layout()
# plt.savefig("figs/eps_MC_cycles_random_ordered_T1.pdf")

# plt.figure()
# plt.plot(np.arange(1, n_cycles), mean_m_abs[0], label="Random init.")
# plt.plot(np.arange(1, n_cycles), mean_m_abs[1], label="Ordered init.")
# plt.xlabel("MC cycles")
# plt.legend()
# plt.ylabel(r"$\langle |m| \rangle$")
# plt.text(x=0.7, y=0.3, s=r"$T=1.0 \, J  / k_B$", transform=plt.gca().transAxes)
# plt.tight_layout()
# plt.savefig("figs/m_abs_MC_cycles_random_ordered_T1.pdf") 

# ------------------------------------ #

# We make histograms for T=1 and T=2.4

energyT24, _ = readfile("datasets/hist_data_L20_initRand_MC100000_T2.4.txt")
L = 20
N = L*L
epsilon = energyT24 / N

plt.figure()
plt.hist(epsilon, density=True, bins=20)
plt.xlabel("$\epsilon$")
plt.ylabel("Density")
plt.text(
    x=0.1, 
    y=0.8, 
    s=f"$T=2.4 \, J  / k_B$\nVar$(\epsilon) ={np.round(np.var(epsilon), 3)}$",
    transform=plt.gca().transAxes)
plt.tight_layout()

plt.savefig("figs/hist_L20_T2.4_initRand.pdf")

energyT1, _ = readfile("datasets/hist_data_L20_initOrdered_MC100000_T1.txt")
L = 20
N = L*L
epsilon = energyT1 / N

plt.figure()
plt.hist(epsilon, density=True, bins=[-2,-1.98,-1.96,-1.94])
plt.xlabel("$\epsilon$")
plt.ylabel("Density")
plt.text(
    x=0.5, 
    y=0.8, 
    s=f"$T=1.0 \, J  / k_B$\nVar$(\epsilon) = {np.round(np.var(epsilon), 5)}$",
    transform=plt.gca().transAxes)
plt.tight_layout()
plt.savefig("figs/hist_L20_T1_initRand.pdf")
