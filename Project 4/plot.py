import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression

sns.set_theme("notebook", "whitegrid", palette="colorblind")

from common import *

#print(plt.rcParams.keys())
params = {
    'legend.fontsize': 9,
    'font.size': 9,
    'figure.figsize': (8.647*cm, 8.0*cm),
    'axes.labelsize': 9,
    'axes.titlesize': 9,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'lines.markersize': 3.0,
    'lines.linewidth': 1.0,
    }

plt.rcParams.update(params)


# -------------------------------------#
# plot to esimtate burn time when T=2.4

# data collected using different initiations
energy_rand, magnetisation_rand = readfile("datasets/data_L20_initRand_MC50000_T2.4.txt")
energy_ordered, magnetisation_ordered = readfile("datasets/data_L20_initOrdered_MC50000_T2.4.txt")

L=20
n_cycles = 50_000
mean_eps = [[],[]]
mean_m_abs = [[],[]]

for k in range(1, n_cycles):
    mean_eps[0].append(np.mean(energy_rand[:k]) / (L*L))
    mean_eps[1].append(np.mean(energy_ordered[:k]) / (L*L))

    mean_m_abs[0].append(np.mean(np.abs(magnetisation_rand[:k])) / (L*L))
    mean_m_abs[1].append(np.mean(np.abs(magnetisation_ordered[:k])) / (L*L))

plt.figure()
plt.plot(np.arange(1, n_cycles), mean_eps[0], '-',label="Random initial state")
plt.plot(np.arange(1, n_cycles), mean_eps[1], '-',label="Ordered initial state")
plt.xlabel("MC cycles")
plt.legend()
plt.ylabel(r"$\langle \epsilon \rangle$")
plt.text(x=0.7, y=0.1, s=r"$T=2.4 \, J / k_B$", transform=plt.gca().transAxes)
plt.tight_layout()
plt.savefig("figs/eps_MC_cycles_random_ordered_T2.4.pdf")

plt.figure()
plt.plot(np.arange(1, n_cycles), mean_m_abs[0], label="Random initial state")
plt.plot(np.arange(1, n_cycles), mean_m_abs[1], label="Ordered initial state")
plt.xlabel("MC cycles")
plt.legend()
plt.ylabel(r"$\langle |m| \rangle$")
plt.text(x=0.7, y=0.1, s=r"$T=2.4 \, J  / k_B$", transform=plt.gca().transAxes)
plt.tight_layout()
plt.savefig("figs/m_abs_MC_cycles_random_ordered_T2.4.pdf")

# ----------------------------------------------- #

# plot to estimate burn time when T=1

energy_rand, magnetisation_rand = readfile("datasets/data_L20_initRand_MC50000_T1.txt")
energy_ordered, magnetisation_ordered = readfile("datasets/data_L20_initOrdered_MC50000_T1.txt")

mean_eps = [[],[]]
mean_m_abs = [[],[]]

for k in range(1, n_cycles):
    mean_eps[0].append(np.mean(energy_rand[:k]) / (L*L))
    mean_eps[1].append(np.mean(energy_ordered[:k]) / (L*L))

    mean_m_abs[0].append(np.mean(np.abs(magnetisation_rand[:k])) / (L*L))
    mean_m_abs[1].append(np.mean(np.abs(magnetisation_ordered[:k])) / (L*L))


plt.figure()
plt.plot(np.arange(1, n_cycles), mean_eps[0], label="Random initial state")
plt.plot(np.arange(1, n_cycles), mean_eps[1], label="Ordered initial state")
plt.xlabel("MC cycles")
plt.legend()
plt.ylabel(r"$\langle \epsilon \rangle$")
plt.text(x=0.7, y=0.7, s=r"$T=1.0 \, J / k_B$", transform=plt.gca().transAxes)
plt.tight_layout()
plt.savefig("figs/eps_MC_cycles_random_ordered_T1.pdf")

plt.figure()
plt.plot(np.arange(1, n_cycles), mean_m_abs[0], label="Random initial state")
plt.plot(np.arange(1, n_cycles), mean_m_abs[1], label="Ordered initial state")
plt.xlabel("MC cycles")
plt.legend()
plt.ylabel(r"$\langle |m| \rangle$")
plt.text(x=0.7, y=0.3, s=r"$T=1.0 \, J  / k_B$", transform=plt.gca().transAxes)
plt.tight_layout()
plt.savefig("figs/m_abs_MC_cycles_random_ordered_T1.pdf") 

# ------------------------------------ #

# We make histograms for T=1 and T=2.4

energyT10, _ = readfile("datasets/hist_data_L20_initRand_MC100000_T10.txt")
L = 20
N = L*L
epsilon = energyT10 / N

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

plt.savefig("figs/hist_L20_T10.0_initRand.pdf")
print("T=10.0")
print(np.var(epsilon))
print(np.mean(epsilon))

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
print("T=2.4")
print(np.var(epsilon))
print(np.mean(epsilon))

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
print("T=1")
print(np.var(epsilon))
print(np.mean(epsilon))

##------------------------------##
temperatures1 = np.round(np.linspace(2.1, 2.4, 40), 5)

temperatures = np.sort(temperatures1)
print(temperatures)
epsilon_values = [[] for _ in range(4)]
m_abs_values = [[] for _ in range(4)]
heat_capacity_values = [[] for _ in range(4)]
susceptibility_values = [[] for _ in range(4)]

epsilon_values_err = [[] for _ in range(4)]
m_abs_values_err = [[] for _ in range(4)]
heat_capacity_values_err = [[] for _ in range(4)]
susceptibility_values_err = [[] for _ in range(4)]

for i, L in enumerate([40, 60, 80, 100]):
    for temp in temperatures:
        
        energy, magnetisation = readfile(f"datasets/data_parallel_L{L}initRand_MC1000000T{temp}.txt")

        m_abs_values[i].append(compute_m_abs(magnetisation, L=L))
        epsilon_values[i].append(compute_epsilon(energy, L=L))
        heat_capacity_values[i].append(compute_heat_capacity(energy, T=temp, L=L))
        susceptibility_values[i].append(compute_susceptibility(magnetisation, T=temp, L=L))

        m_abs_values_err[i].append(compute_m_abs_error(magnetisation, L=L))
        epsilon_values_err[i].append(compute_epsilon_error(energy, L=L))


plt.figure()
plt.plot(temperatures, epsilon_values[0], '-o', label="$L=40$", markersize=3)
plt.plot(temperatures, epsilon_values[1], '-o', label="$L=60$", markersize=3)
plt.plot(temperatures, epsilon_values[2], '-o', label="$L=80$", markersize=3)
plt.plot(temperatures, epsilon_values[3], '-o', label="$L=100$", markersize=3)
plt.xlabel("Temperature $(J/k_B)$")
plt.ylabel(r"$\langle\epsilon\rangle$")
plt.legend()
plt.tight_layout()
plt.savefig("figs/epsilon_temperature.pdf")

plt.figure()
plt.plot(temperatures, m_abs_values[0], '-o', label="$L=40$", markersize=3)
plt.plot(temperatures, m_abs_values[1], '-o', label="$L=60$", markersize=3)
plt.plot(temperatures, m_abs_values[2], '-o', label="$L=80$", markersize=3)
plt.plot(temperatures, m_abs_values[3], '-o', label="$L=100$", markersize=3)
plt.xlabel("Temperature $(J/k_B)$")
plt.ylabel(r"$\langle|m|\rangle$")
plt.tight_layout()
plt.legend()
plt.savefig("figs/m_abs_temperature.pdf")

plt.figure()
plt.plot(temperatures, heat_capacity_values[0], '-o', label="$L=40$", markersize=3)
plt.plot(temperatures, heat_capacity_values[1], '-o', label="$L=60$", markersize=3)
plt.plot(temperatures, heat_capacity_values[2], '-o', label="$L=80$", markersize=3)
plt.plot(temperatures, heat_capacity_values[3], '-o', label="$L=100$", markersize=3)
plt.xlabel("Temperature $(J/k_B)$")
plt.ylabel(r"$C_V$")
plt.legend()
plt.tight_layout()
plt.savefig("figs/heat_capacity_temperature.pdf")

plt.figure()
plt.plot(temperatures, susceptibility_values[0], '-o', label="$L=40$", markersize=3)
plt.plot(temperatures, susceptibility_values[1], '-o', label="$L=60$", markersize=3)
plt.plot(temperatures, susceptibility_values[2], '-o', label="$L=80$", markersize=3)
plt.plot(temperatures, susceptibility_values[3], '-o', label="$L=100$", markersize=3)
plt.xlabel("Temperature $(J/k_B)$")
plt.ylabel(r"$\chi$")
plt.legend()
plt.tight_layout()
plt.savefig("figs/susceptibility_temperature.pdf")

Tc_values = [[] for _ in range(3)]
for i, L in enumerate([40, 60, 80, 100]):
    m = np.argmax(susceptibility_values[i])
    n = np.argmax(heat_capacity_values[i])
    print("L = ", L)
    print("T (combined) = ", (temperatures[m] + temperatures[n]) / 2)
    print("Tc (heat cap.)= ", temperatures[n])
    print("Tc (sus.)= ", temperatures[m])

    Tc_values[0].append((temperatures[m] + temperatures[n]) / 2)
    Tc_values[1].append(temperatures[m])
    Tc_values[2].append(temperatures[n])


L = np.array([40, 60, 80, 100])
Tc_values0 = np.array(Tc_values[0])
reg = LinearRegression()
reg.fit(1/L.reshape(-1, 1), Tc_values0)
regline = reg.predict(1/L.reshape(-1, 1))
print(reg.coef_)
print(reg.intercept_)
plt.figure()
plt.plot(1/L, Tc_values0, 'o', label = "Data")
plt.plot(1/L, regline, 'k-', label="Regression line")
plt.xlabel("$L^{-1}$")
plt.ylabel("$T_c(L)$")
plt.text(
    x=0.5, 
    y=0.15, 
    s=f"$T_c(\infty)={np.round(reg.intercept_, 3)} \, J  / k_B$",
    transform=plt.gca().transAxes)
plt.tight_layout()
plt.savefig("figs/Tc_regression_combined.pdf")

Tc_values1 = np.array(Tc_values[1])
reg = LinearRegression()
reg.fit(1/L.reshape(-1, 1), Tc_values1)
regline = reg.predict(1/L.reshape(-1, 1))
print(reg.coef_)
print(reg.intercept_)

plt.figure()
plt.plot(1/L, Tc_values1, 'o', label = "Data")
plt.plot(1/L, regline, 'k--', label="Regression line")
plt.xlabel("$L^{-1}$")
plt.ylabel("$T_c(L)$")
plt.text(
    x=0.4, 
    y=0.15, 
    s=f"$T_c(L=\infty)={np.round(reg.intercept_, 3)} \, J  / k_B$",
    transform=plt.gca().transAxes)
plt.tight_layout()
plt.savefig("figs/Tc_regression_susceptibility.pdf")
StdErr = np.var(Tc_values1)


Tc_values2 = np.array(Tc_values[2])
reg = LinearRegression()
reg.fit(1/L.reshape(-1, 1), Tc_values2)
regline = reg.predict(1/L.reshape(-1, 1))
print(reg.coef_)
print(reg.intercept_)
plt.figure()
plt.plot(1/L, Tc_values2, 'o', label = "Data")
plt.plot(1/L, regline, 'k-', label="Regression line")
plt.xlabel("$L^{-1}$")
plt.ylabel("T_c(L)")
plt.tight_layout()
plt.savefig("figs/Tc_regression_heat_capacity.pdf")
