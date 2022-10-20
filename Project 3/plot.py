import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import seaborn as sns

sns.set_theme("notebook", "whitegrid")

from plot_functions import *

# plot relative error of both methods
rel_error_RK4_4 = readfile("data/rel_error_RK4_dt_0.0125.txt")
rel_error_RK4_8 = readfile("data/rel_error_RK4_dt_0.00625.txt")
rel_error_RK4_16 = readfile("data/rel_error_RK4_dt_0.003125.txt")
rel_error_RK4_32 = readfile("data/rel_error_RK4_dt_0.0015625.txt")
errors = [rel_error_RK4_4[1:], rel_error_RK4_8[1:], rel_error_RK4_16[1:], rel_error_RK4_32[1:]]
plot_rel_error(errors, "Runge-Kutta 4th order", "Runge-Kutta 4th order")

rel_error_FE_4 = readfile("data/rel_error_FE_dt_0.0125.txt")
rel_error_FE_8 = readfile("data/rel_error_FE_dt_0.00625.txt")
rel_error_FE_16 = readfile("data/rel_error_FE_dt_0.003125.txt")
rel_error_FE_32 = readfile("data/rel_error_FE_dt_0.0015625.txt")
errors = [rel_error_FE_4[1:], rel_error_FE_8[1:], rel_error_FE_16[1:], rel_error_FE_32[1:]]
plot_rel_error(errors, "Forward Euler", "Forward Euler")

# compute convergence rate
error_RK4_4 = readfile("data/error_RK4_dt_0.0125.txt")
error_RK4_8 = readfile("data/error_RK4_dt_0.00625.txt")
error_RK4_16 = readfile("data/error_RK4_dt_0.003125.txt")
error_RK4_32 = readfile("data/error_RK4_dt_0.0015625.txt")
errors = [error_RK4_4[1:], error_RK4_8[1:], error_RK4_16[1:], error_RK4_32[1:]]
rate_RK4 = convergence_rate(errors)
print("Runge-Kutta 4th order error rate = ", np.round(rate_RK4, 3))

error_FE_4 = readfile("data/error_FE_dt_0.0125.txt")
error_FE_8 = readfile("data/error_FE_dt_0.00625.txt")
error_FE_16 = readfile("data/error_FE_dt_0.003125.txt")
error_FE_32 = readfile("data/error_FE_dt_0.0015625.txt")
errors = [error_FE_4[1:], error_FE_8[1:], error_FE_16[1:], error_FE_32[1:]]
rate_FE = convergence_rate(errors)
print("Forward Euler error rate = ", np.round(rate_FE, 3))

# plot positions with Coloumb interaction
p1_x = readfile("data/x_values_RK4_p1_interaction_1_dt_0.001.txt")
p1_y = readfile("data/y_values_RK4_p1_interaction_1_dt_0.001.txt")
p1_z = readfile("data/z_values_RK4_p1_interaction_1_dt_0.001.txt")

p2_x = readfile("data/x_values_RK4_p2_interaction_1_dt_0.001.txt")
p2_y = readfile("data/y_values_RK4_p2_interaction_1_dt_0.001.txt")
p2_z = readfile("data/z_values_RK4_p2_interaction_1_dt_0.001.txt")

p1_v_x = readfile("data/v_x_values_RK4_p1_interaction_1_dt_0.001.txt")
p1_v_z = readfile("data/v_z_values_RK4_p1_interaction_1_dt_0.001.txt")

p2_v_x = readfile("data/v_x_values_RK4_p2_interaction_1_dt_0.001.txt")
p2_v_z = readfile("data/v_z_values_RK4_p2_interaction_1_dt_0.001.txt")

plot_position_time([p1_z, p2_z], time = np.arange(0, 50, 0.001), 
labels =["Particle 1", "Particle 2"], ylab="$z$ values", title="With interaction", save=True)
#plot_position_time(p1_z, np.arange(0,50,0.001), "Particle 1; $z$-values", "With interaction", save=True)
#plot_position_time(p2_z, np.arange(0,50,0.001), "Particle 2; $z$-values", "With interaction", save=True)

plot_xy_trajectory(p1_x, p1_y, p2_x, p2_y, "Particle 1", "Particle 2", "With interaction", "($T=50\mu$s)", save=True)
plot_phase_space(p1_x, p1_v_x, p2_x, p2_v_x, "Particle 1", "Particle 2", "With interaction", "($T=50\mu$s)", 
"$x$ values", "$v_x$ values", save=True)
plot_phase_space(p1_z, p1_v_z, p2_z, p2_v_z, "Particle 1", "Particle 2", "With interaction", "($T=50\mu$s)", 
"$z$ values", "$v_z$ values", save=True)
plot_xyz_trajectory(p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, "Particle 1", "Particle 2", "With interaction", save=True)

# plot positions without Coloumb interaction
p1_x = readfile("data/x_values_RK4_p1_interaction_0_dt_0.001.txt")
p1_y = readfile("data/y_values_RK4_p1_interaction_0_dt_0.001.txt")
p1_z = readfile("data/z_values_RK4_p1_interaction_0_dt_0.001.txt")

p2_x = readfile("data/x_values_RK4_p2_interaction_0_dt_0.001.txt")
p2_y = readfile("data/y_values_RK4_p2_interaction_0_dt_0.001.txt")
p2_z = readfile("data/z_values_RK4_p2_interaction_0_dt_0.001.txt")

p1_v_x = readfile("data/v_x_values_RK4_p1_interaction_0_dt_0.001.txt")
p1_v_z = readfile("data/v_z_values_RK4_p1_interaction_0_dt_0.001.txt")

p2_v_x = readfile("data/v_x_values_RK4_p2_interaction_0_dt_0.001.txt")
p2_v_z = readfile("data/v_z_values_RK4_p2_interaction_0_dt_0.001.txt")

plot_position_time([p1_z], time = np.arange(0, 50, 0.001), 
labels =["Particle 1"], ylab="$z$ values", title="plot_period", save=True, plot_period=True)

plot_position_time([p1_z, p2_z], time = np.arange(0, 50, 0.001), 
labels =["Particle 1", "Particle 2"], ylab="$z$ values", title="No interaction", save=True, plot_period=False)

plot_xy_trajectory(p1_x, p1_y, p2_x, p2_y, "Particle 1", "Particle 2", "No interaction", "($T=50\mu$s)", save=True)
plot_phase_space(p1_x, p1_v_x, p2_x, p2_v_x, 
"Particle 1", "Particle 2", "No interaction", "($T=50\mu$s)", 
"$x$ values", "$v_x$ values", save=True)
plot_phase_space(p1_z, p1_v_z, p2_z, p2_v_z, "Particle 1", "Particle 2", "No interaction", "($T=50\mu$s)",
"$z$ values", "$v_z$ values",  save=True)
plot_xyz_trajectory(p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, "Particle 1", "Particle 2", "No interaction", save=True)

# plot fraction of particles trapped
# should make function for this...
frac_f1 = readfile("data/fraction_interaction_0_dw_0.02_f_0.1_tot_time_500.txt")
frac_f2 = readfile("data/fraction_interaction_0_dw_0.02_f_0.4_tot_time_500.txt")
frac_f3 = readfile("data/fraction_interaction_0_dw_0.02_f_0.7_tot_time_500.txt")
frequency = []
with open("data/fraction_interaction_0_dw_0.02_f_0.1_tot_time_500.txt", 'r') as f:
    for line in f:
        frequency.append(float(line.split()[1]))

plot_frequency_fraction(frequency, [frac_f1, frac_f2, frac_f3],
labels_list=["$f=0.1$", "$f=0.4$", "$f=0.7$"], title="No interaction")

# plt.figure()
# plt.plot(frequency, frac_f1, label = "$f = 0.1$")
# plt.plot(frequency, frac_f2, label = "$f = 0.4$")
# plt.plot(frequency, frac_f3, label = "$f = 0.7$")
# plt.legend()
# plt.xlabel("Frequency $\omega_V$")
# plt.ylabel("Fraction trapped")
# plt.savefig("figs/fraction_frequency_plot.pdf")

cm = 1/2.54
plt.figure(figsize=(12*cm, 12*cm))
x = readfile("data/x_values_resonance_interaction_1_dt_0.05_wV_1.5_f_0.7.txt")
y = readfile("data/y_values_resonance_interaction_1_dt_0.05_wV_1.5_f_0.7.txt")
z = readfile("data/z_values_resonance_interaction_1_dt_0.05_wV_1.5_f_0.7.txt")

x1 = readfile("data/x_values_resonance_interaction_0_dt_0.05_wV_1.5_f_0.7.txt")
y1 = readfile("data/y_values_resonance_interaction_0_dt_0.05_wV_1.5_f_0.7.txt")
z1 = readfile("data/z_values_resonance_interaction_0_dt_0.05_wV_1.5_f_0.7.txt")

plt.plot(x, y, label = "With interaction")
plt.plot(x1, y1, label = "No interaction")
plt.legend()
plt.xlim([-60,60])
plt.ylim([-60,60])
plt.ylabel("$y$ value ($\mu$m)")
plt.xlabel("$x$ value ($\mu$m)")
plt.text(x=-50, y=-50, s="$\omega_V=1.5, f = 0.7$")
plt.plot([x1[0],x[0]], [y1[0], y[0]], 'ro')
plt.plot([x1[-1], x[-1]], [y1[-1], y[-1]], 'k+')
plt.tight_layout()
plt.savefig("figs/xy_resonance_trajectory.pdf")
plt.close()

plt.figure(figsize=(12*cm, 10*cm))
time = np.arange(0, 500, 0.05)
plt.plot(time, z, label = "With interaction")
plt.plot(time, z1, label = "No interaction")
plt.legend()
plt.xlim([0, 42])
plt.ylim([-600, 600])
plt.xlabel("Time ($\mu$s)")
plt.ylabel("$z$ value ($\mu$m)")
plt.text(x=3, y=-500, s="$\omega_V=1.5, f = 0.7$")
plt.tight_layout()
plt.savefig("figs/z_time_resonance.pdf")
plt.close()

#plot_xyz_trajectory(x[:600], y[:600], z[:600], x1[:600], y1[:600], z1[:600], 
#"With interaction, $T=15 \mu s$", "No interaction, $T=20 \mu s$", "$\omega_V=1.5, f = 0.7$", save=True)

frac_0 = readfile("data/fraction_interaction_0_dw_0.005_f_0.1_tot_time_500.txt")
frac_1 = readfile("data/fraction_interaction_1_dw_0.005_f_0.1_tot_time_500.txt")
frequency = []
with open("data/fraction_interaction_0_dw_0.005_f_0.1_tot_time_500.txt", 'r') as f:
    for line in f:
        frequency.append(float(line.split()[1]))

plot_frequency_fraction(
    frequency, [frac_0, frac_1], 
    labels_list=["No interaction", "With interaction"], title="$f = 0.1$")