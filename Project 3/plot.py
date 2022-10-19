import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


from plot_functions import *

# plot relative error of both methods
rel_error_RK4_4 = readfile("data/rel_error_RK4_dt_0.0125.txt")
rel_error_RK4_8 = readfile("data/rel_error_RK4_dt_0.00625.txt")
rel_error_RK4_16 = readfile("data/rel_error_RK4_dt_0.003125.txt")
rel_error_RK4_32 = readfile("data/rel_error_RK4_dt_0.0015625.txt")
errors = [rel_error_RK4_4[1:], rel_error_RK4_8[1:], rel_error_RK4_16[1:], rel_error_RK4_32[1:]]
plot_rel_error(errors, "Runge-Kutta 4th order")

rel_error_FE_4 = readfile("data/rel_error_FE_dt_0.0125.txt")
rel_error_FE_8 = readfile("data/rel_error_FE_dt_0.00625.txt")
rel_error_FE_16 = readfile("data/rel_error_FE_dt_0.003125.txt")
rel_error_FE_32 = readfile("data/rel_error_FE_dt_0.0015625.txt")
errors = [rel_error_FE_4[1:], rel_error_FE_8[1:], rel_error_FE_16[1:], rel_error_FE_32[1:]]
plot_rel_error(errors, "Forward Euler")

# compute convergence rate
error_RK4_4 = readfile("data/error_RK4_dt_0.0125.txt")
error_RK4_8 = readfile("data/error_RK4_dt_0.00625.txt")
error_RK4_16 = readfile("data/error_RK4_dt_0.003125.txt")
error_RK4_32 = readfile("data/error_RK4_dt_0.0015625.txt")
errors = [error_RK4_4[1:], error_RK4_8[1:], error_RK4_16[1:], error_RK4_32[1:]]
rate_RK4 = convergence_rate(errors)
print("Runge-Kutta 4th order error rate = ", rate_RK4)

error_FE_4 = readfile("data/error_FE_dt_0.0125.txt")
error_FE_8 = readfile("data/error_FE_dt_0.00625.txt")
error_FE_16 = readfile("data/error_FE_dt_0.003125.txt")
error_FE_32 = readfile("data/error_FE_dt_0.0015625.txt")
errors = [error_FE_4[1:], error_FE_8[1:], error_FE_16[1:], error_FE_32[1:]]
rate_FE = convergence_rate(errors)
print("Forward Euler error rate = ", rate_FE)

# plot positions with Coloumb interaction
p1_x = readfile("data/x_values_RK4_p1_interaction_1_dt_0.001.txt")
p1_y = readfile("data/y_values_RK4_p1_interaction_1_dt_0.001.txt")
p1_z = readfile("data/z_values_RK4_p1_interaction_1_dt_0.001.txt")

p2_x = readfile("data/x_values_RK4_p2_interaction_1_dt_0.001.txt")
p2_y = readfile("data/y_values_RK4_p2_interaction_1_dt_0.001.txt")
p2_z = readfile("data/z_values_RK4_p2_interaction_1_dt_0.001.txt")

p1_v_x = readfile("data/v_x_values_RK4_p1_interaction_0_dt_0.001.txt")
p1_v_z = readfile("data/v_z_values_RK4_p1_interaction_0_dt_0.001.txt")

p2_v_x = readfile("data/v_x_values_RK4_p2_interaction_0_dt_0.001.txt")
p2_v_z = readfile("data/v_z_values_RK4_p2_interaction_0_dt_0.001.txt")

plot_position_time(p1_z, np.arange(0,50,0.001), "Particle 1; $z$-values", "With interaction", save=True)
plot_position_time(p2_z, np.arange(0,50,0.001), "Particle 2; $z$-values", "With interaction", save=True)

plot_xy_trajectory(p1_x, p1_y, p2_x, p2_y, "Particle 1", "Particle 2", "With interaction", "($T=50\mu$s)", save=True)
plot_phase_space(p1_x, p1_v_x, p2_x, p2_v_x, "Particle 1", "Particle 2", "With interaction", "($T=50\mu$s)", save=True)
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

plot_position_time(p1_z, np.arange(0, 50, 0.001), "Particle 1; $z$-values", "No interaction", save=True)
plot_position_time(p2_z, np.arange(0, 50, 0.001), "Particle 2; $z$-values", "No interaction", save=True)

plot_xy_trajectory(p1_x, p1_y, p2_x, p2_y, "Particle 1", "Particle 2", "No interaction", "($T=50\mu$s)", save=True)
plot_phase_space(p1_x, p1_v_x, p2_x, p2_v_x, "Particle 1", "Particle 2", "No interaction", "($T=50\mu$s)", save=True)
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

# plt.figure()
# plt.plot(frequency, frac_f1, label = "$f = 0.1$")
# plt.plot(frequency, frac_f2, label = "$f = 0.4$")
# plt.plot(frequency, frac_f3, label = "$f = 0.7$")
# plt.legend()
# plt.xlabel("Frequency $\omega_V$")
# plt.ylabel("Fraction trapped")
# plt.savefig("figs/fraction_frequency_plot.pdf")

plt.figure()
x = readfile("data/x_values_resonance_interaction_1_dt_0.05_wV_1.5_f_0.7.txt")
y = readfile("data/y_values_resonance_interaction_1_dt_0.05_wV_1.5_f_0.7.txt")
x1 = readfile("data/x_values_resonance_interaction_0_dt_0.05_wV_1.5_f_0.7.txt")
y1 = readfile("data/y_values_resonance_interaction_0_dt_0.05_wV_1.5_f_0.7.txt")
plt.plot(x[:400], y[:400], label = "With interaction")
plt.plot(x1[:450], y1[:450], label = "No interaction")
plt.legend()
plt.show()

plt.figure()
z = readfile("data/z_values_resonance_interaction_1_dt_0.05_wV_1.5_f_0.7.txt")
z1 = readfile("data/z_values_resonance_interaction_0_dt_0.05_wV_1.5_f_0.7.txt")
time = np.arange(0, 500, 0.05)
plt.plot(time[:400], z[:400])
plt.plot(time[:500], z1[:500])
plt.show()

plt.figure()
frac_0 = readfile("data/fraction_interaction_0_dw_0.005_f_0.1_tot_time_500.txt")
frac_1 = readfile("data/fraction_interaction_1_dw_0.005_f_0.1_tot_time_500.txt")
frequency = []
with open("data/fraction_interaction_0_dw_0.005_f_0.1_tot_time_500.txt", 'r') as f:
    for line in f:
        frequency.append(float(line.split()[1]))

plt.plot(frequency, frac_0, label = "No interaction; $f = 0.1$")
plt.plot(frequency, frac_1, label = "With interaction; $f = 0.1$")
plt.legend()
plt.xlabel("Frequency $\omega_V$")
plt.ylabel("Fraction trapped")
plt.show()