import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


def readfile(filename):

    values = []
    with open(filename, 'r') as f:
        for line in f:
            values.append(float(line.split()[0]))

    return np.array(values)

def plot_xy_trajectory(x1, y1, x2, y2, label1, label2, title1, title2, save=False):
    cm = 1/2.54
    plt.figure(figsize=(12*cm, 12*cm))
    plt.plot(x1, y1, label = label1)
    plt.plot(x2, y2, label = label2)
    plt.title(title1 + " " + title2)
    plt.legend()
    plt.xlabel("$x$ values ($\mu$m)")
    plt.ylabel("$y$ values ($\mu$m)")
    plt.tight_layout()
    if save:
        plt.savefig(f"figs/xy_trajectory_{label1}_{label2}_{title1}.pdf")
        plt.close()


def plot_position_time(position, time, label, title, save=False):
    cm = 1/2.54
    plt.figure(figsize=(12*cm, 10*cm))
    plt.plot(time, position, label = label)
    plt.legend()
    plt.xlabel("Time ($\mu$s)")
    plt.ylabel("Position ($\mu$m)")
    plt.title(title)
    plt.tight_layout()
    if save:
        plt.savefig(f"figs/position_time_{label}_{title}.pdf")
        plt.close()


def plot_phase_space(x1, v1, x2, v2, label1, label2, title1, title2, save=False):
    cm = 1/2.54
    plt.figure(figsize=(12*cm, 12*cm))
    plt.plot(x1, v1, label = label1)
    plt.plot(x2, v2, label = label2)
    plt.title(title1 + " " + title2)
    plt.legend()
    plt.xlabel("$x$ values ($\mu$m)")
    plt.ylabel("$v_x$ values ($\mu$m / $\mu$s)")
    plt.tight_layout()
    if save:
        plt.savefig(f"figs/phase_space_{label1}_{label2}_{title1}.pdf")
        plt.close()


def plot_xyz_trajectory(x1, y1, z1, x2, y2, z2, label1, label2, title, save=False):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(x1, y1, z1, label = label1)
    ax.plot3D(x2, y2, z2, label = label2)
    plt.legend()
    plt.title(title)
    plt.tight_layout()
    if save:
        plt.savefig(f"figs/xyz_trajectory_{label1}_{label2}_{title}.pdf")
        plt.close()


def plot_rel_error(errors_list, title, Nsteps_list = [4000, 8000, 16000, 32000], tot_time=50):
    cm = 1/2.54
    plt.figure(figsize=(12*cm, 10*cm))
    
    for i in range(len(errors_list)):
        N_steps = Nsteps_list[i]
        dt = tot_time/N_steps
        time = np.arange(dt, tot_time, dt)
        rel_error = errors_list[i]
        plt.plot(time, rel_error, label = f"$N = {N_steps}$")
    plt.yscale('log')
    plt.ylabel("Relative error")
    #plt.ylabel("$\dfrac{|r_i - r|}{|r|}$")    
    plt.xlabel("Time ($\mu$s)")
    plt.legend()
    plt.title(title)
    plt.tight_layout()    
    plt.savefig(f"figs/rel_error_{title}.pdf")
    plt.close()


def convergence_rate(errors, Nstep=[4000,8000,16000,32000], tot_time=50):
    r = 0
    for i in range(1,4):
        dt = tot_time / Nstep[i]
        dtt = tot_time / Nstep[i-1]
        r += np.log(np.max(errors[i]) / np.max(errors[i-1])) / np.log(dt/dtt)
    
    return 1/3*r        


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

# plot with interaction
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

# plot without interaction
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

# check that time dependence works

# time = np.arange(0, 50, 0.001)

# p1_z_time_dependent_FE = readfile("data/z_values_time_dependent_p1_f_0.7_w_2.2_FE_dt_0.001.txt")
# p1_z_time_dependent_RK4 = readfile("data/z_values_time_dependent_p1_f_0.7_w_2.2_RK4_dt_0.001.txt")
# p1_z_no_time_dependence_FE = readfile("data/z_values_FE.txt")
# p1_z_no_time_dependence_RK4 = readfile("data/z_values_FE.txt")

# plt.figure()
# plt.plot(time, p1_z_no_time_dependence_FE, '--', label = "FE; No time dependence")
# plt.plot(time, p1_z_no_time_dependence_RK4, ':', label = "RK4; No time dependence")

# plt.plot(time, p1_z_time_dependent_FE,label = "FE; Time-dependent E field")
# plt.plot(time, p1_z_time_dependent_RK4, '-.', label = "RK4; Time-dependent E field")
# plt.legend()
# plt.show()

frac_f1 = readfile("data/fraction_f_0.1_500us.txt")
frac_f2 = readfile("data/fraction_f_0.4_500us.txt")
frac_f3 = readfile("data/fraction_f_0.7_500us.txt")
frequency = []
with open("data/fraction_f_0.1_500us.txt", 'r') as f:
    for line in f:
        frequency.append(float(line.split()[1]))

plt.figure()
plt.plot(frequency, frac_f1, label = "$f = 0.1$")
plt.plot(frequency, frac_f2, label = "$f = 0.4$")
plt.plot(frequency, frac_f3, label = "$f = 0.7$")
plt.legend()
plt.xlabel("Frequency $\omega_V$")
plt.ylabel("Fraction trapped")
plt.savefig("figs/fraction_frequency_plot.pdf")

frac_f1 = readfile("data/fraction_f_0.1_500us_dw_0.005.txt")
frac_f2 = readfile("data/fraction_f_0.4_500us_dw_0.005.txt")
frac_f3 = readfile("data/fraction_f_0.7_500us_dw_0.005.txt")
frequency = []
with open("data/fraction_f_0.1_500us_dw_0.005.txt", 'r') as f:
    for line in f:
        frequency.append(float(line.split()[1]))

plt.figure()
plt.plot(frequency, frac_f1, label = "$f = 0.1$")
plt.plot(frequency, frac_f2, label = "$f = 0.4$")
plt.plot(frequency, frac_f3, label = "$f = 0.7$")
plt.legend()
plt.xlabel("Frequency $\omega_V$")
plt.ylabel("Fraction trapped")
plt.savefig("figs/fraction_frequency_plot_close_up.pdf")

