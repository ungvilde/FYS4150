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

def plot_xyz_trajectory(x1, y1, z1, x2, y2, z2, label1, label2, title):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(x1, y1, z1, label = label1)
    ax.plot3D(x2, y2, z2, label = label2)
    plt.legend()
    plt.title(title)
    plt.tight_layout()
    plt.savefig(f"figs/xyz_trajectory__{label1}_{label2}_{title}.pdf")

def plot_rel_error(errors_list, title, Nsteps_list = [4000, 8000, 16000, 32000], tot_time=50):
    cm = 1/2.54
    plt.figure(figsize=(12*cm, 8*cm))
    
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

error_RK4_4 = readfile("data/rel_error_RK4_dt_0.0125.txt")
error_RK4_8 = readfile("data/rel_error_RK4_dt_0.00625.txt")
error_RK4_16 = readfile("data/rel_error_RK4_dt_0.003125.txt")
error_RK4_32 = readfile("data/rel_error_RK4_dt_0.0015625.txt")

errors = [error_RK4_4[1:], error_RK4_8[1:], error_RK4_16[1:], error_RK4_32[1:]]
plot_rel_error(errors, "Runge-Kutta 4th order")

error_FE_4 = readfile("data/rel_error_FE_dt_0.0125.txt")
error_FE_8 = readfile("data/rel_error_FE_dt_0.00625.txt")
error_FE_16 = readfile("data/rel_error_FE_dt_0.003125.txt")
error_FE_32 = readfile("data/rel_error_FE_dt_0.0015625.txt")

errors = [error_FE_4[1:], error_FE_8[1:], error_FE_16[1:], error_FE_32[1:]]
plot_rel_error(errors, "Forward Euler")

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

plot_xy_trajectory(p1_x, p1_y, p2_x, p2_y, "Particle 1", "Particle 2", "With interaction", "($T=50\mu$s)", save=True)
plot_phase_space(p1_x, p1_v_x, p2_x, p2_v_x, "Particle 1", "Particle 2", "With interaction", "($T=50\mu$s)", save=True)
plot_xyz_trajectory(p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, "Particle 1", "Particle 2", "With interaction")

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

plot_xy_trajectory(p1_x, p1_y, p2_x, p2_y, "Particle 1", "Particle 2", "No interaction", "($T=50\mu$s)", save=True)
plot_phase_space(p1_x, p1_v_x, p2_x, p2_v_x, "Particle 1", "Particle 2", "No interaction", "($T=50\mu$s)", save=True)
plot_xyz_trajectory(p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, "Particle 1", "Particle 2", "No interaction")
