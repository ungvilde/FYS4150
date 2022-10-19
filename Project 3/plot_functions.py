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