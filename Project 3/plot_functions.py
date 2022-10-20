from turtle import color, position
from matplotlib import lines
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
    plt.plot([x1[0],x2[0]], [y1[0], y2[0]], 'ro')
    plt.plot([x1[-1], x2[-1]], [y1[-1], y2[-1]], 'k+')

    #plt.title(title1 + " " + title2)
    plt.legend()
    plt.xlabel("$x$ values ($\mu$m)")
    plt.ylabel("$y$ values ($\mu$m)")
    plt.tight_layout()
    if save:
        plt.savefig(f"figs/xy_trajectory_{label1}_{label2}_{title1}.pdf")
        plt.close()


def plot_position_time(positions, time, labels, ylab, title, save=False, plot_period=False):
    cm = 1/2.54
    plt.figure(figsize=(12*cm, 10*cm))
    n = len(positions)
    for i in range(n):
        plt.plot(time, positions[i], label = labels[i])
        if plot_period and i==0:
            ii = np.argmax(positions[i])
            jj = np.argmax(positions[i][(ii+1):])

            T = time[jj] - time[ii]
            plt.hlines(y=np.max(positions[i]), xmin=time[ii], xmax=time[jj], linestyles='dotted', colors='r')
            plt.plot([time[ii], time[jj]], [np.max(positions[i]), np.max(positions[i])],  "r.")
            plt.text(x = 0, y = np.max(positions[i])+2, s="$T_{num}" + f"= ${T}" + " $\mu$s")
            plt.ylim([-25,25])

    plt.legend(loc=4)
    plt.xlabel("Time ($\mu$s)")
    plt.ylabel(ylab + " ($\mu$m)")
    #plt.title(title)
    plt.tight_layout()
    if save:
        plt.savefig(f"figs/position_time_{ylab}_{title}.pdf")
        plt.close()


def plot_phase_space(x1, v1, x2, v2, label1, label2, title1, title2, xlab, ylab, save=False):
    cm = 1/2.54
    plt.figure(figsize=(12*cm, 12*cm))
    plt.plot(x1, v1, label = label1)
    plt.plot(x2, v2, label = label2)
    plt.plot([x1[0],x2[0]], [v1[0], v2[0]], 'ro')
    plt.plot([x1[-1], x2[-1]], [v1[-1], v2[-1]], 'k+')
    #plt.title(title1 + " " + title2)
    plt.legend(loc=4)
    plt.xlabel(xlab + " ($\mu$m)")
    plt.ylabel(ylab + "  (m/s)")
    plt.tight_layout()
    if save:
        plt.savefig(f"figs/phase_space_{label1}_{label2}_{title1}+{xlab}.pdf")
        plt.close()


def plot_xyz_trajectory(x1, y1, z1, x2, y2, z2, label1, label2, title, save=False):
    cm = 1/2.54
    #ax = plt.figure(figsize=(12*cm, 12*cm))
    #ax = plt.axes(projection='3d')
    ax = plt.figure(figsize=(12*cm, 12*cm)).add_subplot(projection='3d')
    ax.plot3D(x1, y1, z1, label = label1)
    ax.plot3D(x2, y2, z2, label = label2)
    ax.scatter3D([x1[0], x2[0]], [y1[0], y2[0]], [z1[0], z2[0]], marker = 'o', color='r', depthshade=0)
    ax.scatter3D([x1[-1], x2[-1]], [y1[-1], y2[-1]], [z1[-1], z2[-1]], marker = '+', color='k', depthshade=0)
    ax.legend(loc=4)
    #plt.title(title)
    plt.tight_layout()
    if save:
        #plt.show()
        plt.savefig(f"figs/xyz_trajectory_{label1}_{label2}_{title}.pdf")
        plt.close()


def plot_rel_error(errors_list, title, text, Nsteps_list = [4000, 8000, 16000, 32000], tot_time=50):
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
    plt.text(x=0.1, y=0.1, s=text, transform=plt.gca().transAxes)
    #plt.title(title)
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

def plot_frequency_fraction(frequency, fraction_list, labels_list, title):
    n = len(fraction_list)

    cm = 1/2.54
    plt.figure(figsize=(12*cm, 10*cm))
    for i in range(n):
        plt.plot(frequency, fraction_list[i], label = labels_list[i])
    plt.xlabel("Frequency $\omega_V$ (MHz)")
    plt.ylabel("Fraction")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"figs/frequency_fraction_{title}.pdf")
    plt.close()