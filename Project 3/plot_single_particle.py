import numpy as np
import matplotlib.pyplot as plt

def readfile(filename):

    values = []
    with open(filename, 'r') as f:
        for line in f:
            values.append(float(line.split()[0]))

    return np.array(values)

x_vals_FE = readfile("data/x_values_FE.txt")
y_vals_FE = readfile("data/y_values_FE.txt")
z_vals_FE = readfile("data/z_values_FE.txt")

x_vals_RK4 = readfile("data/x_values_RK4.txt")
y_vals_RK4 = readfile("data/y_values_RK4.txt")
z_vals_RK4 = readfile("data/z_values_RK4.txt")

def solve_analytic(t):
    v0 = 25 # init. velocity
    B0 = 9.65 * 10
    V0 = 2.41*10**6 # init. potential
    m = 40.078 # mass
    d = 500 # char. distance
    q = 1 # charge
    w_z = np.sqrt(2*q*V0 / (m * d**2))
    w0 = q*B0 / m

    # initial conditions
    z0 = 20
    y0 = 0
    x0 = 20
    dx0 = 0
    dy0 = v0
    dz0 = 0

    z = z0 * np.cos(w_z * t)

    w_pluss = (w0 + np.sqrt(w0**2 - 2*w_z**2)) / 2
    w_minus = (w0 - np.sqrt(w0**2 - 2*w_z**2)) / 2

    A_pluss = (v0 + w_minus*x0)/(w_minus - w_pluss)
    A_minus = -(v0 + w_pluss*x0)/(w_minus - w_pluss)

    R_pluss = np.abs(A_pluss + A_minus)
    R_minus = np.abs(A_pluss - A_minus)

    x = A_pluss*np.cos(-w_pluss*t) + A_minus*np.cos(-w_minus*t)
    y = A_pluss*np.sin(-w_pluss*t) + A_minus*np.sin(-w_minus*t)

    return x, y, z

x_vals = []
y_vals = []
z_vals = []

dt = 0.001
time = np.arange(0, 50, dt)
for t in time:
    x, y, z = solve_analytic(t)
    x_vals.append(x)
    y_vals.append(y)
    z_vals.append(z)

cm = 1/2.54

plt.figure(figsize=(12*cm, 8*cm))
plt.plot(time, z_vals, 'k', label = "Analytic")
plt.plot(time, z_vals_RK4, 'r--', label = "RK4; dt = 0.001")
plt.plot(time, z_vals_FE, 'b:', label = "FE; dt = 0.001")
plt.legend()
plt.xlabel("$\mu s$")
plt.ylabel("$\mu m$")
plt.title("$z$ values")
plt.tight_layout()
plt.savefig("figs/zvals_FE_analytic_single_particle.pdf")

plt.figure(figsize=(12*cm, 8*cm))
plt.plot(time, x_vals, 'k',label = "Analytic")
plt.plot(time, x_vals_RK4, 'r--', label = "RK4; dt = 0.001")
plt.plot(time, x_vals_FE, 'b:', label = "FE; dt = 0.001")
plt.legend()
plt.xlabel("$\mu s$")
plt.ylabel("$\mu m$")
plt.title("$x$ values")
plt.tight_layout()
plt.savefig("figs/xvals_FE_analytic_single_particle.pdf")

plt.figure(figsize=(12*cm, 8*cm))
plt.plot(time, y_vals, 'k',label = "Analytic")
plt.plot(time, y_vals_RK4, 'r--', label = "RK4; dt = 0.001")
plt.plot(time, y_vals_FE, 'b:', label = "FE; dt = 0.001")
plt.legend()
plt.xlabel("$\mu$s")
plt.ylabel("$\mu$m")
plt.title("$y$ values")
plt.tight_layout()
plt.savefig("figs/yvals_FE_analytic_single_particle.pdf")
