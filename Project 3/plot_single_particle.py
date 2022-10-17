from cProfile import label
from importlib.resources import read_text
import numpy as np
import matplotlib.pyplot as plt
from regex import V0

def readfile(filename):

    values = []
    with open(filename, 'r') as f:
        for line in f:
            values.append(float(line.split()[0]))

    return np.array(values)

x_vals_FE = readfile("x_values.txt")
z_vals_FE = readfile("z_values.txt")

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

    return x, y, z, R_pluss, R_minus

x_vals = []
y_vals = []
z_vals = []

dt = 0.001
time = np.arange(0, 50, dt)
print("N_steps = ", len(time))

for t in time:
    x, y, z, _, _ = solve_analytic(t)
    x_vals.append(x)
    y_vals.append(y)
    z_vals.append(z)

_, _, _, R_pluss, R_minus = solve_analytic(0)

plt.plot(time, z_vals,label = "Analytic")
plt.plot(time, z_vals_FE, '--', label = "FE; dt = 0.001")
plt.legend()
plt.xlabel("$\mu s$")
plt.ylabel("$\mu m$")
#plt.hlines(y =[R_minus, R_pluss], xmin = 0, xmax=49)
plt.show()

plt.plot(time, x_vals,label = "Analytic")
plt.plot(time, x_vals_FE, '--', label = "FE; dt = 0.001")
plt.legend()
plt.xlabel("$\mu s$")
plt.ylabel("$\mu m$")
#plt.hlines(y =[R_minus, R_pluss], xmin = 0, xmax=49)
plt.show()

#plt.plot(time, y_vals)
#plt.plot(time, z_vals)